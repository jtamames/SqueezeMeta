from collections import defaultdict, namedtuple
from math import ceil
from statistics import mean, stdev, quantiles
from datetime import datetime
from sys import stdout, stderr
from dataclasses import dataclass


from tqdm import tqdm
from pysam import AlignmentFile
import plotext as plt

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

@dataclass(slots=True)
class Pair: # close enough to a mutable namedtuple
    """To store pairs mapping to two different references"""
    query: str
    edge: tuple[str, str]
    ref: str
    rstart: int
    rlengthFull: int
    riden: float
    rqcov: float
    nextref: str
    nrstart: int
    nrlengthFull: int
    nriden: float
    nrqcov: float

#@profile
def main(args):
    
    # To calculate the insert size distribution
    insert_sizes = []

    # To keep track of the observed read pairs
    linker_pairs = {}

    eprint()
    eprint(f'Parsing BAM file {args.bamfile}...')

    added_isize = 0
    missing_nm_tag = 0

    # Iterate through the bam file
    with AlignmentFile(args.bamfile, 'rb') as samfile:

        ref_lengths = dict(zip(samfile.references, samfile.lengths))
        mapped_reads = samfile.mapped

        for record in tqdm(samfile.fetch(), total = mapped_reads): # AlignmentFile.fetch() is already iterating through the mapped reads only

            if not record.has_tag('NM'):
                # Some reads had no NM tag or CIGAR string even though they seemed to have reference name and reference start
                missing_nm_tag += 1
                continue

            query, ref, nextref = record.query_name, record.reference_name, record.next_reference_name
            rstart, nrstart = record.reference_start, record.next_reference_start
            
            iden = round(100 * ( 1 - (record.get_tag('NM') / record.query_alignment_length) ), 2) # NM (edit distance to the template) / aligned length
            qcov = round(100 * record.query_alignment_length / record.query_length, 2) # Aligned length / read length

            # Skip if the mate is unmapped
            if record.mate_is_unmapped:
                continue

            # If this pair is mapping to the same reference, store the insert size
            if ref == nextref:
                if added_isize < args.insert_size_reads:
                    insert_sizes.append( int(abs(rstart - nrstart)) )
                    added_isize += 1
        
            # Otherwise check if we haven't seen this pair before, and record it if that's the case
            elif query not in linker_pairs:

                rlengthFull, nrlengthFull  = ref_lengths[ref], ref_lengths[nextref]

                edge = tuple( sorted( (ref, nextref) ) )

                pair = Pair(query = query, edge = edge,
                            ref = ref, rstart = rstart, rlengthFull = rlengthFull, riden = iden, rqcov = qcov,
                            nextref = nextref, nrstart = nrstart, nrlengthFull = nrlengthFull, nriden = 0, nrqcov = 0)

                linker_pairs[query] = pair

            # If we have seen this pair before, record the reverse identity and query coverage
            else:
                linker_pairs[query].nriden = iden
                linker_pairs[query].nrqcov = qcov
                
    eprint()

    # Get insert size distribution and nth percentile
    eprint('Calculating insert size...')
    mean_insert_size = round(mean(insert_sizes))
    std_insert_size = round(stdev(insert_sizes))
    min_insert_size = min(insert_sizes) 
    nth_percentile_insert_size = quantiles(insert_sizes, n=100)[args.insert_size_percentile-1]

    # Make a plot in the terminal
    # plotext wants two iterables for the barplot, one with the bar names and one with the bar counts
    # So we need to manually aggregate our data into bins
    NBINS = 80
 
    binsize = ceil( (nth_percentile_insert_size - min_insert_size + 1) / NBINS) # make sure it's at least 1
    barplot_isize_ranges = [mean([min_insert_size + binsize*i, min_insert_size + binsize*(i+1)]) for i in range(NBINS)]
    barplot_isize_counts = [0]*NBINS
    for isize in insert_sizes:
        if isize > nth_percentile_insert_size:
            continue
        bin_ = (isize - min_insert_size) // binsize
        barplot_isize_counts[bin_] += 1
    eprint()
    plt.simple_bar(barplot_isize_ranges, barplot_isize_counts, title = 'Insert size distribution', width = plt.terminal_width() - 20)
    #plt.show()
    # Normally we'd use plt.show() here, but we hack it to be able to print to stderr
    eprint(plt.active().monitor.matrix.canvas, add_timestamp = False)

    eprint()
    eprint(f'Insert size Min: {min(insert_sizes)} | {args.insert_size_percentile}th percentile: {nth_percentile_insert_size} | Max: {max(insert_sizes)} | Mean: {mean_insert_size} | Std: {std_insert_size}')
    eprint()

    # Go through all the pairs mapping to two different references
    # See if the mapping positions are consistent with the inferred alignment length
    # So they shouldn't map very far away from the beginning/end of our contigs
    eprint('Building scaffolds...')
    nodes = set()
    edges = defaultdict(int)

    for pair in linker_pairs.values():

        fwd_size = min(pair.rstart , pair.rlengthFull  - pair.rstart )
        rev_size = min(pair.nrstart, pair.nrlengthFull - pair.nrstart)
        insert_size = fwd_size + rev_size

        min_iden = min(pair.riden, pair.nriden)
        min_qcov = min(pair.rqcov, pair.nrqcov)

        if insert_size <= nth_percentile_insert_size and min_iden >= args.min_identity and min_qcov >= args.min_query_coverage:
            edges[pair.edge] += 1
            nodes.update(pair.edge)

    eprint()

    # Write the output
    try:
        outfile = None
        if not args.output:
            outfile = stdout
            oname = 'stdout'
        else:
            outfile = open(args.output, 'w')
            oname = args.output
        eprint(f'Writing output to {oname}...')
        for edge, cov in sorted(edges.items(), key = lambda x: x[1], reverse = True):
            print(f'{edge[0]}\t{edge[1]}\t{cov}', file = outfile)
        eprint()
    finally:
        if args.output and outfile:
            outfile.close()

    eprint('Done')
    eprint(f'Mapped reads in BAM file: {mapped_reads}')
    eprint(f'Reads with no NM tag: {missing_nm_tag}')
    eprint(f'Pairs used for insert size calculation: {added_isize}')
    eprint(f'All scaffolding pairs: {len(linker_pairs)}')
    eprint(f'Good scaffolding pairs: {sum(edges.values())}')
    eprint(f'Contigs in BAM file: {len(ref_lengths)}')
    eprint(f'Nodes in scaffold graph: {len(nodes)}')
    eprint(f'Edges in scaffold graph: {len(edges)}')
    eprint()
        

def eprint(*args, add_timestamp = True, **kwargs):
    if args and add_timestamp:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        args = list(args)
        args[0] = f'{ts}\t{args[0]}'
    print(*args, file=stderr, **kwargs)

def parse_args():
    parser = ArgumentParser(description = 'Use paired reads in a BAM file to scaffold contigs, outputs scaffold edges and weights',
                            formatter_class = ArgumentDefaultsHelpFormatter)
    parser.add_argument('bamfile', type = str,
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', type = str, default = None,
                        help = 'Output file (`None` will redirect to stdout)')
    parser.add_argument('-i', '--min-identity', type = float, default = 95,
                        help = 'Ignore pairs in which at least one mate has a percent alignment identity smaller than this threshold')
    parser.add_argument('-c', '--min-query-coverage', type = float, default = 80,
                        help = 'Ignore pairs in which at least one mate has a percent query coverage smaller than this threshold')
    parser.add_argument('-r', '--insert-size-reads', type = int, default = 100000,
                        help = 'Reads used to calculate insert size distribution')
    parser.add_argument('-p', '--insert-size-percentile', type = int, default = 95,
                        help = 'Ignore pairs with an insert size above the nth percentile of the insert size distribution in the input BAM file')
    return parser.parse_args()

if __name__ == '__main__':
    main(parse_args())
