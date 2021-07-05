#!/usr/bin/env python
"""
@author: inodb
"""
import os
import argparse

import pysam
import multiprocessing

UNAMB_PREFIX = "unamb_read_count_"
AMB_PREFIX = "amb_read_count_"

def get_fasta_accs(fastafile):
    out_list = []

    for line in open(fastafile):
        if line.startswith('>'):
            out_list.append(line.split()[0][1:])

    return out_list


def init_count_dict(contigfa, reffa):
    caccs = get_fasta_accs(contigfa)
    raccs = get_fasta_accs(reffa) 

    column_header = [UNAMB_PREFIX + r for r in raccs] + \
                 [AMB_PREFIX + r for r in raccs]

    # create table like dict, key is contig, ambiguous and unambiguous counters
    # for each genome are keys from inner dict
    count_dict = dict((c, dict((ch, 0) for ch in column_header)) for c in caccs) 

    return count_dict, column_header


def is_ambiguous_align(tags, multi_align_tag):
    """Returns whether the read aligns to multiple locations. The
    multi_align_tag depends on mapper. For bowtie2 it is XS."""
    for t in tags:
        if t[0] == multi_align_tag:
            return True
    return False

def extract_read_ref_origin(readname):
    """Extracts the name of the reference the read originated from.
    
    Read names are expected to look like this:

    >gi|29611500|ref|NC_004703.1|_Bacteroides_thetaiotaomicron_VPI-5482_nr0_+_R1
    
    The first part of the read contains the sequence name (as found in the
    database) from which it originates. The second part identifies the read
    number (e.g. nr0). The + or - identifies the strand from which it was
    sequenced and the last bit tells if it was an R1 or R2 read."""
    return "_".join(readname.split("_")[:-3])
    

def count_contigs_per_genome(bamfile, count_dict, multi_align_tag='XS'):
    bamh = pysam.Samfile(bamfile)

    for read in bamh.fetch():
        # get the contig row this read is mapping to
        # all contig keys should already be in count_dict
        contigrow = count_dict[bamh.getrname(read.tid)]

        # get reference from read name
        ref_origin = extract_read_ref_origin(read.qname)

        # all columns should already be in count_dict initialized with 0 values
        # so just add the count depending on whether the read alignes once or
        # multiple times
        if is_ambiguous_align(read.tags, multi_align_tag):
            contigrow[AMB_PREFIX + ref_origin] += 1
        else:
            contigrow[UNAMB_PREFIX + ref_origin] += 1

    return count_dict
        
    
def print_count_dict(count_dict, column_header):
    print(("contig" + "\t%s" * len(column_header)) % tuple(column_header))

    for contig in count_dict:
        print(("%s" + "\t%i" * len(column_header)) % tuple([contig] + \
            [count_dict[contig][ch] for ch in column_header]))

 
def parallel_count_contigs_per_genome(args):
    bamfile, count_dict = args
    return count_contigs_per_genome(bamfile, count_dict)

def sum_count_dicts(cd1, cd2, column_header):
    assert(len(cd1) > 0 and len(cd2) > 0)

    for c in cd1:
        for ch in column_header:
            cd1[c][ch] += cd2[c][ch]

    return cd1


def main(contigfa, reffa, bamfiles, max_n_processors):
    count_dict, column_header = init_count_dict(contigfa, reffa)

    # Determine counts from bamfiles in parallel
    #
    # NOTE: no need to copy count_dict, multiprocessing does that by itself
    count_args = [(bf, count_dict) for bf in bamfiles]
    n_processes = min(multiprocessing.cpu_count(), max_n_processors)
    pool = multiprocessing.Pool(processes=n_processes)
    poolrv = pool.map(parallel_count_contigs_per_genome, count_args)

    # Sum results
    for rv in poolrv:
        count_dict = sum_count_dicts(count_dict, rv, column_header)

    # print results
    print_count_dict(count_dict, column_header)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("contigfa", help="Contigs fasta file")
    parser.add_argument("reffa", help="Reference fasta file")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser.add_argument('-m','--max_n_processors',type=int, default=1,
        help="Specify the maximum number of processors to use, if absent, all present processors will be used.") 

    args = parser.parse_args()

    for bf in args.bamfiles:
        if not os.path.isfile(bf + ".bai"):
            raise Exception

    main(args.contigfa, args.reffa, args.bamfiles, args.max_n_processors)
