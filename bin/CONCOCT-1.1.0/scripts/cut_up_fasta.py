#!/usr/bin/env python
"""Cut up fasta file in non-overlapping or overlapping parts of equal length.

Optionally creates a BED-file where the cutup contigs are specified in terms
of the original contigs. This can be used as input to concoct_coverage_table.py.
"""
from __future__ import print_function
import argparse
from Bio import SeqIO


def cut_up_fasta(fastfiles, chunk_size, overlap, merge_last, bedoutfile):
    if bedoutfile:
        bedoutfile_fh = open(bedoutfile, 'w')
    for ff in fastfiles:
        for record in SeqIO.parse(ff, "fasta"):
            if (not merge_last and len(record.seq) > chunk_size) or (merge_last and len(record.seq) >= 2 * chunk_size):
                i = 0
                for split_seq in chunks(record.seq, chunk_size, overlap, merge_last):
                    print(">{}.concoct_part_{}\n{}".format(record.id, i, split_seq))
                    if bedoutfile:
                        print("{0}\t{2}\t{3}\t{0}.concoct_part_{1}".format(record.id, i, chunk_size*i, chunk_size*i+len(split_seq)),
                              file=bedoutfile_fh)
                    i = i + 1
            else:
                print(">{}.concoct_part_0\n{}".format(record.id, record.seq))
                if bedoutfile:
                    print("{0}\t0\t{1}\t{0}.concoct_part_0".format(record.id, len(record.seq)),
                          file=bedoutfile_fh)

    if bedoutfile:
        bedoutfile_fh.close()

def chunks(l, n, o, merge_last):
    """ Yield successive n-sized chunks from l with given overlap o between the
    chunks.
    """
    assert n > o

    if not merge_last:
        for i in range(0, len(l), n - o):
            yield l[i:i + n]
    else:
        for i in range(0, len(l) - n + 1, n - o):
            yield l[i:i + n] if i + n + n - o <= len(l) else l[i:]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "contigs", nargs="+", help="Fasta files with contigs\n")
    parser.add_argument("-c", "--chunk_size", default=1999, type=int, help="Chunk size\n")
    parser.add_argument("-o", "--overlap_size", default=1900, type=int, help="Overlap size\n")
    parser.add_argument("-m", "--merge_last", default=False, action="store_true", help="Concatenate final part to last contig\n")
    parser.add_argument("-b", "--bedfile", default=None, help="BEDfile to be created with exact regions of the original"
                                                              " contigs corresponding to the newly created contigs")
    args = parser.parse_args()
    cut_up_fasta(args.contigs, args.chunk_size, args.overlap_size, args.merge_last, args.bedfile)
