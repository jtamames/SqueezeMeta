#!/usr/bin/env python
"""
With contigs cutup with cut_up_fasta.py as input, sees to that the consequtive
parts of the original contigs are merged.

prints result to stdout.

@author: alneberg
"""
from __future__ import print_function
import sys
import os
import argparse
from collections import defaultdict, Counter
import re

CONTIG_PART_EXPR = re.compile("(.*)\.concoct_part_([0-9]*)")

def original_contig_name_special(s):
    n = s.split(".")[-1]
    try:
        original_id, part_nr = CONTIG_PART_EXPR.match(s).group(1,2)
        return original_id, part_nr
    except AttributeError: # No matches for concoct_part regex
        return s, 0

def main(args):
    all_seqs = {}
    all_originals = defaultdict(dict)
    first = True
    with open(args.cutup_clustering_result, 'r') as ifh:
        for line in ifh:
            if first:
                first=False
                if 'contig_id' not in line:
                    sys.stderr.write(("ERROR! The term 'contig_id' was not found on the first row. Please make sure that there "
                    "is a header line before continuing. Exiting\n"))
                    sys.exit(-1)
                continue
            line = line.strip()
            contig_id, cluster_id = line.split(',')
            original_contig_name, part_id = original_contig_name_special(contig_id)

            all_originals[original_contig_name][part_id] = cluster_id

    merged_contigs_stack = []

    sys.stdout.write("contig_id,cluster_id\n")
    for original_contig_id, part_ids_d in all_originals.items():
        if len(part_ids_d) > 1:
            c = Counter(part_ids_d.values())
            cluster_id = c.most_common(1)[0][0]
            c_string = [(a,b) for a, b in c.items()]
            if len(c.values()) > 1:
                sys.stderr.write("No consensus cluster for contig {}: {}\t Chosen cluster: {}\n".format(original_contig_id, c_string, cluster_id))
        else:
            cluster_id = list(part_ids_d.values())[0]

        sys.stdout.write("{},{}\n".format(original_contig_id, cluster_id))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cutup_clustering_result", help=("Input cutup clustering result."))
    args = parser.parse_args()

    main(args)
