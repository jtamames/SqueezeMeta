#!/usr/bin/env python
"""
@author: inodb
"""
import sys
import os
import argparse
import subprocess
import errno
from signal import signal, SIGPIPE, SIG_DFL

from Bio import SeqIO
from Bio.SeqUtils import GC

from BCBio import GFF

TAXONOMY = ('phylum', 'class', 'order', 'family', 'genus', 'species')

def get_gc_and_len_dict(fastafile):
    """Creates a dictionary with the fasta id as key and GC and length as keys
    for the inner dictionary."""
    out_dict = {}

    for rec in SeqIO.parse(fastafile, "fasta"):
        out_dict[rec.id] = {}
        out_dict[rec.id]["GC"] = GC(rec.seq)
        out_dict[rec.id]["length"] = len(rec.seq)

    return out_dict


def get_bedcov_dict(bedcoverage):
    """Uses the BEDTools genomeCoverageBed histogram output to determine mean
    coverage and percentage covered for each contig.
    
    Returns dict with fasta id as key and percentage covered and cov_mean as
    keys for the inner dictionary."""
    out_dict = {}

    # Check if given argument is a file, otherwise use the content of the
    # variable
    if os.path.isfile(bedcoverage):
        fh = open(bedcoverage)
    else:
        fh = bedcoverage.split('\n')[:-1]

    for line in fh:
        cols = line.split()

        try:
            d = out_dict[cols[0]]
        except KeyError:
            d = {}
            out_dict[cols[0]] = d

        if int(cols[1]) == 0:
            d["percentage_covered"] = 100 - float(cols[4]) * 100.0
        else:
            d["cov_mean"] = d.get("cov_mean", 0) + int(cols[1]) * float(cols[4])

    return out_dict

def get_gff_dict(gfffile):
    """Creates a dictionary with product information from given gff file.
    
    Returns dictionary. Dictionary key is the contig id, values are products for the contig."""
    out_dict = {}

    for rec in GFF.parse(gfffile):

        # Add features if there are any
        if rec.features > 0:
            gff_info = None
                
            # Add all features
            # Features are separated by ,
            # example:
            # featuretype;product;product,featuretype;product
            # or
            # CDS;protein3;protein31,CDS;protein3
            for f in rec.features:
                if len(f.qualifiers['product']) > 0:
                    # if gff_info is None, do not add ',' separator
                    try:
                        gff_info += ",%s" % ";".join([f.type] + f.qualifiers['product'])
                    except TypeError:
                        gff_info = ";".join([f.type] + f.qualifiers['product'])

            # Test if there were any features with a product
            if gff_info == None:
                gff_info = "N/A"
        else:
            gff_info = "N/A"

        out_dict[rec.id] = gff_info

    return out_dict

def print_sample_columns(t):
    sys.stdout.write(("\tcov_mean_sample_%s" * len(t)) % t)
    sys.stdout.write(("\tpercentage_covered_%s" * len(t)) % t)

def print_input_table(fastadict, bedcovdicts, taxonomydict=None, gffdict=None, samplenames=None):
    """Writes the input table for Probin to stdout. See hackathon google
    docs."""

    # Header
    sys.stdout.write(("%s" + "\t%s" * 9) % (('contig', 'length', 'GC') + TAXONOMY + ('gff_info',)))
    if samplenames == None:
        # Use index if no sample names given in header
        print_sample_columns(tuple(range(len(bedcovdicts))))
    else:
        # Use given sample names in header
        assert(len(samplenames) == len(bedcovdicts))
        print_sample_columns(tuple(samplenames))
    sys.stdout.write("\n")

    # Content
    assert(len(fastadict) > 0)
    for acc in fastadict:
        # fasta stats
        sys.stdout.write("%s\t%s\t%s"  %
            (
                acc,
                fastadict[acc]['length'],
                fastadict[acc]['GC']
            )
        )

        # taxonomy
        if taxonomydict != None:
            for t in TAXONOMY:
                try:
                    sys.stdout.write("\t%s" % taxonomydict[acc][t])
                except KeyError:
                    sys.stdout.write("\tN/A")
        else:
            sys.stdout.write("\tN/A" * len(TAXONOMY))

        # gff data
        if gffdict:
            try:
                sys.stdout.write("\t%s" % gffdict[acc])
            except KeyError:
                sys.stdout.write("\tN/A")
        else:
            sys.stdout.write("\tN/A")

        # Print mean
        for bcd in bedcovdicts:
            try:
                # Print cov mean
                sys.stdout.write("\t%f" % (bcd[acc]["cov_mean"]))
            except KeyError:
                # No reads mapped to this contig
                sys.stdout.write("\t0")
                bcd[acc] = {"percentage_covered":0}
        

        # Print percentage covered
        for bcd in bedcovdicts:
            # If no 0 bases with 0 coverage, then all bases in the contig are covered
            sys.stdout.write("\t%f" % bcd[acc].get("percentage_covered", 100))

        sys.stdout.write("\n")


def get_taxonomy_dict(taxonomyfile):
    """Creates a dictionary from given taxonomy file.
    
    The dictionary has contig name as key for the outer dictionary and genus
    and species as keys for the inner dictionary."""
    outdict = {}

    for line in open(taxonomyfile):
        cols = line.split(',')

        # Should be 7 columns. Contig name, Phylum, Class, Order, Family,
        # Genus, Species.
        assert(len(cols) == len(TAXONOMY) + 1)

        outdict[cols[0]] = dict(list(zip(TAXONOMY, cols[1:-1] + [cols[-1].rstrip('\n')])))

    return outdict


def generate_input_table(fastafile, bamfiles, taxonomyfile=None, gfffile=None,
    samplenames=None, isbedfiles=False):
    """Reads input files into dictionaries then prints everything in the table
    format required for running ProBin."""
    bedcovdicts = []
    
    # Determine coverage information from bam file using BEDTools
    for i, bf in enumerate(bamfiles):
        if isbedfiles == False:
            p = subprocess.Popen(["genomeCoverageBed", "-ibam", bf], stdout=subprocess.PIPE)
            out, err = p.communicate()
            if p.returncode != 0:
                sys.stderr.write(out)
                raise Exception('Error with genomeCoverageBed')
            else:
                bedcovdicts.append(get_bedcov_dict(out))
        else:
            bedcovdicts.append(get_bedcov_dict(bf))
                

    # Determine annotations from gfffile
    gffdict = get_gff_dict(gfffile) if gfffile != None else None

    # Determine taxonomies from taxonomyfile
    taxonomydict = get_taxonomy_dict(taxonomyfile) if taxonomyfile != None else None

    print_input_table(get_gc_and_len_dict(fastafile), bedcovdicts,
        taxonomydict=taxonomydict, gffdict=gffdict, samplenames=samplenames)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="Contigs fasta file")
    parser.add_argument("bamfiles", nargs='+', help="BAM files with mappings to contigs")
    parser.add_argument("--taxonomyfile", help="The taxonomy for the contigs. Has"
    " 7 columns. Contig name, Phylum, Class, Order, Family, Genus, Species.")
    parser.add_argument("--gfffile", default=None, help="GFF file with features info.")
    parser.add_argument("--samplenames", default=None, help="File with sample names, one line each. Should be same nr as bamfiles.")
    parser.add_argument("--isbedfiles", action='store_true',
        help="The bamfiles argument are outputs of genomeCoverageBed, not the actual bam file. Skips running genomeCoverageBed from within this script.")
    args = parser.parse_args()

    # Get sample names
    if args.samplenames != None:
        samplenames = [ s[:-1] for s in open(args.samplenames).readlines() ]
        if len(samplenames) != len(args.bamfiles):
            raise Exception("Nr of names in samplenames should be equal to nr of given bamfiles")
    else:
        samplenames=None
    
    # ignore broken pipe error when piping output
    # http://newbebweb.blogspot.pt/2012/02/python-head-ioerror-errno-32-broken.html
    signal(SIGPIPE,SIG_DFL)

    generate_input_table(args.fastafile, args.bamfiles,
        taxonomyfile=args.taxonomyfile, gfffile=args.gfffile,
        samplenames=samplenames, isbedfiles=args.isbedfiles)
