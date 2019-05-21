#!/usr/bin/python3

"""
Part of the SqueezeMeta distribution. 25/03/2018 Original version,
                            (c) Fernando Puente-Sánchez, CNB-CSIC.

Generate tabular outputs from SqueezeMeta results.

USAGE: make-tables.py PROJECT_NAME OUTPUT_DIRECTORY

OPTIONS:
    --trusted_functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables
    --ignore_unclassified: Ignore ORFs without assigned functions in
        TPM calculation
    --write_parsed_tax: Write tables with the per-rank taxonomy of each
        ORF and contig
"""

from os.path import abspath, dirname, realpath
from os import mkdir
import argparse

from collections import defaultdict

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.append('{}/../lib/'.format(utils_home))
from utils import parse_conf_file, parse_orf_table, parse_tax_table, parse_contig_table, aggregate_tax_abunds, normalize_abunds, TAXRANKS, TAXRANKS_SHORT 


def main(args):
    ### Get result files paths from SqueezeMeta_conf.pl
    perlVars = parse_conf_file(args.project_path)
    nokegg, nocog, nopfam, doublepass = map(int, [perlVars['$nokegg'], perlVars['$nocog'], perlVars['$nopfam'], perlVars['$doublepass']])

    ### Create output dir.
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise

    ### Calculate tables and write results.
    prefix = args.output_dir + '/' + perlVars['$projectname'] + '.'
    
    sampleNames, orfs, kegg, cog, pfam = parse_orf_table(perlVars['$mergedfile'], nokegg, nocog, nopfam, args.trusted_functions, args.ignore_unclassified)

    # Round aggregated functional abundances.
    # We can have non-integer abundances bc of the way we split counts in ORFs with multiple KEGGs.
    # We round the aggregates to the closest integer for convenience.
    kegg['abundances'] = {k: a.round().astype(int) for k,a in kegg['abundances'].items()}
    cog['abundances']  = {k: a.round().astype(int) for k,a in  cog['abundances'].items()}
    pfam['abundances'] = {k: a.round().astype(int) for k,a in pfam['abundances'].items()}

    #write_results(sampleNames, orfs['tpm'], prefix + 'orf.tpm.tsv')
    if not nokegg:
        write_results(sampleNames, kegg['abundances'], prefix + 'KO.abund.tsv')
        write_results(sampleNames, kegg['tpm'], prefix + 'KO.tpm.tsv')
    if not nocog:
        write_results(sampleNames, cog['abundances'], prefix + 'COG.abund.tsv')
        write_results(sampleNames, cog['tpm'], prefix + 'COG.tpm.tsv')
    if not nopfam:
        write_results(sampleNames, pfam['abundances'], prefix + 'PFAM.abund.tsv')
        write_results(sampleNames, pfam['tpm'], prefix + 'PFAM.tpm.tsv')

    fun_prefix = perlVars['$fun3tax_blastx'] if doublepass else perlVars['$fun3tax']
    orf_tax, orf_tax_wranks = parse_tax_table(fun_prefix + '.wranks')
    orf_tax_nofilter, orf_tax_nofilter_wranks = parse_tax_table(fun_prefix + '.noidfilter.wranks')

    # Add ORFs not present in the input tax file.
    unclass_list = ['Unclassified' for rank in TAXRANKS_SHORT]
    unclass_list_wranks = ['{}_Unclassified' for rank in TAXRANKS_SHORT]
    for orf in orfs['abundances']:
        if orf not in orf_tax:
            assert orf not in orf_tax_wranks
            assert orf not in orf_tax_nofilter
            assert orf not in orf_tax_nofilter_wranks
            orf_tax[orf] = unclass_list
            orf_tax_wranks[orf] = unclass_list_wranks
            orf_tax_nofilter[orf] = unclass_list
            orf_tax_nofilter_wranks[orf] = unclass_list_wranks

    
    orf_tax_prokfilter, orf_tax_prokfilter_wranks = {}, {}

    for orf in orf_tax:
        tax = orf_tax[orf]
        tax_nofilter = orf_tax_nofilter[orf]
        if 'Bacteria' in (tax[0], tax_nofilter[0]) or 'Archaea' in (tax[0], tax_nofilter[0]): # We check both taxonomies.
            orf_tax_prokfilter       [orf] = tax
            orf_tax_prokfilter_wranks[orf] = orf_tax_wranks[orf]
        else:
            orf_tax_prokfilter       [orf] = tax_nofilter
            orf_tax_prokfilter_wranks[orf] = orf_tax_nofilter_wranks[orf]
    
    contig_abunds, contig_tax, contig_tax_wranks = parse_contig_table(perlVars['$contigtable'])

    if args.write_parsed_tax:
        write_results(TAXRANKS, orf_tax, prefix + 'orf.tax.allfilter.tsv')
        write_results(TAXRANKS, orf_tax_nofilter, prefix + 'orf.tax.nofilter.tsv')
        write_results(TAXRANKS, orf_tax_prokfilter, prefix + 'orf.tax.prokfilter.tsv')
        write_results(TAXRANKS, contig_tax, prefix + 'contig.tax.tsv')

    for idx, rank in enumerate(TAXRANKS):
        tax_abunds_orfs = aggregate_tax_abunds(orfs['abundances'], orf_tax_prokfilter, idx)
        write_results(sampleNames, tax_abunds_orfs, prefix + '{}.prokfilter.abund.orfs.tsv'.format(rank))

        tax_abunds_contigs = aggregate_tax_abunds(contig_abunds, contig_tax_wranks, idx)
        write_results(sampleNames, tax_abunds_contigs, prefix + '{}.allfilter.abund.tsv'.format(rank))
        write_results(sampleNames, normalize_abunds(tax_abunds_contigs, 100), prefix + '{}.allfilter.percent.tsv'.format(rank))



def write_results(sampleNames, rowDict, outname):
    with open(outname, 'w') as outfile:
        outfile.write('\t{}\n'.format('\t'.join(sampleNames)))
        for row in sorted(rowDict):
            outfile.write('{}\t{}\n'.format(row, '\t'.join(map(str, rowDict[row]))))



def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate SqueezeMeta results into tables', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('--trusted_functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--ignore_unclassified', action='store_true', help='Ignore ORFs without assigned functions in TPM calculation')
    parser.add_argument('--write_parsed_tax', action='store_true', help='Write tables with the per-rank taxonomy of each ORF and contig')

    return parser.parse_args()




if __name__ == '__main__':
    main(parse_args())
