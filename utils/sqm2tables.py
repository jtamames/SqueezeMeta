#!/usr/bin/env python3

"""
Part of the SqueezeMeta distribution. 25/03/2018 Original version,
                            (c) Fernando Puente-Sánchez, CNB-CSIC.

Generate tabular outputs from SqueezeMeta results.

USAGE: sqm2tables.py [-h] project_path output_dir
                     [--trusted-functions] [--ignore-unclassified]
                     [--sqm2anvio] [--force-overwrite]

OPTIONS:
    --trusted-functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables
    --ignore-unclassified: Ignore ORFs without assigned functions in
        TPM calculation
    --sqm2anvio: Write the required files for sqm2anvio
    --force-overwrite: Write results even if the output directory
        already exists
"""

from os.path import abspath, dirname, realpath
from os import mkdir, listdir
from sys import exit
import argparse

from collections import defaultdict

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.insert(0, '{}/../lib/'.format(utils_home))
data_dir = '{}/../data'.format(utils_home)

from utils import parse_conf_file, parse_orf_table, parse_tax_table, parse_contig_table, parse_bin_table, parse_tax_string, read_orf_names, aggregate_tax_abunds, normalize_abunds, write_orf_seqs, write_contig_seqs, write_row_dict, TAXRANKS 


def main(args):
    ### Get result files paths from SqueezeMeta_conf.pl
    perlVars = parse_conf_file(args.project_path, override = {'$projectdir': args.project_path})
    nokegg, nocog, nopfam, doublepass = map(int, [perlVars['$nokegg'], perlVars['$nocog'], perlVars['$nopfam'], perlVars['$doublepass']])

    ### Create output dir.
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise
        elif args.sqm2anvio or args.force_overwrite: # We know what we are doing.
            pass
        else:
            print('\nThe directory {} already exists. Please remove it or use a different output name.\n'.format(args.output_dir))
            exit(1)

    ### Calculate tables and write results.
    prefix = args.output_dir + '/' + perlVars['$projectname'] + '.'
    
    ### Functions
    if not args.sqm2anvio:
        # Were custom annotation databases used in this project?
        methods = [f.split('.')[-1] for f in listdir(perlVars['$resultpath']) if len(f.split('.'))>2 and f.split('.')[-2] == 'fun3']
        customMethods = [method for method in methods if method not in ('kegg', 'cog', 'pfam', 'wranks')]

        # Parse ORF table.
        sampleNames, orfs, kegg, cog, pfam, custom = parse_orf_table(perlVars['$mergedfile'], nokegg, nocog, nopfam, args.trusted_functions, args.ignore_unclassified, customMethods, data_dir)

        # Round aggregated functional abundances.
        # We can have non-integer abundances bc of the way we split counts in ORFs with multiple KEGGs.
        # We round the aggregates to the closest integer for convenience.
        kegg['abundances'] = {k: a.round().astype(int) for k,a in kegg['abundances'].items()}
        cog['abundances']  = {k: a.round().astype(int) for k,a in  cog['abundances'].items()}
        pfam['abundances'] = {k: a.round().astype(int) for k,a in pfam['abundances'].items()}

        #write_row_dict(sampleNames, orfs['tpm'], prefix + 'orf.tpm.tsv')
        if not nokegg:
            write_row_dict(['Name', 'Path'], kegg['info'], prefix + 'KO.names.tsv')
            write_row_dict(sampleNames, kegg['abundances'], prefix + 'KO.abund.tsv')
            write_row_dict(sampleNames, kegg['bases'], prefix + 'KO.bases.tsv')
            write_row_dict(sampleNames, kegg['tpm'], prefix + 'KO.tpm.tsv')
            if 'copyNumber' in kegg:
                write_row_dict(sampleNames, kegg['copyNumber'], prefix + 'KO.copyNumber.tsv')
        if not nocog:
            write_row_dict(['Name', 'Path'], cog['info'], prefix + 'COG.names.tsv')
            write_row_dict(sampleNames, cog['abundances'], prefix + 'COG.abund.tsv')
            write_row_dict(sampleNames, cog['bases'], prefix + 'COG.bases.tsv')
            write_row_dict(sampleNames, cog['tpm'], prefix + 'COG.tpm.tsv')
            if 'copyNumber' in cog:
                write_row_dict(sampleNames, cog['copyNumber'], prefix + 'COG.copyNumber.tsv')
                write_row_dict(sampleNames, {'COG0468': cog['coverages']['COG0468']}, prefix + 'RecA.tsv')
        if not nopfam:
            write_row_dict(sampleNames, pfam['abundances'], prefix + 'PFAM.abund.tsv')
            write_row_dict(sampleNames, pfam['bases'], prefix + 'PFAM.bases.tsv')
            write_row_dict(sampleNames, pfam['tpm'], prefix + 'PFAM.tpm.tsv')
            if 'copyNumber' in pfam:
                write_row_dict(sampleNames, pfam['copyNumber'], prefix + 'PFAM.copyNumber.tsv')
        for method, d in custom.items():
            write_row_dict(['Name'], d['info'], prefix + method + '.names.tsv')
            write_row_dict(sampleNames, d['abundances'], prefix + method + '.abund.tsv')
            write_row_dict(sampleNames, d['bases'], prefix + method + '.bases.tsv')
            write_row_dict(sampleNames, d['tpm'], prefix + method + '.tpm.tsv')
            if 'copyNumber' in d:
                write_row_dict(sampleNames, d['copyNumber'], prefix + method + '.copyNumber.tsv')
 
    else:
        # Not super beautiful code. Just read the orf names and create a fake orf dict
        # since we need to know the names of all the orfs to create the taxonomy output.
        orfs = {'abundances': read_orf_names(perlVars['$mergedfile'])}

    ### Taxonomy.
    fun_prefix = perlVars['$fun3tax_blastx'] if doublepass else perlVars['$fun3tax']
    orf_tax, orf_tax_wranks = parse_tax_table(fun_prefix + '.wranks')
    orf_tax_nofilter, orf_tax_nofilter_wranks = parse_tax_table(fun_prefix + '.noidfilter.wranks')

    # Add ORFs not present in the input tax file.
    unclass_list, unclass_list_wranks = parse_tax_string('n_Unclassified')
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

    write_row_dict(TAXRANKS, orf_tax, prefix + 'orf.tax.allfilter.tsv')
    write_row_dict(TAXRANKS, contig_tax, prefix + 'contig.tax.tsv')

    if not args.sqm2anvio:
        fna_blastx = perlVars['$fna_blastx'] if doublepass else None
        write_orf_seqs(orfs['abundances'].keys(), perlVars['$aafile'], fna_blastx, perlVars['$rnafile'], perlVars['$trnafile'] + '.fasta', prefix + 'orf.sequences.tsv')
        write_contig_seqs(perlVars['$contigsfna'], prefix + 'contig.sequences.tsv')

        write_row_dict(TAXRANKS, orf_tax_nofilter, prefix + 'orf.tax.nofilter.tsv')
        write_row_dict(TAXRANKS, orf_tax_prokfilter, prefix + 'orf.tax.prokfilter.tsv')

        ### Bins
        if not int(perlVars['$nobins']):
            bin_tpm, bin_tax, bin_tax_wranks = parse_bin_table(perlVars['$bintable'])
            write_row_dict(TAXRANKS, bin_tax, prefix + 'bin.tax.tsv')

        for idx, rank in enumerate(TAXRANKS):
            tax_abunds_orfs = aggregate_tax_abunds(orfs['abundances'], orf_tax_prokfilter_wranks, idx)
            write_row_dict(sampleNames, tax_abunds_orfs, prefix + '{}.prokfilter.abund.tsv'.format(rank))
            #write_row_dict(sampleNames, normalize_abunds(tax_abunds_orfs, 100), prefix + '{}.prokfilter.percent.tsv'.format(rank))

            tax_abunds_contigs = aggregate_tax_abunds(contig_abunds, contig_tax_wranks, idx)
            write_row_dict(sampleNames, tax_abunds_contigs, prefix + '{}.allfilter.abund.tsv'.format(rank))
            #write_row_dict(sampleNames, normalize_abunds(tax_abunds_contigs, 100), prefix + '{}.allfilter.percent.tsv'.format(rank))



def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate SqueezeMeta results into tables', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('--trusted-functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--ignore-unclassified', action='store_true', help='Ignore ORFs without assigned functions in TPM calculation')
    parser.add_argument('--sqm2anvio', action='store_true', help='Write the required files for sqm2anvio')
    parser.add_argument('--force-overwrite', action='store_true', help='Write results even if the output directory already exists')

    return parser.parse_args()



if __name__ == '__main__':
    main(parse_args())
