#!/usr/bin/env python3

"""
Part of the SqueezeMeta distribution. 10/05/2021.
    (c) Fernando Puente-Sánchez, 2019-2020, CNB-CSIC / 2021 SLU.

Generate tabular outputs from SqueezeMeta results.

USAGE: sqm2tables.py [-h] project_path output_dir
                     [--trusted-functions] [--ignore-unclassified]
                     [--sqm2anvio] [--force-overwrite] [--doc]

OPTIONS:
    --trusted-functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables
    --ignore-unclassified: Ignore reads without assigned functions in
        TPM calculation
    --sqm2anvio: Write the required files for sqm2anvio
    --force-overwrite: Write results even if the output directory
        already exists
    --doc: show this documentation
"""

from os.path import abspath, dirname, realpath, isfile
from os import mkdir, listdir
from sys import exit, argv
import argparse

from collections import defaultdict

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.insert(0, '{}/../lib/'.format(utils_home))
data_dir = '{}/../data'.format(utils_home)

from utils import parse_conf_file, parse_mappingstat, parse_orf_table, parse_tax_table, parse_contig_table, parse_contig_tax, parse_bin_table, parse_tax_string, read_orf_names, aggregate_tax_abunds, normalize_abunds, write_orf_seqs, write_contig_seqs, write_row_dict, TAXRANKS, TAXRANKS_SHORT 


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
        elif not any(f.endswith('.superkingdom.allfilter.abund.tsv') for f in listdir(args.output_dir)):
            pass
        else:
            print('\nThe directory {} already contains results. Please remove it or use a different output name.\n'.format(args.output_dir))
            exit(1)

    prefix = args.output_dir + '/' + perlVars['$projectname'] + '.'

    ### Load total reads and bases
    total_reads, total_bases = parse_mappingstat(perlVars['$mappingstat'])
    
    ### Functions
    if not args.sqm2anvio:
        # Were custom annotation databases used in this project?
        methods = [f.split('.')[-1] for f in listdir(perlVars['$resultpath']) if len(f.split('.'))>2 and f.split('.')[-2] == 'fun3']
        customMethods = [method for method in methods if method not in ('kegg', 'cog', 'pfam', 'wranks')]

        # Parse ORF table.
        sampleNames, orfs, kegg, cog, pfam, custom, noCDSorfs, noCDScontigs = parse_orf_table(perlVars['$mergedfile'], total_reads, total_bases, nokegg, nocog, nopfam,
                                                                                              args.trusted_functions, args.ignore_unclassified, customMethods, data_dir)

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
            write_row_dict(sampleNames, kegg['coverages'], prefix + 'KO.cov.tsv')
            write_row_dict(sampleNames, kegg['tpm'], prefix + 'KO.tpm.tsv')
            if 'copyNumber' in kegg:
                write_row_dict(sampleNames, kegg['copyNumber'], prefix + 'KO.copyNumber.tsv')
        if not nocog:
            write_row_dict(['Name', 'Path'], cog['info'], prefix + 'COG.names.tsv')
            write_row_dict(sampleNames, cog['abundances'], prefix + 'COG.abund.tsv')
            write_row_dict(sampleNames, cog['bases'], prefix + 'COG.bases.tsv')
            write_row_dict(sampleNames, cog['coverages'], prefix + 'COG.cov.tsv')
            write_row_dict(sampleNames, cog['tpm'], prefix + 'COG.tpm.tsv')
            if 'copyNumber' in cog:
                write_row_dict(sampleNames, cog['copyNumber'], prefix + 'COG.copyNumber.tsv')
                write_row_dict(sampleNames, {'COG0468': cog['coverages']['COG0468']}, prefix + 'RecA.tsv')
        if not nopfam:
            write_row_dict(sampleNames, pfam['abundances'], prefix + 'PFAM.abund.tsv')
            write_row_dict(sampleNames, pfam['bases'], prefix + 'PFAM.bases.tsv')
            write_row_dict(sampleNames, pfam['coverages'], prefix + 'PFAM.cov.tsv')
            write_row_dict(sampleNames, pfam['tpm'], prefix + 'PFAM.tpm.tsv')
            if 'copyNumber' in pfam:
                write_row_dict(sampleNames, pfam['copyNumber'], prefix + 'PFAM.copyNumber.tsv')
        for method, d in custom.items():
            write_row_dict(['Name'], d['info'], prefix + method + '.names.tsv')
            write_row_dict(sampleNames, d['abundances'], prefix + method + '.abund.tsv')
            write_row_dict(sampleNames, d['bases'], prefix + method + '.bases.tsv')
            write_row_dict(sampleNames, d['coverages'], prefix + method + '.cov.tsv')
            write_row_dict(sampleNames, d['tpm'], prefix + method + '.tpm.tsv')
            if 'copyNumber' in d:
                write_row_dict(sampleNames, d['copyNumber'], prefix + method + '.copyNumber.tsv')
 
    else:
        # Not super beautiful code. Just read the orf names and create a fake orf dict
        # since we need to know the names of all the orfs to create the taxonomy output.
        orfs = {'abundances': read_orf_names(perlVars['$mergedfile'])}
        noCDSorfs = noCDScontigs = set()

    ### Taxonomy.
    fun_prefix = perlVars['$fun3tax_blastx'] if doublepass else perlVars['$fun3tax']
    orf_tax, orf_tax_wranks = parse_tax_table(fun_prefix + '.wranks', noCDS = noCDSorfs)
    orf_tax_nofilter, orf_tax_nofilter_wranks = parse_tax_table(fun_prefix + '.noidfilter.wranks', noCDS = noCDSorfs)

    contig_abunds, _, _ = parse_contig_table(perlVars['$contigtable'])
    contig_tax, contig_tax_wranks = parse_contig_tax(perlVars['$interdir'] + '/09.' + perlVars['$projectname'] + '.contiglog', noCDS = noCDScontigs)
    contig_tax_nofilter, contig_tax_nofilter_wranks = parse_contig_tax(perlVars['$interdir'] + '/09.' + perlVars['$projectname'] + '.contiglog.noidfilter', noCDS = noCDScontigs)
    
    # Add ORFs/contigs not present in the input tax file.
    
    def add_features(abunds, tax, tax_wranks, tax_nofilter, tax_nofilter_wranks, noCDS):
        for feat in abunds:
            if feat not in tax:
                assert feat not in tax_wranks
                assert feat not in tax_nofilter
                assert feat not in tax_nofilter_wranks
                if feat in noCDS:
                    unclass_list, unclass_list_wranks = parse_tax_string('n_No CDS', emptyClassString = 'No CDS')
                else:
                    unclass_list, unclass_list_wranks = parse_tax_string('n_Unclassified', emptyClassString = 'Unclassified')
                tax[feat] = unclass_list
                tax_wranks[feat] = unclass_list_wranks
                tax_nofilter[feat] = unclass_list
                tax_nofilter_wranks[feat] = unclass_list_wranks

    add_features(orfs['abundances'], orf_tax, orf_tax_wranks, orf_tax_nofilter, orf_tax_nofilter_wranks, noCDSorfs)
    add_features(contig_abunds, contig_tax, contig_tax_wranks, contig_tax_nofilter, contig_tax_nofilter_wranks, noCDScontigs)
    
    orf_tax_prokfilter, orf_tax_prokfilter_wranks = {}, {}
    contig_tax_prokfilter, contig_tax_prokfilter_wranks = {}, {}

    def calc_prokfilter(tax, tax_wranks, tax_nofilter, tax_nofilter_wranks, tax_prokfilter, tax_prokfilter_wranks):
        for feat in tax:
            tax_af = tax[feat]
            tax_nf = tax_nofilter[feat]
            if 'Bacteria' in (tax_af[0], tax_nf[0]) or 'Archaea' in (tax_af[0], tax_nf[0]): # We check both taxonomies.
                tax_prokfilter       [feat] = tax_af
                tax_prokfilter_wranks[feat] = tax_wranks[feat]
            else:
                tax_prokfilter       [feat] = tax_nf
                tax_prokfilter_wranks[feat] = tax_nofilter_wranks[feat]

    calc_prokfilter(orf_tax, orf_tax_wranks, orf_tax_nofilter, orf_tax_nofilter_wranks, orf_tax_prokfilter, orf_tax_prokfilter_wranks)
    calc_prokfilter(contig_tax, contig_tax_wranks, contig_tax_nofilter, contig_tax_nofilter_wranks, contig_tax_prokfilter, contig_tax_prokfilter_wranks)
    
    write_row_dict(TAXRANKS, orf_tax, prefix + 'orf.tax.allfilter.tsv')
    write_row_dict(TAXRANKS, orf_tax_prokfilter, prefix + 'orf.tax.prokfilter.tsv')
    write_row_dict(TAXRANKS, orf_tax_nofilter, prefix + 'orf.tax.nofilter.tsv')
    write_row_dict(TAXRANKS, contig_tax, prefix + 'contig.tax.allfilter.tsv')
    write_row_dict(TAXRANKS, contig_tax_prokfilter, prefix + 'contig.tax.prokfilter.tsv')
    write_row_dict(TAXRANKS, contig_tax_nofilter, prefix + 'contig.tax.nofilter.tsv')


    if not args.sqm2anvio:
        fna_blastx = perlVars['$fna_blastx'] if doublepass else None
        write_orf_seqs(orfs['abundances'].keys(), perlVars['$aafile'], fna_blastx, perlVars['$rnafile'], perlVars['$trnafile'] + '.fasta', prefix + 'orf.sequences.tsv')
        write_contig_seqs(perlVars['$contigsfna'], prefix + 'contig.sequences.tsv')

        write_row_dict(TAXRANKS, orf_tax_nofilter, prefix + 'orf.tax.nofilter.tsv')
        write_row_dict(TAXRANKS, orf_tax_prokfilter, prefix + 'orf.tax.prokfilter.tsv')

        ### Bins
        if not int(perlVars['$nobins']) and isfile(perlVars['$bintable']):
            bin_tpm, bin_tax, bin_tax_wranks = parse_bin_table(perlVars['$bintable'])
            write_row_dict(TAXRANKS, bin_tax, prefix + 'bin.tax.tsv')

        for idx, rank in enumerate(TAXRANKS):
            unmapped_str = ';'.join([rs+'_Unmapped' for i, rs in enumerate(TAXRANKS_SHORT) if i <= idx])

            tax_abunds_contigs = aggregate_tax_abunds(contig_abunds, contig_tax_wranks, idx)
            tax_abunds_contigs[unmapped_str] = total_reads - sum(tax_abunds_contigs.values())
            write_row_dict(sampleNames, tax_abunds_contigs, prefix + '{}.allfilter.abund.tsv'.format(rank))

            tax_abunds_contigs = aggregate_tax_abunds(contig_abunds, contig_tax_prokfilter_wranks, idx)
            tax_abunds_contigs[unmapped_str] = total_reads - sum(tax_abunds_contigs.values())
            write_row_dict(sampleNames, tax_abunds_contigs, prefix + '{}.prokfilter.abund.tsv'.format(rank))

            tax_abunds_contigs = aggregate_tax_abunds(contig_abunds, contig_tax_nofilter_wranks, idx)
            tax_abunds_contigs[unmapped_str] = total_reads - sum(tax_abunds_contigs.values())
            write_row_dict(sampleNames, tax_abunds_contigs, prefix + '{}.nofilter.abund.tsv'.format(rank))




def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate SqueezeMeta results into tables', epilog='Fernando Puente-Sánchez (CNB-SLU) 2021\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('--trusted-functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--ignore-unclassified', action='store_true', help='Ignore ORFs without assigned functions in TPM calculation')
    parser.add_argument('--sqm2anvio', action='store_true', help='Write the required files for sqm2anvio')
    parser.add_argument('--force-overwrite', action='store_true', help='Write results even if the output directory already exists')
    parser.add_argument('--doc', action='store_true', help='Show documentation')

    return parser.parse_args()



if __name__ == '__main__':
    if '--doc' in argv: # hack so we can pass only --doc without getting an error for not providing the required positional arguments
        print(__doc__)
    else:
        main(parse_args())
