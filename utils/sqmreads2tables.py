#!/usr/bin/env python3
"""
Part of the SqueezeMeta distribution. 25/03/2018 Original version,
                            (c) Fernando Puente-Sánchez, CNB-CSIC.

Generate tabular outputs from sqm_reads results.

USAGE: sqm_reads2tables.py [-h] project_path output_dir
                     [--trusted-functions] [--ignore-unclassified]

OPTIONS:
    --trusted-functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables
    --force-overwrite: Write results even if the output directory
        already exists.
"""

from os.path import abspath, dirname, realpath
from os import mkdir, listdir
from sys import exit
import argparse

from collections import defaultdict
from pandas import DataFrame

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.insert(0, '{}/../lib/'.format(utils_home))
from utils import parse_tax_string

TAXFILTERS = ('nofilter', 'allfilter', 'prokfilter')
TAXRANKS   = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
FUNMETHODS = {'kegg': 'KO', 'cogs': 'COG'}

def main(args):
    ### Create output dir.
    try:
       mkdir(args.output_dir)
    except OSError as e:
        if e.errno != 17:
            raise
        elif args.force_overwrite:
            pass
        else:
            print('\nThe directory {} already exists. Please remove it or use a different output name.\n'.format(args.output_dir))
            exit(1)


    ### Project name and samples.
    project_name = args.project_path.strip('/').split('/')[-1]
    output_prefix = project_name #args.output_dir.strip('/').split('/')[-1]
    samples = defaultdict(int)
    with open('{}/{}.out.mappingstat'.format(args.project_path, project_name)) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            sample, file_, total_reads, reads_with_hits_to_nr = line.strip().split('\t')
            samples[sample] += int(total_reads)
 

    ### Parse taxonomy.

    def tax_file_to_dict(path, out_dict):
        with open(path) as infile:
            for line in infile:
                fields = line.strip().split('\t')
                # .replace('_nofilter') so that reads from allfilter and nofilter have the same name, but reads from fwd and rev don't.
                read = path.replace('_nofilter', '') + fields[0]
                tax_string = fields[1] if len(fields) > 1 else 'n_Unclassified'
                tax, tax_wranks = parse_tax_string(tax_string)
                out_dict[read] = tax_wranks


    tax_dict = {filt: {rank: {sample: defaultdict(int) for sample in samples} for rank in TAXRANKS} for filt in TAXFILTERS}

    for sample in samples:

        read_tax = {'nofilter': {}, 'allfilter': {}, 'prokfilter': {}}

        ### Parse nofilter taxonomy.
        nofilter_tax_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith('.tax_nofilter.wranks')]
        for tax_file in nofilter_tax_files:
            path = '{}/{}/{}'.format(args.project_path, sample, tax_file)
            tax_file_to_dict(path, read_tax['nofilter'])

        ### Parse taxonomy with filters.
        allfilter_tax_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith('.tax.wranks')]
        for tax_file in allfilter_tax_files:
            path = '{}/{}/{}'.format(args.project_path, sample, tax_file)
            tax_file_to_dict(path, read_tax['allfilter'])

        assert read_tax['nofilter'].keys() == read_tax['allfilter'].keys()

        ### Generate taxonomy with filters only for prokaryotes.
        for read in read_tax['nofilter']:
           if 'k_Bacteria' in (read_tax['nofilter'][read][0], read_tax['allfilter'][read][0]) or 'k_Archaea' in (read_tax['nofilter'][read][0], read_tax['allfilter'][read][0]):
               read_tax['prokfilter'][read] = read_tax['allfilter'][read]
           else:
               read_tax['prokfilter'][read] = read_tax['nofilter'][read]
              
        ### Aggregate counts from the same taxa. 
        for filt in TAXFILTERS:
            for read, tax in read_tax[filt].items():
                for i, rank in enumerate(TAXRANKS):
                    tax_dict[filt][rank][sample][tax[i]] += 1
   

    ### Add unclassified and write results.
    for filt in TAXFILTERS:
        if filt == 'nofilter':
            continue
        for i, rank in enumerate(TAXRANKS):
            dict_to_write = tax_dict[filt][rank]
            for sample, taxa in dict_to_write.items():
                classified_reads = sum(taxa.values())
                total_reads = samples[sample]
                dict_to_write[sample][parse_tax_string('n_Unclassified')[1][i]] += (total_reads - classified_reads)
            DataFrame.from_dict(dict_to_write).fillna(0).to_csv('{}/{}.{}.{}.abund.tsv'.format(args.output_dir, output_prefix, rank, filt), sep='\t')



    ### Parse functions.
    fun_dict = {method: {sample: defaultdict(float) for sample in samples} for method in FUNMETHODS}
    for sample in samples:
        for method in FUNMETHODS:
            fun_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith(method)]
            for fun_file in fun_files:
                with open('{}/{}/{}'.format(args.project_path, sample, fun_file)) as infile:
                    infile.readline()
                    infile.readline()
                    for line in infile:
                        fields = line.strip().split('\t')
                        while len(fields) < 3:
                            fields.append('Unclassified')
                        funs = fields[2] if args.trusted_functions else fields[1]
                        funs = funs.split(';')
                        for fun in funs:
                            fun_dict[method][sample][fun] += 1/len(funs) # Split the counts in multi-kegg annotations.

    for method, method_name in FUNMETHODS.items():
        dict_to_write = fun_dict[method]
        for sample, funs in dict_to_write.items():
            classified_reads = sum(funs.values())
            total_reads = samples[sample]
            dict_to_write[sample]['Unclassified'] += (total_reads - classified_reads)
        DataFrame.from_dict(dict_to_write).fillna(0).to_csv('{}/{}.{}.abund.tsv'.format(args.output_dir, output_prefix, method_name), sep='\t')




def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate SqueezeMeta results into tables', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('--trusted-functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--force-overwrite', action='store_true', help='Write results even if the output directory already exists')

    return parser.parse_args()



if __name__ == '__main__':
    main(parse_args())

