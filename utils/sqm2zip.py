#!/usr/bin/env python3

"""
Part of the SqueezeMeta distribution. 02/03/2023.
    (c) Fernando Puente-Sánchez, 2019-2020, CNB-CSIC / 2021-2023 SLU.

Package the essential files from a SqueezeMeta project into a single zip file for easy transfer.
The resulting zip file can be then be loaded directly into SQMtools.

USAGE: sqm2tables.py [-h] project_path output_dir
                     [--trusted-functions] [--ignore-unclassified]
                     [--force-overwrite] [--doc]

OPTIONS:
    --trusted-functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables. This will be ignored
        if the project_path/results/tables directory already exists
    --ignore-unclassified: Ignore reads without assigned functions in
        TPM calculation. This will be ignored if the
        project_path/results/tables directory already exists
    --force-overwrite: Write results even if the output file
        already exists
    --doc: show this documentation
"""

from os.path import abspath, realpath, dirname, exists
from os import listdir
from sys import exit, argv
from subprocess import call
import argparse
from zipfile import ZipFile

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.insert(0, '{}/../lib/'.format(utils_home))

from utils import parse_conf_file

def main(args):

    ### Parse project variables
    perlVars = parse_conf_file(args.project_path, override = {'$projectdir': args.project_path})
    project_name = perlVars['$projectname']
    output = f'{args.output_dir}/{project_name}.zip' 

    ### Check if the output file exists
    if exists(output):
        if args.force_overwrite:
            call(['rm','-r', output])
        else:
            print('\nThe file {} already exists. Please remove it or use a different output name.\n'.format(output))
            exit(1)

    ### Check whether we need to create the tables for this project
    if not exists(f'{args.project_path}/results/tables/{project_name}.phylum.nofilter.abund.tsv'):
        print(f'Generating tabular outputs for project in {args.project_path}')
        sqm2tables_args = ['python3', f'{utils_home}/sqm2tables.py', args.project_path, f'{args.project_path}/results/tables']
        if args.trusted_functions:
            sqm2tables_args.append('--trusted-functions')
        if args.ignore_unclassified:
            sqm2tables_args.append('--ignore-unclassified')
        call(sqm2tables_args)

    ### Which files will we include?
    target_files = ['creator.txt',
                    'SqueezeMeta_conf.pl',
                    f'results/10.{project_name}.mappingstat',
                    f'results/13.{project_name}.orftable',
                    f'results/19.{project_name}.contigtable']

    target_files += [f'results/tables/{f}' for f in listdir(f'{args.project_path}/results/tables')]

    if not int(perlVars['$nobins']) and exists(perlVars['$bintable']):
        target_files += [f'results/18.{project_name}.bintable', f'intermediate/18.{project_name}.contigsinbins']

    # ... and go for it
    with ZipFile(output, 'w') as outzip:
        for f in target_files:
            outzip.write(f'{args.project_path}/{f}', arcname=f)






def parse_args():
    parser = argparse.ArgumentParser(description='Package a SqueezeMeta project into a zip file', epilog='Fernando Puente-Sánchez (CNB-SLU) 2023\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('--trusted-functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--ignore-unclassified', action='store_true', help='Ignore ORFs without assigned functions in TPM calculation')
    parser.add_argument('--force-overwrite', action='store_true', help='Write results even if the output directory already exists')
    parser.add_argument('--doc', action='store_true', help='Show documentation')
    args = parser.parse_args()
    return args



if __name__ == '__main__':
    if '--doc' in argv: # hack so we can pass only --doc without getting an error for not providing the required positional arguments
        print(__doc__)
    else:
        main(parse_args())

