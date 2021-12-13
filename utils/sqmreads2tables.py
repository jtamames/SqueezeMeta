#!/usr/bin/env python3
"""
Part of the SqueezeMeta distribution. 22/07/2021.
    (c) Fernando Puente-Sánchez, 2019-2020, CNB-CSIC / 2021 SLU.

Generate tabular outputs from sqm_reads.pl or sqm_longreads.pl results.

USAGE: sqm_reads2tables.py [-h] project_path output_dir [-q "QUERY"]
                     [--trusted-functions] [--ignore-unclassified]
                     [--doc]

OPTIONS:
    -q/--query: Optional query for filtering your results (see below)
    --trusted-functions: Include only ORFs with highly trusted KEGG and
        COG assignments in aggregated functional tables
    --force-overwrite: Write results even if the output directory
        already exists
    --doc: Show this documentation

QUERY SYNTAX:

- Please enclose query strings within double brackets.

- Queries are combinations of relational operations in the form of
  <SUBJECT> <OPERATOR> <VALUE> (e.g. "PHYLUM == Bacteroidetes")
  joined by logical operators (AND, OR).

- Values are case-sensitive
  
- Parentheses can be used to group operations together.

- The "AND" and "OR" logical operators can't appear together in the same
  expression. Parentheses must be used to separate them into different
  expressions. e.g:
   "GENUS == Escherichia OR GENUS == Prevotella AND FUN CONTAINS iron"
        would not be valid. Parentheses must be used to write either:
   "(GENUS == Escherichia OR GENUS == Prevotella)" AND FUN CONTAINS iron"
        -> reads from either Escherichia or Prevotella which contain
           annotations related to iron.
   "GENUS == Escherichia OR (GENUS == Prevotella AND FUN CONTAINS iron)"
        -> reads splits from Escherichia, and splits of Prevotella which
           contains annotations related to iron.
           
- Another example query would be:
  "(PHYLUM == Bacteroidetes
    OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria])
    AND FUN CONTAINS iron"
    - This would select all the reads assigned to either the
      Bacteroidetes phylum or the Alphaproteobacteria or
      Gammaproteobacteria classes, that also contain the substring
      "iron" in the functional annotations.

- Possible subjects are:
    - FUN: search within all the combined databases used for functional
           annotation.
    - FUNH: search within the KEGG BRITE and COG functional hierarchies
         (e.g. "FUNH CONTAINS Carbohydrate metabolism" will select all
          the splits containing a gene associated with the broad
          "Carbohydrate metabolism" category)
    - SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES:
          search within the taxonomic annotation at the requested
          taxonomic rank.
"""

from os.path import abspath, dirname, realpath
from os import mkdir, listdir
from sys import exit, argv
import argparse

from collections import defaultdict
from pandas import DataFrame

from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.insert(0, '{}/../lib/'.format(utils_home))
data_dir = ('{}/../data'.format(utils_home))

from utils import parse_tax_string
from splitFilter import DictFilter

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
            sample, file_, total_reads, reads_with_hits_to_nr, *total_hits = line.strip().split('\t') # *_ since longreads output will have one more column
            if total_hits:
                longreads = True
            else:
                longreads = False
            samples[sample] += int(total_reads)
 

    # Is there any custom annotation method apart from kegg and COG?
    custom_methods = [f.split('.')[-1].replace('fun', '') for f in listdir(args.project_path) if 'fun' in f.split('.')[-1] and not f.endswith('funcog') and not f.endswith('funkegg')]
    for method in custom_methods:
        FUNMETHODS[method] = method


    tax_dict = {filt: {rank: {sample: defaultdict(int) for sample in samples} for rank in TAXRANKS} for filt in TAXFILTERS}
    fun_dict = {method: {sample: defaultdict(float) for sample in samples} for method in FUNMETHODS}
    fun_info = {}

    # Parse fun names and hierarchy paths for kegg/cog.
    for method, method_name in FUNMETHODS.items():
        if method not in ('kegg', 'cogs'):
            continue
        function_info = 'keggfun2.txt' if method == 'kegg' else 'coglist.txt'
        with open('{}/{}'.format(data_dir, function_info)) as infile:
            infile.readline() # Burn headers.
            fun_info[method] = {}
            for line in infile:
                if method == 'kegg':
                    fun_id, gene_name, fun_name, path = line.strip().split('\t')
                else:
                    line = line.strip().split('\t')
                    if len(line) == 3:
                        fun_id, fun_name, path = line
                    else: # UGH!
                        fun_id, fun_name = line
                        path = '{} (path not available)',format(fun_id)
                fun_info[method][fun_id] = (fun_name, path)

    # Parse function names for extra methods.
    for method, method_name in FUNMETHODS.items():
        if method in ('kegg', 'cogs'):
            continue
        method_name = FUNMETHODS[method]
        fun_info[method] = {}
        if method in ('kegg', 'cogs'):
            continue
        with open('{}/{}.out.allreads.fun{}'.format(args.project_path, project_name, method)) as infile:
            infile.readline() # Burn headers.
            infile.readline()
            for line in infile:
                line = line.strip('\n').split('\t') # Explicitly strip just '\n' so I don't remove tabs when there are empty fields.
                fun_id, fun_name = line[0], line[-1]
                fun_info[method][fun_id] = fun_name
 

    ### Process samples
    def tax_file_to_dict(path, out_dict):
        with open(path) as infile:
            for line in infile:
                fields = line.strip().split('\t')
                if not longreads:
                     read = 'FILE={}|READ={}'.format(path.rsplit('.', 2)[0], fields[0]) # store the path + name of the original file without suffixes
                                                                                        # so standard and nofilter look the same, but fwd and rev don't
                else:
                     read = fields[0]
                if not longreads:
                    tax_string = fields[1] if len(fields) > 1 else 'n_Unclassified'
                else:
                    tax_string = fields[1] if fields[1] else 'n_Unclassified'
                tax, tax_wranks = parse_tax_string(tax_string)
                assert read not in out_dict
                out_dict[read] = tax_wranks
    
    # Go anti go!
    filterOK = []
    for sample in samples:

        read_tax = {'nofilter': {}, 'allfilter': {}, 'prokfilter': {}}
        tax_reads = set()

        # Parse nofilter taxonomy
        if not longreads:
            nofilter_tax_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith('.tax_nofilter.wranks')]
        else:
            nofilter_tax_files = ['readconsensus_nofilter.txt']
        for tax_file in nofilter_tax_files:
            path = '{}/{}/{}'.format(args.project_path, sample, tax_file)
            tax_file_to_dict(path, read_tax['nofilter'])

        # Parse taxonomy with filters.
        if not longreads:
            allfilter_tax_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith('.tax.wranks')]
        else:
            allfilter_tax_files = ['readconsensus.txt']
        for tax_file in allfilter_tax_files:
            path = '{}/{}/{}'.format(args.project_path, sample, tax_file)
            tax_file_to_dict(path, read_tax['allfilter'])
        assert read_tax['nofilter'].keys() == read_tax['allfilter'].keys()

        # Generate taxonomy with filters only for prokaryotes.
        for read in read_tax['nofilter']:
            if 'k_Bacteria' in (read_tax['nofilter'][read][0], read_tax['allfilter'][read][0]) or 'k_Archaea' in (read_tax['nofilter'][read][0], read_tax['allfilter'][read][0]):
                read_tax['prokfilter'][read] = read_tax['allfilter'][read]
            else:
                read_tax['prokfilter'][read] = read_tax['nofilter'][read]

        for read in read_tax['prokfilter']:
            tax_reads.add(read)       

        read_fun = {method: {} for method in FUNMETHODS}
        fun_reads = set()

        # Parse functions.
        for method in FUNMETHODS:
            fun_files = [f for f in listdir('{}/{}'.format(args.project_path, sample)) if f.endswith(method)]
            for fun_file in fun_files:
                path = '{}/{}/{}'.format(args.project_path, sample, fun_file)
                with open(path) as infile:
                    infile.readline()
                    infile.readline()
                    for line in infile:
                        fields = line.strip().split('\t')
                        while len(fields) < 3:
                            fields.append('Unclassified')
                        if not longreads:
                            read = 'FILE={}|READ={}'.format(path.rsplit('.', 1)[0], fields[0]) ## store the path + name of the original file (without the .kegg .cogs... suffixes)
                        else:
                            read = fields[0]
                        funs = fields[2] if args.trusted_functions else fields[1]
                        funs = funs.split(';')
                        assert read not in read_fun[method]
                        read_fun[method][read] = funs
                        fun_reads.add(read)
        if longreads:
            # Get all fun reads even if they are not annotated
            with open('{}/{}.out.allreads'.format(args.project_path, project_name)) as infile:
                infile.readline()
                infile.readline()
                for line in infile:
                    s, _, read, *__ = line.split('\t')
                    if s == sample:
                        if s not in fun_reads:
                            fun_reads.add(read) # With this, fun_reads has ALL the reads/ORFs

        # Propagate unclassified funs
        for read in fun_reads:
            for method in FUNMETHODS:
                if read not in read_fun[method]:
                    read_fun[method][read] = ['Unclassified']

        # Propagate reads from fun to tax and vice-versa (shortreads only)
        if not longreads:
            for filt in TAXFILTERS:
                for read in fun_reads:
                    if read not in read_tax[filt]:
                        read_tax[filt][read] = parse_tax_string('n_Unclassified')[1]
                        tax_reads.add(read)
            for method in FUNMETHODS:
                for read in tax_reads:
                    if read not in read_fun[method]:
                       read_fun[method][read] = ['Unclassified']
                       fun_reads.add(read)

        # Tax_fun name equivalence (for longreads)
        tax2fun = defaultdict(list)
        fun2tax = {}
        for read in fun_reads:
            tread = read.rsplit('_', 1)[0] if longreads else read
            tax2fun[tread].append(read)
            fun2tax[read] = tread


        # Select reads.
        if args.query:
            sfilter = DictFilter(read_tax, 'prokfilter', read_fun, fun_info, tax2fun, fun2tax)
            print('- Query tree is:')
            sfilter.print_tree(args.query)
            gr = sfilter.run(args.query) # here we also take care of fun-tax names equivalence if working with longreads
            good_tax_reads = set()
            good_fun_reads = set()
            for read in gr:
                if read in tax_reads:
                    good_tax_reads.add(read)
                if read in fun_reads:
                    good_fun_reads.add(read)
            if gr:
                filterOK.append(True)
            else:
                filterOK.append(False)
        else:
            good_tax_reads = tax_reads
            good_fun_reads = fun_reads
            filterOK.append(True)

        # Aggregate counts from the same taxa. 
        for filt in TAXFILTERS:
            for read, tax in read_tax[filt].items():
                if read not in good_tax_reads:
                    continue
                for i, rank in enumerate(TAXRANKS):
                    tax_dict[filt][rank][sample][tax[i]] += 1
        for i, rank in enumerate(TAXRANKS): # And add unclassified
            for filt in TAXFILTERS:
                taxa = tax_dict[filt][rank][sample]
                if not args.query:
                    uString = parse_tax_string('n_Unclassified')[1][i]
                    classified_reads = sum(taxa.values()) - taxa[uString]
                    total_reads = samples[sample]
                    taxa[uString] = (total_reads - classified_reads)
            
                        
        # Aggregate counts from the same function.
        for method in FUNMETHODS:
            for read, funs in read_fun[method].items():
                if read not in good_fun_reads:
                    continue
                for fun in funs:
                    fun_dict[method][sample][fun] += 1/len(funs) # Split the counts in multi-kegg annotations.
            functions = fun_dict[method][sample]
            classified_reads = sum(functions.values())
            if not args.query:
                if not longreads:
                    classified_reads = sum(functions.values()) - functions['Unclassified']
                    total_reads = samples[sample]
                    functions['Unclassified'] = (total_reads - classified_reads)
                else:
                    pass # we already counted all reads

    ### Test whether we got something
    if not any(filterOK):
        print('\nYour query generated no positive results\n')
        exit(0)

    ### Write tax results.
    for filt in TAXFILTERS:
        for i, rank in enumerate(TAXRANKS):
            dict_to_write = tax_dict[filt][rank]
            DataFrame.from_dict(dict_to_write).fillna(0).to_csv('{}/{}.{}.{}.abund.tsv'.format(args.output_dir, output_prefix, rank, filt), sep='\t')


    ### Write fun results.
    for method, method_name in FUNMETHODS.items():
        dict_to_write = fun_dict[method]
        DataFrame.from_dict(dict_to_write).fillna(0).to_csv('{}/{}.{}.abund.tsv'.format(args.output_dir, output_prefix, method_name), sep='\t')


    # Write function names and hierarchy paths
    for method, method_name in FUNMETHODS.items():
        allFuns = sorted({fun for sample in fun_dict[method] for fun in fun_dict[method][sample]})
        info = fun_info[method]
        with open('{}/{}.{}.names.tsv'.format(args.output_dir, output_prefix, method_name), 'w') as outfile:
            outfile.write('\tName\tPath\n')
            for fun in allFuns:
                if fun == 'Unclassified':
                    continue
                if method in ('kegg', 'cogs'):
                    if fun in info:
                        outfile.write('{}\t{}\t{}\n'.format(fun, info[fun][0], info[fun][1]))
                    else:
                        outfile.write('{}\t{} (name not available)\t{} (path not available)\n'.format(fun, fun, fun))
                else:
                    outfile.write('{}\t{}\t\n'.format(fun, info[fun]))

            if method in ('kegg', 'cogs'):
               outfile.write('Unclassified\tUnclassified\tUnclassified\n')
            else:
               outfile.write('Unclassified\tUnclassified\t\n')
            



def parse_args():
    parser = argparse.ArgumentParser(description='Aggregate SqueezeMeta results into tables', epilog='Fernando Puente-Sánchez (CNB-SLU) 2021\n')
    parser.add_argument('project_path', type=str, help='Base path of the SqueezeMeta project')
    parser.add_argument('output_dir', type=str, help='Output directory')
    parser.add_argument('-q', '--query', type=str, nargs='+', help='Query for filtering the results')
    parser.add_argument('--trusted-functions', action='store_true', help='Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables')
    parser.add_argument('--force-overwrite', action='store_true', help='Write results even if the output directory already exists')
    parser.add_argument('--doc', action='store_true', help='Show documentation')
    args = parser.parse_args()
    if args.query:
        args.query = ' '.join(args.query)
    return args



if __name__ == '__main__':
    if '--doc' in argv: # hack so we can pass only --doc without getting an error for not providing the required positional arguments
        print(__doc__)
    else:
        main(parse_args())

