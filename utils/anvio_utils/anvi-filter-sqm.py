#!/usr/bin/env python3

"""
Select the anvi'o splits (which are themselves contigs, or parts of
long contigs) that fulfill a set of conditions, and visualize them
using anvi-interactive.

USAGE:

anvi-filter-sqm.py [-h] -p PROFILE_DB -c CONTIGS_DB -t TAXONOMY
      -q "QUERY" [-o OUTPUT_DIR] [-m MAX_SPLITS] [--enforce-clustering]
      [--extra-anvio-args "EXTRA_ANVIO_ARGS"] [-s {safe,yolo}]
     

OPTIONS:

- The "-o/--output_dir" parameter controls the output directory for
  the filtered anvi'o databases. It defaults to "filteredDB".

- The "-m/--max-splits" parameter controls the maximum number of splits
  that will be loaded into anvi'o. If the provided query returns a
  higher number of splits, the program will stop. By default it is set
  to 25,000, larger values may make the anvi'o interface to respond
  slowly. Setting --max-splits to 0 will allow an arbitrarily large
  number of splits to be loaded.

- By default, splits are only clustered based on their taxonomy. Passing
  the "--enforce-clustering" flag will make anvi'o perform an
  additional clustering based on abundances across samples and sequence
  composition.

- Extra arguments for anvi-interactive can be passed with the parameter
  --extra-anvio-args "extra args". Extra arguments must be surrounded
  by quotes.
     e.g. --extra-anvio-args "--taxonomic-level t_phylum --title Parrot"
     
- By default, the script uses an in-house method to subset the anvi'o
  databases. It's ~5x quicker than using anvi-split in anvio5, and works well for
  us. However, the night is dark and full of bugs, so if you feel that
  your anvi'o view is missing some information, you can call the script
  with "-s safe" parameter. This will call anvi-split which should be
  much safer than our hacky solution.



QUERY SYNTAX:

- Please enclose query strings within double brackets.

- Queries are combinations of relational operations in the form of
  <SUBJECT> <OPERATOR> <VALUE> (e.g. "PHYLUM == Bacteroidetes")
  joined by logical operators (AND, OR).

- Parentheses can be used to group operations together.

- The "AND" and "OR" logical operators can't appear together in the same
  expression. Parentheses must be used to separate them into different
  expressions. e.g:
   "GENUS == Escherichia OR GENUS == Prevotella AND FUN CONTAINS iron"
        would not be valid. Parentheses must be used to write either:
   "(GENUS == Escherichia OR GENUS == Prevotella)" AND FUN CONTAINS iron"
        -> splits from either Escherichia or Prevotella which contain
           ORFs related to iron.
   "GENUS == Escherichia OR (GENUS == Prevotella AND FUN CONTAINS iron)"
        -> all splits from Escherichia, and splits of Prevotella which
           contains ORFs related to iron.

- Another example query would be:
  "(PHYLUM == Bacteroidetes
    OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria])
    AND FUN CONTAINS iron AND Sample1 > 1"
    - This would select all the anvi'o splits assigned to either the
      Bacteroidetes phylum or the Alphaproteobacteria or
      Gammaproteobacteria classes, that also contain the substring
      "iron" in the functional annotations of any of their ORFs, and
      whose anvi'o abundance (mean coverage of a split divided by
      overall sample mean coverage) in Sample1 is higher than 1.

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
    - <SAMPLE_NAME>: search within the abundances in the
          sample named <SAMPLE_NAME>
        (e.g. if you have two samples named "Sample1" and "Sample2", 
         the query string "Sample1 > 0.5 AND Sample2 > 0.5" would
         return the anvi'o splits with an anvi'o abundance higher than
         0.5 in both samples)

- Posible relational operators are "==", "!=", ">=", "<=", ">", "<",
                      "IN", "NOT IN", "CONTAINS", "DOES NOT CONTAIN"

"""

import argparse
from os.path import abspath, dirname, realpath, exists
from sys import path, stderr
import glob
utils_home = abspath(dirname(realpath(__file__)))
path.append('{}/../../lib/'.format(utils_home))
import dendropy
import sqlite3
from collections import defaultdict
from sys import argv
from subprocess import call
from os import devnull
DEVNULL = open(devnull, 'wb')
from splitFilter import SplitFilter
from subsetAnvio import subset_anvio


OUTTREE = 'Taxonomy.nwk'
COLLECTION_NAME = 'SqueezeMeta'

def main(args):
    
    ### Check that the output directory does not exist.
    outdir = args.output_dir.rstrip('/')
    if exists(outdir):
        print('\nThe directory {} already exists. Please remove it or use a different output name.\n'.format(args.output_dir))
        exit(1)
    if exists(outdir+'_YOLO'): # Due to an aborted run.
        call(['rm', '-r', outdir+'_YOLO'])

    ### Check that anvi'o can be found in PATH
    missingAnvi = call(['which', 'anvi-self-test'], stdout=DEVNULL)
    if missingAnvi: # ecode of which was 1
        print('We can\'t find anvi\'o in your PATH. Is the environment activated?')
        exit(1)

    ### Load data
    sfilter = SplitFilter(args.contigs_db, args.profile_db, args.taxonomy)
    splits, samples, contigTax = sfilter.get_parsed_annotations()
    print('')
    print('- {} splits found in anvi\'o profile database.'.format(len(splits)))
    print('')
    if samples:
        print('- Sample names are:')
        print('\t{}'.format('\t'.join(samples)))
        print('')
    else:
        print('- You provided a blank profile. Queries based on sample abundances will not work')

    ### Run filter
    print('- Query tree is:')
    sfilter.print_tree(args.query)
    print('')
    goodSplits = sfilter.run(args.query)
    print('- {} splits fulfilling the required conditions.'.format(len(goodSplits)))
    if not goodSplits:
        print ('Exiting')
        exit(1)
    if args.max_splits and len(goodSplits) > args.max_splits:
        print('')
        print('Too many splits satisfying the provided query.')
        print('Loading too many splits can make anvi\'o slow to respond. Please consider providing a narrower query.')
        print('Alternatively, use --max-splits NUMBER to increase this threshold or --max-splits 0 to disable it completely')
        print('')
        exit(1) 
    doubleParents = makeTaxTree(goodSplits, contigTax, OUTTREE)
    # assert not doubleParents # There shouldn't be any, but comment for being lenient. The makeTaxTree function warns the user if double parents are found.

    ### Run anvi'o pipeline
    call(['anvi-delete-collection', '-p', args.profile_db, '-C', COLLECTION_NAME], stdout=DEVNULL, stderr=DEVNULL)

    with open('good_splits.tsv', 'w') as outfile:
        for split in goodSplits:
            outfile.write('{}\t{}\n'.format(split, COLLECTION_NAME))

    if args.split_mode == 'safe' or not samples: # Yolo parser expects non-blank profiles.
        cdb = args.contigs_db
        pdb = args.profile_db

    else: # yolo!
        try:
            subset_anvio(goodSplits, args.contigs_db, args.profile_db, outdir+'_YOLO')
        except sqlite3.OperationalError:
            stderr.write('\nSomething went wrong while using the yolo parser. Most likely SQLite\'s SQLITE_MAX_VARIABLE_NUMBER is too low. Please try again adding the "-s safe" parameter.\n\n')
            raise
        cdb = outdir + '_YOLO/CONTIGS.db'
        pdb = outdir + '_YOLO/PROFILE.db'

    call(['anvi-import-collection', '-c', cdb, '-p', pdb, '-C', COLLECTION_NAME, 'good_splits.tsv'])

    splitCommand = ['anvi-split', '-c', cdb, '-p', pdb, '-C', COLLECTION_NAME, '-o', outdir]
    if args.enforce_clustering:
        splitCommand.append('--enforce-hierarchical-clustering')
    else:
        splitCommand.append('--skip-hierarchical-clustering')
    call(splitCommand)
    
    if args.split_mode == 'yolo' and samples:
      call(['rm', '-r', outdir + '_YOLO']) # leave no evidence
 
    if not samples: # Anvi-split creates no profile file if called with a blank profile. We re-generate it.
        profileCommand = ['anvi-profile', '-o', '{}/{}/BLANK'.format(outdir, COLLECTION_NAME), '-c', '{}/{}/CONTIGS.db'.format(outdir, COLLECTION_NAME),
                          '-S','Blank', '--blank-profile']
        if args.enforce_clustering:
            if args.enforce_clustering:
                profileCommand.append('--enforce-hierarchical-clustering')
            else:
                profileCommand.append('--skip-hierarchical-clustering')
        call(profileCommand)
        blankFiles = glob.glob('{}/{}/BLANK/*'.format(outdir, COLLECTION_NAME))
        call(['mv'] + blankFiles + ['{}/{}/'.format(outdir, COLLECTION_NAME)])
        call(['rmdir', '{}/{}/BLANK'.format(outdir, COLLECTION_NAME)])

    
    # Move files before calling anvi-interactive
    call(['mv', 'good_splits.tsv', outdir])
    call(['mv', OUTTREE, outdir])

    ### Run anvi interactive

    interactiveCommand = ['anvi-interactive', '-c', '{}/{}/CONTIGS.db'.format(outdir, COLLECTION_NAME), '-p', '{}/{}/PROFILE.db'.format(outdir, COLLECTION_NAME)]
    if samples:
         interactiveCommand.extend(['-t', '{}/{}'.format(outdir, OUTTREE)])
    else:
        print('anvi-interactive fails to display the custom taxonomy tree when no SAM files are present, so we won\'t include it')

    if args.extra_anvio_args:
        interactiveCommand.extend([x.strip() for x in args.extra_anvio_args.split(' ')])

    print('Running anvi-interactive with:')
    print(' '.join(interactiveCommand))

    call(interactiveCommand)


def makeTaxTree(splits, contigTax, outname):
    RANK_PREFIXES = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    ### Create namespace and node collection
    names = set()
    contigTax2 = {}
    # Consider only contigs with at least one split fultilling the requested conditions.
    contigsWithSplits = {split.rsplit('_split', 1)[0] for split in splits}
    contigTax = {k:v for k,v in contigTax.items() if k in contigsWithSplits} 
    for k, v in contigTax.items():
        for i, p in enumerate(RANK_PREFIXES): # Add rank prefixes
            v[i] = '{}_{}'.format(p, v[i])
            v[i] = v[i].replace('(','[').replace(')',']') # Parentheses will break newick
        contigTax2[k] = v
    [names.update([k]+v) for k,v in contigTax2.items()]
    names.update(splits) # We want to have the contigs AND the splits in our tree
    nodes = {name: dendropy.Node() for name in names}
    taxa = []
    for name, node in nodes.items():
        taxon = dendropy.Taxon(name)
        node.taxon = taxon
        taxa.append(taxon)
    namespace = dendropy.TaxonNamespace()
    namespace.add_taxa(taxa)

    ### Create and populate tree
    tree = dendropy.Tree(taxon_namespace=namespace)
    parents={}

    for split in splits:
        contig = split.rsplit('_split', 1)[0]
        tax = contigTax2[contig]
        tree.seed_node.add_child(nodes[tax[0]])
        for i in range(1,len(tax)):
            nodes[tax[i-1]].add_child(nodes[tax[i]])
            if tax[i] not in parents:
                parents[tax[i]] = set([tax[i-1]])
            else:
                parents[tax[i]].add(tax[i-1])
              
        nodes[tax[-1]].add_child(nodes[contig])
        nodes[contig].add_child(nodes[split])

    # All nodes should have only one parent!
    doubleParents = False
    for p in parents:
        if len(parents[p])>1:
            if not doubleParents:
                doubleParents = True
                print('Some nodes in your taxonomy have more than one parent!!')
                print('This is due to a bug that should have been corrected some time ago')
                print('Maybe you ran the project with an older version of SqueezeMeta?')
                print('You will be unable to load the taxonomy tree in anvi\'o for this particular query')
            #print(p, parents[p])
  

    with open(outname, 'w') as outfile:
        outfile.write(tree.as_string('newick').replace('\'',''))

    return doubleParents


def parse_args():
    parser = argparse.ArgumentParser(description='Filter SqueezeMeta results and visualize them with anvi\'o', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('-p', '--profile-db', type=str, required=True, help='Anvi\'o profile database')
    parser.add_argument('-c', '--contigs-db', type=str, required=True, help='Anvi\'o contigs database')
    parser.add_argument('-t', '--taxonomy', type=str, required = True, help='SqueezeMeta contigs taxonomy')
    parser.add_argument('-q', '--query', type=str, required=True, nargs='+')
    parser.add_argument('-o', '--output-dir', type=str, default='filteredDB', help='Output directory')
    parser.add_argument('-m', '--max-splits', type=int, default=25000, help='Maximum number of splits to visualize')
    parser.add_argument('--enforce-clustering', action='store_true', help='Hierarchically cluster splits based on abundances and composition')
    parser.add_argument('--extra-anvio-args', type=str, help='Extra arguments for anvi-interactive')
    parser.add_argument('-s', '--split-mode', type=str, choices = ['safe', 'yolo'], default='yolo', help='How should we try to subset the database?')
    args = parser.parse_args()
    args.query = ' '.join(args.query)
    return args
 


if __name__=='__main__':
    main(parse_args())
