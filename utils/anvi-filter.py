#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Select the anvi'o splits (which are themselves contigs, or parts of
long contigs) that fulfill a set of conditions, and visualize them
using anvi-interactive.

USAGE:

anvio-view.py -p PROFILE_DB -c CONTIGS_DB -t CONTIGS_TAXONOMY -q QUERY
     [--max-splits int ] [--enforce-clustering]
     [--extra-anvio-args "extra args"]

QUERY SYNTAX:

- Queries are combinations of relational operations in the form of
  <SUBJECT> <OPERATOR> <VALUE> (e.g. "PHYLUM == Bacteroidetes""
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
    - FUNH: search within the KEGG BRITE functional hierarchy
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

- The "--max-splits" parameter controls the maximum number of splits
  that will be loaded into anvi\'o. If the provided query returns a
  higher number of splits, the program will stop. By default it is set
  to 25,000, larger values may make the anvi\'o interface to respond
  slowly. Setting --max-splits to 0 will allow an arbitrarily large
  number of splits to be loaded.

- By default, splits are only clustered based on their taxonomy. Passing
  the "--enforce-clustering" flag will make anvi\'o perform an
  additional clustering based on abundances across samples and sequence
  composition.

- Extra arguments for anvi-interactive can be passed with the parameter
  --extra-anvio-args "extra args". Extra arguments must be surrounded
  by quotes.
     e.g. --extra-anvio-args "--taxonomic-level t_phylum --title Parrot"
"""

import argparse
from os.path import abspath, dirname, realpath
from sys import path
utils_home = abspath(dirname(realpath(__file__)))
path.append('{}/../lib/'.format(utils_home))
import dendropy
import sqlite3
from collections import defaultdict
from sys import argv
from subprocess import call
from os import devnull
DEVNULL = open(devnull, 'wb')
from splitFilter import SplitFilter

OUTTREE = 'Taxonomy.nwk'
COLLECTION_NAME = 'SqueezeMeta'

def main(args):
    missingAnvi = call(['which', 'anvi-self-test'], stdout=DEVNULL)
    if missingAnvi: # ecode of which was 1
        print('We can\'t find anvi\'o in your PATH. Is the environment activated?')
        exit(-1)
    sfilter = SplitFilter(args.contigs_db, args.profile_db, args.taxonomy)
    splits, samples, contigTax = sfilter.get_parsed_annotations()
    print('')
    print('- {} splits found in anvi\'o profile database.'.format(len(splits)))
    print('')
    print('- Sample names are:')
    print('\t{}'.format('\t'.join(samples)))
    print('')
    print('- Query tree is:')
    sfilter.print_tree(args.query)
    print('')
    goodSplits = sfilter.run(args.query)
    print('- {} splits fulfilling the required conditions.'.format(len(goodSplits)))
    if not goodSplits:
        print ('Exiting')
        exit(0)
    if args.max_splits and len(goodSplits) > args.max_splits:
        print('')
        print('Too many splits satisfying the provided query.')
        print('Loading too many splits can make anvi\'o slow to respond. Please consider providing a narrower query.')
        print('Alternatively, use --max-splits NUMBER to increase this threshold or --max-splits 0 to disable it completely')
        print('')
        exit(0) 
    removedSplits = makeTaxTree(goodSplits, contigTax, OUTTREE)
    goodSplits = goodSplits - removedSplits
    with open('tempCollection.tsv', 'w') as outfile:
        for split in goodSplits:
            outfile.write('{}\t{}\n'.format(split, COLLECTION_NAME))
    call(['rm', '-r', 'tempCol']) # we should do the rm and the delete-collection after anvi-interactive is killed manually, but it doesn't seem to work
    call(['anvi-delete-collection', '-p', args.profile_db, '-C', COLLECTION_NAME])
    call(['anvi-import-collection', '-c', args.contigs_db, '-p', args.profile_db, '-C', COLLECTION_NAME, 'tempCollection.tsv'])

    splitCommand = ['anvi-split', '-c', args.contigs_db, '-p', args.profile_db, '-C', COLLECTION_NAME, '-o', 'tempCol']
    if args.enforce_clustering:
        splitCommand.append('--enforce-hierarchical-clustering')
    else:
        splitCommand.append('--skip-hierarchical-clustering')
    call(splitCommand)

    interactiveCommand = ['anvi-interactive', '-c', 'tempCol/{}/CONTIGS.db'.format(COLLECTION_NAME), '-p', 'tempCol/{}/PROFILE.db'.format(COLLECTION_NAME), '-t', OUTTREE]

    if args.extra_anvio_args:
        interactiveCommand.extend([x.strip() for x in args.extra_anvio_args.split(' ')])

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

    removedSplits = set() # This shouldn't be needed but since we have taxonomy problems do it for now.

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
    for p in parents:
        if len(parents[p])>1:
            print(p, parents[p])

    with open(outname, 'w') as outfile:
        outfile.write(tree.as_string('newick').replace('\'',''))

    return removedSplits


def parse_args():
    parser = argparse.ArgumentParser(description='Filter SqueezeMeta results and visualize them with anvi\'o', epilog='Fernando Puente-Sánchez (CNB) 2019\n')
    parser.add_argument('-p', '--profile-db', type=str, required=True, help='Anvi\'o profile database')
    parser.add_argument('-c', '--contigs-db', type=str, required=True, help='Anvi\'o contigs database')
    parser.add_argument('-t', '--taxonomy', type=str, required = True, help='SqueezeMeta contigs taxonomy')
    parser.add_argument('-q', '--query', type=str, required=True, nargs='+')
    parser.add_argument('-m', '--max-splits', type=int, default=25000, help='Maximum number of splits to visualize')
    parser.add_argument('--enforce-clustering', action='store_true', help='Hierarchically cluster splits based on abundances and composition')
    parser.add_argument('--extra-anvio-args', type=str, help='Extra arguments for anvi-interactive')
    args = parser.parse_args()
    args.query = ' '.join(args.query)
    return args
 


if __name__=='__main__':
    main(parse_args())
