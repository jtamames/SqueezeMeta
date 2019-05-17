#!/usr/bin/python

"""
Select the anvi'o splits (which are themselves contigs, or parts of long contigs) that fulfill a set of conditions, and visualize them using anvi-interactive.

USAGE:
anvio-view.py "<QUERY>"

QUERY SYNTAX:
- Queries are combinations of relational operations in the form of <SUBJECT> <OPERATOR> <VALUE> (e.g. "PHYLUM == Bacteroidetes"" joined by logical operators (AND, OR).
- Parentheses can be used to group operations together.
- The "AND" and "OR" logical operators can't appear together into the same expression. Parentheses must be used to separate them into different expressions.
    - eg. "GENUS == Escherichia OR GENUS == Prevotella AND FUN CONTAINS iron" would not be valid. Parentheses must be used to write either:
        - "(GENUS == Escherichia OR GENUS == Prevotella)" AND FUN CONTAINS iron" -> splits from either Escherichia or Prevotella which contain ORFs related to iron.
        - "GENUS == Escherichia OR (GENUS == Prevotella AND FUN CONTAINS iron)" -> all splits from Escherichia, and splits of Prevotella which contains ORFs related to iron.

- Another example query would be: "(PHYLUM == Bacteroidetes OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria]) AND FUN CONTAINS iron AND Sample1 > 1"
    - This would select all the anvi'o splits assigned to either the Bacteroidetes phylum or the Alphaproteobacteria or Gammaproteobacteria classes,
      that also contain the substring "iron" in the functional annotations of any of their ORFs,
      and whose anvi'o abundance (mean coverage of a split divided by overall sample mean coverage) in Sample1 is higher than 1.

- Possible subjects are:
    - FUN: search within all the combined databases used for functional annotation.
    - SUPERKINGDOM, PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES: search within the taxonomic annotation at the requested taxonomic rank.
    - <SAMPLE_NAME>: search within the abundances in the sample named <SAMPLE_NAME>
        (e.g. if you have two samples named "Sample1" and "Sample2", the query string "Sample1 > 0.5 AND Sample2 > 0.5" would return the anvi'o splits with an anvi'o abundance higher than 0.5 in both samples)

- Posible relational operators are "==", "!=", ">=", "<=", ">", "<", "IN", "NOT IN", "CONTAINS", "DOES NOT CONTAIN"
"""


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

PROJECT_PATH = '/media/disk8/fer/Hadza'
PROJECT_NAME = 'Hadza3'

#QUERY   = '(PHYLUM IN [Bacteroidetes] AND CLASS!=Cytophagia) OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria] OR (HADZA3_SRR1927149>0.5 AND HADZA3_SRR1927149 < 1)'
#QUERY   = '((PHYLUM IN [Bacteroidetes] AND CLASS!=Cytophagia) OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria]) AND (HADZA3_SRR1927149>0.5 AND HADZA3_SRR1927149 < 1)'
#QUERY   = 'FUN CONTAINS iron OR FUN CONTAINS nitrogen'
#INFILE  = '/home/fer/anviotest/iron.txt'
TAXFILE = '{}/tables/{}.contig.tax.tsv'.format(PROJECT_PATH, PROJECT_NAME)
OUTTREE = '{}/temp.nwk'.format(PROJECT_PATH)
CONTIGS = '{}/contigs.db'.format(PROJECT_PATH)
PROFILE = '{}/SAMPLES-MERGED/PROFILE.db'.format(PROJECT_PATH)

def main():
    QUERY = argv[1]
    assert QUERY
    missingAnvi = call(['which', 'anvi-self-test'], stdout=DEVNULL)
    if missingAnvi: # ecode of which was 1
        print('We can\'t find anvi\'o in your PATH. Is the environment activated?')
        exit(-1)
    sfilter = SplitFilter(CONTIGS, PROFILE, TAXFILE)
    splits, samples, contigTax = sfilter.get_parsed_annotations()
    print('- {} splits found in anvi\'o database'.format(len(splits)))
    print('')
    print('- Sample names are:')
    print('\t{}'.format('\t'.join(samples)))
    print('')
    print('- Query tree is:')
    sfilter.print_tree(QUERY)
    print('')
    goodSplits = sfilter.run(QUERY)
    print('- {} splits fulfilling the required conditions'.format(len(goodSplits)))
    if not goodSplits:
        print ('Exiting')
        exit(0)
    removedSplits = makeTaxTree(goodSplits, contigTax, OUTTREE)
    goodSplits = goodSplits - removedSplits
    with open('tempCollection.tsv', 'w') as outfile:
        for split in goodSplits:
            outfile.write('{}\ttempCol\n'.format(split))
    call(['rm', '-r', 'tempCol']) # we should do the rm and the delete-collection after anvi-interactive is killed manually, but it doesn't seem to work
    call(['anvi-delete-collection', '-p', PROFILE, '-C', 'tempCol'])
    call(['anvi-import-collection', '-c', CONTIGS, '-p', PROFILE, '-C', 'tempCol', 'tempCollection.tsv'])
    call(['anvi-split', '-c', CONTIGS, '-p', PROFILE, '-C', 'tempCol', '-o', 'tempCol'])
    call(['anvi-interactive', '-c', 'tempCol/tempCol/CONTIGS.db', '-p', 'tempCol/tempCol/PROFILE.db', '-t', OUTTREE])

def makeTaxTree(splits, contigTax, outname):
    RANK_PREFIXES = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    # Create namespace and node collection
    names = set()
    contigTax2 = {}
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

    # Create and populate tree
    tree = dendropy.Tree(taxon_namespace=namespace)
    parents={}

    removedSplits = set() # This shouldn't be needed but since we have taxonomy problems do it for now.

    for split in splits:
        contig = split.rsplit('_split', 1)[0]
        tax = contigTax2[contig]
        if tax[-1] == 's_Firmicutes bacterium' or tax[4] == 'f_Clostridia bacterium [no family in NCBI]': # Weird taxonomy, find solution, avoid for now!
            print(contig)
            removedSplits.add(split)
            continue
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



if __name__=='__main__':
    main()
