###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import random

import dendropy

from biolib.newick import parse_label, create_label

"""Perform non-parametric bootstrapping on multiple sequence alignment."""

def bootstrap_support(input_tree, replicate_trees, output_tree):
    """ Calculate support for tree with replicates covering the same taxon set.

    Parameters
    ----------
    input_tree : str
      Tree inferred from complete data.
    replicate_trees : iterable
      Files containing replicate trees.
    output_tree: str
      Name of output tree with support values.
    """
    # read tree and bootstrap replicates as unrooted, and
    # calculate bootstrap support
    orig_tree = dendropy.Tree.get_from_path(input_tree, 
                                            schema='newick', 
                                            rooting="force-unrooted", 
                                            preserve_underscores=True)
    orig_tree.bipartitions = True
    orig_tree.encode_bipartitions()

    rep_trees = dendropy.TreeArray(taxon_namespace=orig_tree.taxon_namespace,
                                    is_rooted_trees=False,
                                    ignore_edge_lengths=True,
                                    ignore_node_ages=True,
                                    use_tree_weights=False)

    rep_trees.read_from_files(files=replicate_trees,
                                schema='newick',
                                rooting="force-unrooted",
                                preserve_underscores=True,
                                taxon_namespace=orig_tree.taxon_namespace)

    rep_trees.summarize_splits_on_tree(orig_tree,
                                       is_bipartitions_updated=True,
                                       add_support_as_node_attribute=True,
                                       support_as_percentages=True)

    for node in orig_tree.internal_nodes():
        if node.label:
            support, taxon, aux_info = parse_label(node.label)
            node.label = create_label(node.support, taxon, aux_info)
        else:
            node.label = str(int(node.support))

    orig_tree.write_to_path(output_tree, 
                            schema='newick', 
                            suppress_rooting=True, 
                            unquoted_underscores=True)


def bootstrap_alignment(msa, output_file, frac=1.0):
    """Bootstrap multiple sequence alignment.

    True bootstrapping requires subsampling an alignment,
    with replacement, to construct new alignments
    with the same length as the input alignment. The
    'frac' parameter allows shorter bootstrap alignments
    to be generated in order to reduce computational
    demands.

    Parameters
    ----------
    msa : d[seq_id] -> seq
      Full multiple sequence alignment.
    output_file : str
      File to write bootstrapped alignment.
    frac : float
      Fraction of alignment to subsample.
    """
    alignment_len = len(msa[list(msa.keys())[0]])
    sample_len = int(alignment_len * frac)
    cols = [random.randint(0, alignment_len - 1) for _ in range(sample_len)]

    fout = open(output_file, 'w')
    for seq_id, seq in msa.items():
        fout.write('>' + seq_id + '\n')
        for col in cols:
            fout.write(seq[col])
        fout.write('\n')
    fout.close()
