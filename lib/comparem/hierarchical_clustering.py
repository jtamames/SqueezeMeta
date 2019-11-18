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



import sys
import logging
from collections import defaultdict

import numpy as np
from scipy.cluster import hierarchy

class HierarchicalCluster(object):
    """Perform hierarchical clustering on pairwise value files."""
    
    def __init__(self):
        """Initialization."""
        
        self.logger = logging.getLogger('timestamp')

    def _parse_data(self, 
                        pairwise_value_file, 
                        name_col1, 
                        name_col2, 
                        value_col,
                        similarity,
                        max_sim_value):
        """Initialization.

        Parameters
        ----------
        pairwise_value_file : str
            File with pairwise similarity or dissimilarity values.
        name_col1 : int
            Index of first column with genome names.
        name_col2 : int
            Index of second column with genome names.
        value_col : int
            Index of column with similarity or dissimilarity values.
        similarity : boolean
            Flag indicating file contain similarity values.
        max_sim_value : float   
            Maximum value of similarity scores.
            
        Returns
        -------
        TBD
        """
        
        upper_triangle = defaultdict(lambda : defaultdict(float))
        genome_ids = set()
        with open(pairwise_value_file) as f:
            f.readline()

            genomes = set()
            for line in f:
                line_split = line.rstrip().split('\t')
                
                genome_idA = line_split[name_col1]
                genome_idB = line_split[name_col2]
                value = float(line_split[value_col])
                     
                if similarity:
                    if max_sim_value < value:
                        self.logger.error('Maximum similarity score %f is less than identified pairwise value %f.' % (max_sim_value, value))
                        sys.exit(-1)
                    value = max_sim_value - value
                    
                upper_triangle[genome_idA][genome_idB] = value
                genome_ids.add(genome_idA)
                genome_ids.add(genome_idB)
                
        # sort by row size
        genome_labels = []  
        for genome_id in sorted(upper_triangle, key=lambda genome_id: len(upper_triangle[genome_id]), reverse=True):
            genome_labels.append(genome_id)
            genome_ids.remove(genome_id)
            
        # add in last row with no entries
        if len(genome_ids) != 1:
            self.logger.error("Pairwise value file appears is incorrectly formatted.")
            sys.exit(-1)
        genome_labels.append(genome_ids.pop())
 
        diss_vector = []
        for i in range(0, len(genome_labels)):
            genome_idA = genome_labels[i]
            for j in range(i+1, len(genome_labels)):
                genome_idB = genome_labels[j]
                diss_vector.append(upper_triangle[genome_idA][genome_idB])
        
        return diss_vector, genome_labels
                
    def _save_newick(self, node, newick, parentdist, leaf_names):
        if node.is_leaf():
            genome_id = leaf_names[node.id]
            return "%s:%.2f%s" % (genome_id, parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = self._save_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = self._save_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick
        
    def run(self, pairwise_value_file,
                    method, 
                    similarity,
                    max_sim_value,
                    name_col1,
                    name_col2,
                    value_col,
                    output_tree):
        """Perform hierarchical clustering on pairwise value files.

        Parameters
        ----------
        pairwise_value_file : str
            File with pairwise similarity or dissimilarity values.
        method : str
            Clustering method to use.
        similarity : boolean
            Flag indicating file contain similarity values.
        max_sim_value : float   
            Maximum value of similarity scores.
        name_col1 : int
            Index of first column with genome names.
        name_col2 : int
            Index of second column with genome names.
        value_col : int
            Index of column with similarity or dissimilarity values.
        """
        
        diss_vector, genome_labels = self._parse_data(pairwise_value_file, 
                                                        name_col1, 
                                                        name_col2, 
                                                        value_col, 
                                                        similarity, 
                                                        max_sim_value)
        
        clusters = hierarchy.linkage(diss_vector, method=method)

        tree = hierarchy.to_tree(clusters)
        newick_str = self._save_newick(tree, "", tree.dist, genome_labels)
        
        fout = open(output_tree, 'w')
        fout.write(newick_str + '\n')
        fout.close()
        
      