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
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import logging
import tempfile
from collections import defaultdict

from comparem.aai_calculator import AAICalculator

from biolib.common import concatenate_files, make_sure_path_exists


class Classify(object):
    """Classify genomes based on AAI to reference genomes."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def run(self, query_gene_file,
                    target_gene_file,
                    sorted_hit_table, 
                    evalue_threshold, 
                    per_iden_threshold, 
                    per_aln_len_threshold,
                    num_top_targets,
                    taxonomy_file,
                    keep_rbhs,
                    output_dir):
        """Classify genomes based on AAI to reference genomes.

        Parameters
        ----------
        query_gene_file : str
            File with all query genes in FASTA format.
        target_gene_file : str
            File with all target genes in FASTA format.
        sorted_hit_table : str
            Sorted table indicating genes with sequence similarity.
        evalue_threshold : float
            Evalue threshold used to define a homologous gene.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        num_top_targets : int
            Number of top scoring target genomes to report per query genome.
        taxonomy_file : str
            File indicating taxonomic identification of all target genomes.
        keep_rbhs : boolean
            Flag indicating if RBH should be written to file.
        output_dir : str
            Directory to store AAI results.
        """
        
        # read taxonomic identification of each genome
        taxonomy = {}
        if taxonomy_file:
            for line in open(taxonomy_file):
                genome_id, taxa_str = line.rstrip().split('\t')
                taxonomy[genome_id] = taxa_str

        # calculate AAI between query and target genomes
        aai_output_dir = os.path.join(output_dir, 'aai')
        make_sure_path_exists(aai_output_dir)
        aai_calculator = AAICalculator(self.cpus)
        aai_output_file, rbh_output_file = aai_calculator.run(query_gene_file,
                                                                target_gene_file,
                                                                sorted_hit_table,
                                                                evalue_threshold,
                                                                per_iden_threshold,
                                                                per_aln_len_threshold,
                                                                keep_rbhs,
                                                                aai_output_dir)

        # determine matches to each query genomes
        aai_results_file = os.path.join(aai_output_dir, 'aai_summary.tsv')
        with open(aai_results_file) as f:
            f.readline()
            
            hits = defaultdict(list)
            for line in f:
                line_split = line.rstrip().split('\t')
                query_id = line_split[0]
                target_id = line_split[2]
                aai = float(line_split[5])
                of = float(line_split[7])
                
                hits[query_id].append([target_id, aai, of])
                
        # report top matches
        results_file = os.path.join(output_dir, 'classify.tsv')
        fout = open(results_file, 'w')
        fout.write('Query Id\tTarget Id\tAAI\tOF\tScore')
        if taxonomy:
            fout.write('\tTarget Taxonomy')
        fout.write('\n')
             
        for query_id, cur_hits in hits.items():
            cur_hits.sort(key=lambda x: x[1], reverse=True)
            for i in range(0, min(num_top_targets, len(cur_hits))):
                data = [query_id] + cur_hits[i]
                fout.write('%s\t%s\t%.2f\t%.2f' % tuple(data))
                
                aai = data[2]
                of = data[3]
                fout.write('\t%.2f' % (aai+of))
                
                target_id = cur_hits[i][0]
                if target_id in taxonomy:
                    fout.write('\t%s' % taxonomy[target_id])
                
                fout.write('\n')
        fout.close()
        
        return results_file
