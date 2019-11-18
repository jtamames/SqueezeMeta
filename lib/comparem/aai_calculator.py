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

import os
import time
import logging
from collections import defaultdict

from numpy import mean, std

import biolib.seq_io as seq_io
from biolib.parallel import Parallel
from biolib.common import make_sure_path_exists, remove_extension, concatenate_files


class AAICalculator(object):
    """Calculate AAI between all pairs of genomes."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus

    def _genome_offsets(self, sorted_hit_table):
        """Read blast table to determine byte offsets of hits for each genome.

        Parameters
        ----------
        sorted_hit_table : str
            File containing sorted blast hits.

        Returns
        -------
        dict : d[genome_id] -> (start_pos, end_pos)
           Start and end byte offsets of hits for each genome in blast table.
        """

        offset_table = defaultdict(dict)
        # with open(sorted_hit_table, 'rb', 512 * (10 ** 6)) as f:
        with open(sorted_hit_table) as f: # Open in read mode with default buffering, FPS, NOV 25 2019
            cur_query_genome = None
            cur_target_genome = None
            start_pos = 0
            end_pos = 0
            for line in f:
                hit = line.split('\t')
                query_genome = hit[0]
                target_genome = hit[2]

                if target_genome != cur_target_genome or query_genome != cur_query_genome:
                    if cur_query_genome:
                        offset_table[cur_query_genome][cur_target_genome] = (start_pos, end_pos)

                    cur_query_genome = query_genome
                    cur_target_genome = target_genome
                    start_pos = end_pos

                end_pos += len(line)

            offset_table[cur_query_genome][cur_target_genome]  = (start_pos, end_pos)

        return offset_table

    def _valid_hits(self, hit_table_stream, 
                            offset_table,
                            evalue_threshold,
                            per_identity_threshold, 
                            per_aln_len_threshold,
                            query_genome_id, 
                            target_genome_id):
        """Identify best hits from a genome meeting the specified criteria.

        Hits from genes within query genome are identified which
        satisfy the percent identity threshold, percent alignment
        length threshold, and are to the specified target genome.
        For each gene, the hit with the highest bitscore is identified.

        Parameters
        ----------
        hit_table_stream : stream
            Stream to table with blast hits.
        offset_table : d[genome_id] -> (start_pos, end_pos)
           Start and end byte offsets of hits for each genome in blast table.
        evalue_threshold : float
            Evalue threshold used to define a homologous gene.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        query_genome_id : str
            Unique id of genome to obtained hits for.
        target_genome_id : str
            Unique id of genome to considered hits to.
            
        Returns
        -------
        dict : d[query_id] -> list of lists with blast hit information
           Hits from query genome meeting specified criteria (can be multiple hits if equal bitscores).
        """

        # get valid hits for genome
        hits = {}

        if query_genome_id not in offset_table:
            # proteins in query genome failed to hit any target proteins
            return hits
            
        if target_genome_id not in offset_table[query_genome_id]:
            # proteins in query genome failed to hit any proteins in target genome
            return hits
  
        start_pos, end_pos = offset_table[query_genome_id][target_genome_id]
        hit_table_stream.seek(start_pos)
        while hit_table_stream.tell() < end_pos:
            hit = hit_table_stream.readline().split('\t')

            perc_iden = float(hit[4])
            if perc_iden < per_identity_threshold:
                continue
                
            evalue = float(hit[12])

            query_id = hit[0] + '~' + hit[1]
            query_coverage = int(hit[9]) - int(hit[8]) + 1
            per_aln_len = query_coverage * 100.0 / self.gene_lengths[query_id]
            if per_aln_len < per_aln_len_threshold:
                continue

            target_genome = hit[2]
            target_id = target_genome + '~' + hit[3]
            if target_genome_id and target_genome != target_genome_id:
                continue

            bitscore = float(hit[13])

            prev_hit = hits.get(query_id, None)
            if not prev_hit:
                hits[query_id] = [[target_id, perc_iden, per_aln_len, evalue, bitscore]]
            elif bitscore > prev_hit[0][4]:
                # for each gene, keep the hit with the highest bitscore
                hits[query_id] = [[target_id, perc_iden, per_aln_len, evalue, bitscore]]
            elif bitscore == prev_hit[0][4]:
                hits[query_id].append([target_id, perc_iden, per_aln_len, evalue, bitscore])

        return hits

    def _producer(self, genome_id_list):
        """Identify reciprocal best blast hits between pairs of genomes.

        Parameters
        ----------
        genome_info_pairs : ((genome_idA, # genes), (genome_idB, # genes))
            Identifier of genomes to process.
        """

        #hit_table_stream = open(self.sorted_hit_table, 'rb', 128 * (10 ** 6))
        hit_table_stream = open(self.sorted_hit_table) # Open in read mode with default buffering, FPS, NOV 25 2019 

        # get genome ID and number of genes in genomes to process
        query_genome_id, genomes_to_process = genome_id_list
        
        if self.keep_rbhs:
            fout_stats = open(os.path.join(self.output_dir, query_genome_id + '.rbh.tsv'), 'w')
            fout_stats.write('Genome A\tGenome B\tPercent Identity\tPercent Alignment Length\te-value\tbitscore\n')

        # determing RBHs 
        results = []           
        for cur_genome_id in genomes_to_process:
            hits = self._valid_hits(hit_table_stream, 
                            self.offset_table,
                            self.evalue_threshold,
                            self.per_identity_threshold, 
                            self.per_aln_len_threshold,
                            query_genome_id, 
                            cur_genome_id)
                                    
            cur_hits = self._valid_hits(hit_table_stream, 
                                        self.offset_table,
                                        self.evalue_threshold,
                                        self.per_identity_threshold, 
                                        self.per_aln_len_threshold,
                                        cur_genome_id, 
                                        query_genome_id)

            # report reciprocal best blast hits
            per_identity_hits = []
            for query_id, hit_stats in hits.items():
                bRBH = False
                for query_hit in hit_stats:
                    target_id, per_identA, per_aln_lenA, evalueA, bitscoreA = query_hit
                    
                    for target_hit in cur_hits.get(target_id, []):
                        cur_target_id, per_identB, per_aln_lenB, evalueB, bitscoreB = target_hit
                        if query_id != cur_target_id:
                            continue

                        # take average of statistics in both blast directions as
                        # the results will be similar, but not identical
                        per_ident = 0.5 * (per_identA + per_identB)
                        per_identity_hits.append(per_ident)

                        per_aln_len = 0.5 * (per_aln_lenA + per_aln_lenB)
                        evalue = 0.5 * (evalueA + evalueB)
                        bitscore = 0.5 * (bitscoreA + bitscoreB)

                        if self.keep_rbhs: 
                            fout_stats.write('%s\t%s\t%.2f\t%.2f\t%.2g\t%.2f\n' % (query_id, target_id, per_ident, per_aln_len, evalue, bitscore))

                        # keep only one reciprocal best hit per gene
                        bRBH = True
                        break
                 
                    if bRBH: # keep only one reciprocal best hit per gene
                        break
                                  
            mean_per_identity_hits = 0
            if len(per_identity_hits) > 0:
                mean_per_identity_hits = mean(per_identity_hits)

            std_per_identity_hits = 0
            if len(per_identity_hits) >= 2:
                std_per_identity_hits = std(per_identity_hits)

            num_genesA = self.query_gene_count[query_genome_id]
            num_genesB = self.target_gene_count[cur_genome_id]
            num_rbhs = len(per_identity_hits)
            of = num_rbhs * 100.0 / min(num_genesA, num_genesB)
            
            results.append((query_genome_id,
                                num_genesA,
                                cur_genome_id,
                                num_genesB,
                                num_rbhs,
                                mean_per_identity_hits,
                                std_per_identity_hits,
                                of))
                    
        
        if self.keep_rbhs:
            fout_stats.close()
            
        return results

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : tuple
            Summary statistics for a genome pair.
        consumer_data : list
            Summary statistics of amino acid identity between genome pairs.

        Returns
        -------
        consumer_data
            Summary statistics of amino acid identity between genome pairs.
        """

        if consumer_data == None:
            # setup structure for consumed data
            consumer_data = []

        self.processed_paired += len(produced_data)
        consumer_data.extend(produced_data)

        return consumer_data

    def _progress(self, processed_genomes, total_genomes):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_genomes : int
            Number of genomes processed.
        total_genomes : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return '  Finished processing %d of %d (%.2f%%) pairs.' % (self.processed_paired, 
                                                                    self.num_pairs, 
                                                                    float(self.processed_paired) * 100 / self.num_pairs)

    def run(self, query_gene_file,
                    target_gene_file,
                    sorted_hit_table, 
                    evalue_threshold, 
                    per_iden_threshold, 
                    per_aln_len_threshold,
                    keep_rbhs,
                    output_dir):
        """Calculate amino acid identity (AAI) between pairs of genomes.

        Parameters
        ----------
        query_gene_file : str
            File with all query genes in FASTA format.
        target_gene_file : str or None
            File with all target genes in FASTA format, or None if performing a reciprocal AAI calculation.
        sorted_hit_table : str
            Sorted table indicating genes with sequence similarity.
        evalue_threshold : float
            Evalue threshold used to define a homologous gene.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        keep_rbhs : boolean
            Flag indicating if RBH should be written to file.
        output_dir : str
            Directory to store AAI results.
        """

        self.sorted_hit_table = sorted_hit_table
        self.evalue_threshold = evalue_threshold
        self.per_identity_threshold = per_iden_threshold
        self.per_aln_len_threshold = per_aln_len_threshold
        self.keep_rbhs = keep_rbhs
        self.output_dir = output_dir

        # calculate length of genes and number of genes in each genome
        self.logger.info('Calculating length of genes.')
        self.gene_lengths = {}
        self.query_gene_count = defaultdict(int)
        query_genomes = set()
        for seq_id, seq in seq_io.read_fasta_seq(query_gene_file):
            if seq[-1] == '*':
                self.gene_lengths[seq_id] = len(seq) - 1
            else:
                self.gene_lengths[seq_id] = len(seq)
                
            genome_id = seq_id[0:seq_id.find('~')]
            self.query_gene_count[genome_id] += 1
            query_genomes.add(genome_id)
            
        self.target_gene_count = defaultdict(int)
        target_genomes = set()
        if target_gene_file:
            for seq_id, seq in seq_io.read_fasta_seq(target_gene_file):
                if seq[-1] == '*':
                    self.gene_lengths[seq_id] = len(seq) - 1
                else:
                    self.gene_lengths[seq_id] = len(seq)
                    
                genome_id = seq_id[0:seq_id.find('~')]
                self.target_gene_count[genome_id] += 1
                target_genomes.add(genome_id)
        else:
            self.target_gene_count = self.query_gene_count

        # get byte offset of hits from each genome
        self.logger.info('Indexing sorted hit table.')
        self.offset_table = self._genome_offsets(self.sorted_hit_table)

        # calculate AAI between each pair of genomes in parallel
        if target_genomes:
            # compare query genomes to target genomes
            self.num_pairs = len(query_genomes) * len(target_genomes)
            self.logger.info('Calculating AAI between %d query and %d target genomes:' % (len(query_genomes), len(target_genomes)))
        else:
            # compute pairwise values between target genomes
            ng = len(query_genomes)
            self.num_pairs = (ng*ng - ng) / 2
            self.logger.info('Calculating AAI between all %d pairs of genomes:' % self.num_pairs)
            
        if self.num_pairs == 0:
            self.logger.warning('No genome pairs identified.')
            return

        genome_id_lists = []
        query_genomes = list(query_genomes)
        target_genomes = list(target_genomes)
        for i in range(0, len(query_genomes)):
            genome_idI = query_genomes[i]
            
            if target_genomes:
                genome_id_list = target_genomes
            else:
                genome_id_list = []
                for j in range(i + 1, len(query_genomes)):
                    genome_idJ = query_genomes[j]
                    genome_id_list.append(genome_idJ)

            genome_id_lists.append((genome_idI, genome_id_list))

        self.processed_paired = 0
        parallel = Parallel(self.cpus)
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None
        consumer_data = parallel.run(self._producer, self._consumer, genome_id_lists, progress_func)

        # write results for each genome pair
        self.logger.info('Summarizing AAI results.')
        aai_summay_file = os.path.join(output_dir, 'aai_summary.tsv')
        fout = open(aai_summay_file, 'w')
        fout.write('Genome A\tGenes in A\tGenome B\tGenes in B\t# orthologous genes\tMean AAI\tStd AAI\tOrthologous fraction (OF)\n')

        for data in consumer_data:
            fout.write('%s\t%d\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f\n' % data)

        fout.close()

        # concatenate RBH files
        rbh_output_file = None
        if self.keep_rbhs:
            self.logger.info('Concatenating RBH files.')
            rbh_files = []
            for genome_id in query_genomes:
                rbh_files.append(os.path.join(self.output_dir, genome_id + '.rbh.tsv'))
                
            rbh_output_file = os.path.join(self.output_dir, 'rbh.tsv')
            concatenate_files(rbh_files, rbh_output_file, common_header=True)
            
            for f in rbh_files:
                os.remove(f)
                
        return aai_summay_file, rbh_output_file
