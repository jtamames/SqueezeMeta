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
import logging
import ntpath
import operator
from collections import defaultdict

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.parallel import Parallel

from numpy import mean, std


class LgtCodon(object):
    """Calculate codon usage of genes in genomes."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def _manhattan(self, gene_codon_usage, genome_codon_usage):
        """Calculate Manhattan distance.

        Parameters
        ----------
        gene_codon_usage : dict[gene_id][codon] -> count
            Occurrence of codon in each gene.
        genome_codon_usage : dict[gene_id][codon] -> count
            Occurrence of codon across the genome.

        Returns
        -------
        dict : d[gene_id] ->  [distance, std. dev. from mean]
             Manhattan distance and standard deviation from the mean for each gene.
        """

        dist = {}
        genome_sum_codon = sum(genome_codon_usage.values())
        for gene_id, codon in gene_codon_usage.items():
            d = 0
            gene_sum_codon = sum(codon.values())
            for di in genome_codon_usage:
                d += abs(codon.get(di, 0) * 100.0 / gene_sum_codon - genome_codon_usage.get(di, 0) * 100.0 / genome_sum_codon)
            dist[gene_id] = [d]

        # model all distances as a normal distribution
        m = mean(list(dist.values()))
        s = std(list(dist.values()))

        # calculate standard deviations from the mean
        for gene_id, d in dist.items():
            dist[gene_id].append((d - m) / s)

        return dist

    def codon_usage(self, seqs, genome_id):
        """Calculate codon usage statistics for genome.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        genome_id : str
            Unique id of genome used to create output file.
        """

        # calculate codon usage for each gene and the genome as a while
        gene_codon_usage = defaultdict(lambda: defaultdict(int))
        genome_codon_usage = defaultdict(int)
        gc = {}
        for gene_id, seq in seqs.items():
            gc[gene_id] = seq_tk.gc(seq)

            for i in range(0, len(seq) - 3, 3):
                codon = seq[i:i + 3].upper()
                if 'N' not in codon:
                    gene_codon_usage[gene_id][codon] += 1
                    genome_codon_usage[codon] += 1

        # calculate Manhattan distance for each gene
        manhattan_dist = self._manhattan(gene_codon_usage, genome_codon_usage)

        # report dinucleotide usage of each gene
        codon_set_sorted = sorted(genome_codon_usage.keys())
        gene_ids_sorted = sorted(list(manhattan_dist.items()), key=operator.itemgetter(1), reverse=True)

        output_file = os.path.join(self.output_dir, genome_id + '.codon_usage.tsv')
        fout = open(output_file, 'w')

        fout.write('Gene Id\tGC\tLength (bp)\t# codons\tManhattan distance\tDeviations from mean')
        for di in codon_set_sorted:
            fout.write('\t' + di)
        fout.write('\n')

        genome_gc = seq_tk.gc(''.join(list(seqs.values())))
        genome_sum_codon = sum(genome_codon_usage.values())
        fout.write('%s\t%.2f\t%d\t%d' % ('<complete genome>', genome_gc * 100.0, sum([len(x) for x in list(seqs.values())]), genome_sum_codon))
        fout.write('\t%.1f\t%.1f' % (0, 0))
        for di in codon_set_sorted:
            fout.write('\t%.2f' % (genome_codon_usage.get(di, 0) * 100.0 / genome_sum_codon))
        fout.write('\n')

        for gene_id, dists in gene_ids_sorted:
            codons = gene_codon_usage[gene_id]
            sum_codons = sum(codons.values())
            fout.write('%s\t%.2f\t%d\t%d' % (gene_id, gc[gene_id] * 100, len(seqs[gene_id]), sum_codons))
            fout.write('\t%.2f\t%.2f' % (dists[0], dists[1]))

            for di in codon_set_sorted:
                fout.write('\t%.2f' % (codons.get(di, 0) * 100.0 / sum_codons))
            fout.write('\n')
        fout.close()

    def _producer(self, gene_file):
        """Calculates codon usage statistics of a genome.

        Parameters
        ----------
        gene_file : str
            Fasta file containing amino acid sequences.
        """

        genome_id = ntpath.basename(gene_file)
        genome_id = genome_id.replace('.genes.fna', '')
        genome_id = os.path.splitext(genome_id)[0]

        seqs = seq_io.read_fasta(gene_file)
        self.codon_usage(seqs, genome_id)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of genomes processed.
        total_items : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return '    Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self, gene_files, output_dir):
        """Calculate codon usage over genes with a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes in nucleotide space.
        output_dir : str
            Directory to store results.
        """

        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.logger.info('Calculating codon usage for each genome.')
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None

        parallel = Parallel(self.cpus)
        parallel.run(self._producer, None, gene_files, progress_func)
