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
from collections import defaultdict, namedtuple

import biolib.seq_io as seq_io
from biolib.parallel import Parallel

from numpy import mean


class CodonUsage(object):
    """Calculate codon usage over a set of genomes."""

    def __init__(self, cpus=1, keep_ambiguous=False, stop_codon_only=False):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        keep_ambiguous: boolean
            Keep codons with ambiguous bases.
        stop_codon_only: boolean
            Only calculate usage at the last/stop codon position.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus
        self.keep_ambiguous = keep_ambiguous
        self.stop_codon_only = stop_codon_only
        self.stop_codons = ['TAG', 'TAA', 'TGA', 'UAG', 'UAA', 'UGA']

    def codon_usage(self, seqs):
        """ Calculate codon usage within sequences.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.

        Returns
        -------
        dict : dict[codon] -> count
            Occurrence of each codon.
        """

        codon_usage = defaultdict(int)
        gene_length = defaultdict(list)
        for _seq_id, seq in seqs.items():
            if self.stop_codon_only:
                codon = seq[-3:].upper()
                if self.keep_ambiguous or 'N' not in codon:
                    codon_usage[codon] += 1
                    gene_length[codon].append(len(seq))
            else:
                for i in range(0, len(seq) - 3, 3):
                    codon = seq[i:i + 3].upper()
                    if self.keep_ambiguous or 'N' not in codon:
                        codon_usage[codon] += 1

        # get average gene length for each stop codon
        if gene_length:
            for codon, seq_lens in gene_length.items():
                gene_length[codon] = mean(seq_lens)

        return codon_usage, gene_length

    def _producer(self, gene_file):
        """Calculates codon usage of a genome.

        This function is intended to be used as a producer
        within a producer/consumer multiprocessing framework.
        It calculates the codon usage for a single genome
        and returns the results for consumption by the
        consumer function.

        Parameters
        ----------
        gene_file : str
            Fasta file containing nucleotide sequences.

        Returns
        -------
        str
           Unique identifier of genome.
        dict : d[codon] -> count
            Occurrence of each codon.
        dict : d[codon] -> length
            Average length of genes for a given stop codon.
        """

        genome_id = ntpath.basename(gene_file)
        genome_id = genome_id.replace('.genes.fna', '')
        genome_id = os.path.splitext(genome_id)[0]

        seqs = seq_io.read_fasta(gene_file)
        codon_usage, gene_length = self.codon_usage(seqs)

        return (genome_id, codon_usage, gene_length)

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

        This function is intended to be used as a
        consumer within a producer/consumer multiprocessing
        framework. It stores the codon usage for each
        genome into a dictionary and determines the set
        of codons observed across all genomes.

         Parameters
        ----------
        produced_data : list -> [genome_id, codon_usage]
            Unique id of a genome followed by a dictionary
            indicating its codon usage.
        consumer_data : namedtuple
            Set of codons observed across all genomes (codon_set),
            along with the codon usage of each genome (genome_codon_usage),
            and the average length of genes for stop codons (mean_gene_length).

        Returns
        -------
        consumer_data
            The consumer data structure or None must be returned
        """

        if consumer_data == None:
            # setup data to be returned by consumer
            ConsumerData = namedtuple('ConsumerData', 'codon_set genome_codon_usage mean_gene_length')
            consumer_data = ConsumerData(set(), dict(), dict())

        genome_id, codon_usage, mean_gene_length = produced_data

        consumer_data.codon_set.update(list(codon_usage.keys()))
        consumer_data.genome_codon_usage[genome_id] = codon_usage
        consumer_data.mean_gene_length[genome_id] = mean_gene_length

        return consumer_data

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

        return '  Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self, gene_files):
        """Calculate codon usage over a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes in nucleotide space.

        Returns
        -------
        dict of dict : d[genome_id][codon] -> count
           Codon usage of each genome.
        set
           Set with all identified codons.
        dict of dict : d[genome_id][codon] -> length
            Mean length of genes for each stop codon.
        """

        self.logger.info('Calculating codon usage for each genome.')
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None

        parallel = Parallel(self.cpus)
        consumer_data = parallel.run(self._producer, self._consumer, gene_files, progress_func)

        return consumer_data.genome_codon_usage, consumer_data.codon_set, consumer_data.mean_gene_length
