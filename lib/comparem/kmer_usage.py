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
from collections import namedtuple

import biolib.seq_io as seq_io
from biolib.genomic_signature import GenomicSignature
from biolib.parallel import Parallel


class KmerUsage(object):
    """Calculate kmer usage over a set of genomes.

    The implementation for calculating genomic signatures
    is not optimized for speed. As such, this class is
    useful for k <= 8.
    """

    def __init__(self, k, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.k = k
        self.cpus = cpus

        self.logger.info('Calculating unique kmers of size k = %d.' % self.k)
        self.signatures = GenomicSignature(self.k)

    def _producer(self, genome_file):
        """Calculates kmer usage of a genome.

        Parameters
        ----------
        genome_file : str
            Fasta file containing genomic sequences.

        Returns
        -------
        str
           Unique identifier of genome.
        dict : d[kmer] -> count
            Occurrence of each kmer.
        """

        genome_id = ntpath.basename(genome_file)
        genome_id = os.path.splitext(genome_id)[0]

        seqs = seq_io.read_fasta(genome_file)
        kmer_usage = self.signatures.calculate(seqs)

        return (genome_id, kmer_usage)

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : list -> [genome_id, kmer_usage]
            Unique id of a genome followed by a dictionary
            indicating its kmer usage.
        consumer_data : namedtuple
            Set of kmers observed across all genomes (kmer_set),
            along with the kmer usage of each genome (genome_kmer_usage).

        Returns
        -------
        consumer_data
            The consumer data structure or None must be returned
        """

        if consumer_data == None:
            # setup data to be returned by consumer
            ConsumerData = namedtuple('ConsumerData', 'kmer_set genome_kmer_usage')
            consumer_data = ConsumerData(set(), dict())

        genome_id, kmer_usage = produced_data

        consumer_data.kmer_set.update(list(kmer_usage.keys()))
        consumer_data.genome_kmer_usage[genome_id] = kmer_usage

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

    def run(self, genome_files):
        """Calculate kmer usage over a set of genomes.

        Parameters
        ----------
        genome_files : list
            Fasta files containing genomic sequences in nucleotide space.

        Returns
        -------
        dict of dict : d[genome_id][kmer] -> count
           Kmer usage of each genome.
        set
           Set with all identified kmers.
        """

        self.logger.info('Calculating kmer usage for each genome.')
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None

        parallel = Parallel(self.cpus)
        consumer_data = parallel.run(self._producer, self._consumer, genome_files, progress_func)

        return consumer_data.genome_kmer_usage, consumer_data.kmer_set
