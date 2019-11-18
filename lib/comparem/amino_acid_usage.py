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


class AminoAcidUsage(object):
    """Calculate amino acid usage over a set of genomes."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def amino_acid_usage(self, seqs):
        """ Calculate amino acid usage within sequences.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.

        Returns
        -------
        dict : dict[aa] -> count
            Occurrence of each amino acid.
        """
        aa_usage = defaultdict(int)

        for _seqId, seq in seqs.items():
            for aa in seq:
                if aa != '*':
                    aa = aa.upper()
                    aa_usage[aa] += 1

        return aa_usage

    def _producer(self, gene_file):
        """Calculates amino acid usage of a genome.

        Parameters
        ----------
        gene_file : str
            Fasta file containing amino acid sequences.

        Returns
        -------
        str
           Unique identifier of genome.
        dict : dict[aa] -> count
            Occurrence of each amino acid.
        """

        genome_id = ntpath.basename(gene_file)
        genome_id = genome_id.replace('.genes.faa', '')
        genome_id = os.path.splitext(genome_id)[0]

        seqs = seq_io.read_fasta(gene_file)
        aa_usage = self.amino_acid_usage(seqs)

        return [genome_id, aa_usage]

    def _consumer(self, produced_data, consumer_data):
        """Consume results from producer processes.

         Parameters
        ----------
        produced_data : list -> [genome_id, aa_usage]
            Unique id of a genome followed by a dictionary
            indicating its amino acid usage.
        consumer_data : namedtuple
            Set of amino acids observed across all genomes (aa_set),
            along with the amino acid usage of each genome (genome_aa_usage).

        Returns
        -------
        consumer_data
            The consumer data structure or None must be returned
        """

        if consumer_data == None:
            # setup data to be returned by consumer
            ConsumerData = namedtuple('ConsumerData', 'aa_set genome_aa_usage')
            consumer_data = ConsumerData(set(), dict())

        genome_id, aa_usage = produced_data

        consumer_data.aa_set.update(list(aa_usage.keys()))
        consumer_data.genome_aa_usage[genome_id] = aa_usage

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
        """Calculate amino acid usage over a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes.

        Returns
        -------
        dict of dict : dict[genome_id][aa] -> count
           Amino acid usage of each genome.
        set
           Set with all identified amino acids.
        """

        self.logger.info('Calculating amino acid usage for each genome:')
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None

        parallel = Parallel(self.cpus)
        consumer_data = parallel.run(self._producer, self._consumer, gene_files, progress_func)

        return consumer_data.genome_aa_usage, consumer_data.aa_set
