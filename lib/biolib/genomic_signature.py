#!/usr/bin/env python

###############################################################################
#
# calculateBoundsDeltaGC.py - find confidence intervals for GC distribution
#
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

#from string import maketrans # Commented by FPS NOV 15 2019, we now use str.maketrans instead
import logging

import numpy as np


class GenomicSignature(object):
    """Genomic signature for arbitrary k.

    This implementation for calculating genomic signatures
    is not optimized for speed or memory. As such, this class
    is best used for k <= 10.
    """

    def __init__(self, k):
        """Initialize."""

        self.logger = logging.getLogger()

        self.k = k
        self.compl = str.maketrans('ACGT', 'TGCA') # FPS NOV 15 2019
        self.kmer_cols, self.kmer_index = self._identify_kmers()

    def _identify_kmers(self):
        """Determine unique kmers of a given length."""

        # determine all mers of a given length
        base_words = ("A", "C", "G", "T")
        mers = ["A", "C", "G", "T"]
        for _ in range(1, self.k):
            working_list = []
            for mer in mers:
                for char in base_words:
                    working_list.append(mer + char)
            mers = working_list

        # pare down kmers based on lexicographical ordering
        kmer_set = set()
        for mer in mers:
            kmer = self._lexicographically_lowest(mer)
            kmer_set.add(kmer)

        canonical_kmer_list = list(kmer_set)
        canonical_kmer_list.sort()

        # create mapping from kmers to their canonical order position
        kmer_index = {}
        for index, kmer in enumerate(canonical_kmer_list):
            kmer_index[kmer] = index
            kmer_index[self._rev_comp(kmer)] = index

        return canonical_kmer_list, kmer_index

    def _lexicographically_lowest(self, seq):
        """Return the lexicographically lowest form of a sequence."""
        rseq = self._rev_comp(seq)
        if(seq < rseq):
            return seq
        return rseq

    def _rev_comp(self, seq):
        """Return the reverse complement of a sequence."""
        return seq.translate(self.compl)[::-1]

    def calculate(self, seqs):
        """Calculate genomic signature of sequences.

        Parameters
        ----------
        seqs : d[seq_id] -> seq
            Sequences indexed by sequence id.

        Returns
        -------
        list
            Count of each kmer in the set of sequences.
        """

        sig = [0] * len(self.kmer_cols)
        for seq in list(seqs.values()):
            tmp_seq = seq.upper()
            for i in range(0, len(tmp_seq) - self.k + 1):
                try:
                    kmer_index = self.kmer_index[tmp_seq[i:i + self.k]]
                    sig[kmer_index] += 1
                except KeyError:
                    # unknown kmer due to an ambiguous character
                    pass
                    
        total_kmers = sum(sig)
        for i, c in enumerate(sig):
            sig[i] = float(c)/total_kmers
                    
        return sig

    def canonical_order(self):
        """Canonical order of kmers."""
        return self.kmer_cols

    def seq_signature(self, seq):
        """Calculate genomic signature of a sequence.

        Parameters
        ----------
        seq : str
            Sequences.

        Returns
        -------
        list
            Count of each kmer in the canonical order.
        """

        tmp_seq = seq.upper()

        sig = [0] * len(self.kmer_cols)
        for i in range(0, len(tmp_seq) - self.k + 1):
            try:
                kmer_index = self.kmer_index[tmp_seq[i:i + self.k]]
                sig[kmer_index] += 1
            except KeyError:
                # unknown kmer due to an ambiguous character
                pass

        return sig

    def manhattan(self, sig1, sig2):
        """Calculate Manhattan distance between genomic signatures.

        Parameters
        ----------
        sig1 : list of kmer counts in canonical order
            First genomic signature.
        sig2 : list of kmer counts in canonical order
            Second genomic signature.

        Returns
        -------
        float
            Manhattan distance between signatures.
        """
        return np.sum(np.abs(np.array(sig1) - np.array(sig2)))

    def read(self, kmer_profile_file):
        """Read genomic signatures.

        Parameters
        ----------
        kmer_profile_file : str
            Name of file to read.

        Returns
        -------
        dict : d[seq_id] -> counts in canonical order
            Count of each kmer.
        """

        sig = {}
        with open(kmer_profile_file) as f:
            next(f)
            for line in f:
                line_split = line.split('\t')
                sig[line_split[0]] = np.array([float(x) for x in line_split[1:]])

        return sig
