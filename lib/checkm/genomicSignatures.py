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

import sys
import multiprocessing as mp
from string import maketrans
import logging

import numpy as np

from checkm.util.seqUtils import readFasta


class GenomicSignatures(object):
    def __init__(self, K, threads):
        self.logger = logging.getLogger()

        self.K = K
        self.compl = maketrans('ACGT', 'TGCA')
        self.kmerCols, self.kmerToCanonicalIndex = self.__makeKmerColNames()

        self.totalThreads = threads

    def __makeKmerColNames(self):
        """Work out unique kmers."""

        # determine all mers of a given length
        baseWords = ("A", "C", "G", "T")
        mers = ["A", "C", "G", "T"]
        for _ in range(1, self.K):
            workingList = []
            for mer in mers:
                for char in baseWords:
                    workingList.append(mer + char)
            mers = workingList

        # pare down kmers based on lexicographical ordering
        retList = []
        for mer in mers:
            kmer = self.__lexicographicallyLowest(mer)
            if kmer not in retList:
                retList.append(kmer)

        sorted(retList)

        # create mapping from kmers to their canonical order position
        kmerToCanonicalIndex = {}
        for index, kmer in enumerate(retList):
            kmerToCanonicalIndex[kmer] = index
            kmerToCanonicalIndex[self.__revComp(kmer)] = index

        return retList, kmerToCanonicalIndex

    def __lexicographicallyLowest(self, seq):
        """Return the lexicographically lowest form of this sequence."""
        rseq = self.__revComp(seq)
        if(seq < rseq):
            return seq
        return rseq

    def __revComp(self, seq):
        """Return the reverse complement of a sequence."""
        # build a dictionary to know what letter to switch to
        return seq.translate(self.compl)[::-1]

    def __calculateResults(self, queueIn, queueOut):
        """Calculate genomic signature of sequences in parallel."""
        while True:
            seqId, seq = queueIn.get(block=True, timeout=None)
            if seqId == None:
                break

            sig = self.seqSignature(seq)

            queueOut.put((seqId, sig))

    def __storeResults(self, seqFile, outputFile, totalSeqs, writerQueue):
        """Store genomic signatures to file."""

        # write header
        fout = open(outputFile, 'w')
        fout.write('Sequence Id')
        for kmer in self.canonicalKmerOrder():
            fout.write('\t' + kmer)
        fout.write('\n')

        numProcessedSeq = 0
        while True:
            seqId, sig = writerQueue.get(block=True, timeout=None)
            if seqId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedSeq += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) sequences.' % (numProcessedSeq, totalSeqs, float(numProcessedSeq) * 100 / totalSeqs)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

            fout.write(seqId)
            fout.write('\t' + '\t'.join(map(str, sig)))
            fout.write('\n')

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

        fout.close()

    def canonicalKmerOrder(self):
        return self.kmerCols

    def seqSignature(self, seq):
        sig = [0] * len(self.kmerCols)

        tmp_seq = seq.upper()

        numMers = len(tmp_seq) - self.K + 1
        for i in range(0, numMers):
            try:
                kmerIndex = self.kmerToCanonicalIndex[tmp_seq[i:i + self.K]]
                sig[kmerIndex] += 1  # Note: a numpy array would be slow here due to this single element increment
            except KeyError:
                # unknown kmer (e.g., contains a N)
                pass

        # normalize
        sig = np.array(sig, dtype=float)
        sig /= np.sum(sig)

        return sig

    def calculate(self, seqFile, outputFile):
        """Calculate genomic signature of each sequence."""

        self.logger.info('  Determining tetranucleotide signature of each sequence.')

        # process each sequence in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        seqs = readFasta(seqFile)

        for seqId, seq in seqs.iteritems():
            workerQueue.put((seqId, seq))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        try:
            calcProc = [mp.Process(target=self.__calculateResults, args=(workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__storeResults, args=(seqFile, outputFile, len(seqs), writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

    def distance(self, sig1, sig2):
        return np.sum(np.abs(sig1 - sig2))

    def read(self, tetraProfileFile):
        sig = {}
        with open(tetraProfileFile) as f:
            next(f)
            for line in f:
                lineSplit = line.split('\t')
                sig[lineSplit[0]] = np.array([float(x) for x in lineSplit[1:]])

        return sig
