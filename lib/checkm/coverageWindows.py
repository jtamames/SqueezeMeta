###############################################################################
#
# coverageWindows.py - calculate coverage of windows within sequences
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
import os
import multiprocessing as mp
import logging

import pysam

from numpy import zeros


class ReadLoader:
    """Callback for counting aligned reads with pysam.fetch"""

    def __init__(self, refLength, bAllReads, minAlignPer, maxEditDistPer):
        self.bAllReads = bAllReads
        self.minAlignPer = minAlignPer
        self.maxEditDistPer = maxEditDistPer

        self.numReads = 0
        self.numMappedReads = 0
        self.numDuplicates = 0
        self.numSecondary = 0
        self.numFailedQC = 0
        self.numFailedAlignLen = 0
        self.numFailedEditDist = 0
        self.numFailedProperPair = 0

        self.coverage = zeros(refLength)

    def __call__(self, read):
        self.numReads += 1

        if read.is_unmapped:
            pass
        elif read.is_duplicate:
            self.numDuplicates += 1
        elif read.is_secondary:
            self.numSecondary += 1
        elif read.is_qcfail:
            self.numFailedQC += 1
        elif read.alen < self.minAlignPer * read.rlen:
            self.numFailedAlignLen += 1
        elif read.opt('NM') > self.maxEditDistPer * read.rlen:
            self.numFailedEditDist += 1
        elif not self.bAllReads and not read.is_proper_pair:
            self.numFailedProperPair += 1
        else:
            self.numMappedReads += 1

            # Note: the alignment length (alen) is used instead of the
            # read length (rlen) as this bring the calculated coverage
            # in line with 'samtools depth' (at least when the min
            # alignment length and edit distance thresholds are zero).
            self.coverage[read.pos:read.pos + read.alen] += 1.0


class CoverageStruct():
    def __init__(self, seqLen, mappedReads, coverage):
        self.seqLen = seqLen
        self.mappedReads = mappedReads
        self.coverage = coverage


class CoverageWindows():
    """Calculate coverage of all sequences."""
    def __init__(self, threads):
        self.logger = logging.getLogger()

        self.totalThreads = threads

    def run(self, binFiles, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize):
        """Calculate coverage of full sequences and windows."""

        # make sure BAM file is sorted
        if not os.path.exists(bamFile + '.bai'):
            self.logger.error('  [Error] BAM file is not sorted: ' + bamFile + '\n')
            sys.exit()

        # calculate coverage of each BAM file
        self.logger.info('  Calculating coverage of windows.')
        coverageInfo = mp.Manager().dict()
        coverageInfo = self.__processBam(bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, coverageInfo)

        return coverageInfo

    def __processBam(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, coverageInfo):
        """Calculate coverage of sequences in BAM file."""

        # determine coverage for each reference sequence
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        bamfile = pysam.Samfile(bamFile, 'rb')
        refSeqIds = bamfile.references
        refSeqLens = bamfile.lengths

        # populate each thread with reference sequence to process
        # Note: reference sequences are sorted by number of mapped reads
        # so it is important to distribute reads in a sensible way to each
        # of the threads
        refSeqLists = [[] for _ in range(self.totalThreads)]
        refLenLists = [[] for _ in range(self.totalThreads)]

        threadIndex = 0
        incDir = 1
        for refSeqId, refLen in zip(refSeqIds, refSeqLens):
            refSeqLists[threadIndex].append(refSeqId)
            refLenLists[threadIndex].append(refLen)

            threadIndex += incDir
            if threadIndex == self.totalThreads:
                threadIndex = self.totalThreads - 1
                incDir = -1
            elif threadIndex == -1:
                threadIndex = 0
                incDir = 1

        for i in range(self.totalThreads):
            workerQueue.put((refSeqLists[i], refLenLists[i]))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        try:
            workerProc = [mp.Process(target=self.__workerThread, args=(bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__writerThread, args=(coverageInfo, len(refSeqIds), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None, None, None, None, None, None, None, None, None, None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

        return coverageInfo

    def __workerThread(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, windowSize, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            seqIds, seqLens = queueIn.get(block=True, timeout=None)
            if seqIds == None:
                break

            bamfile = pysam.Samfile(bamFile, 'rb')

            for seqId, seqLen in zip(seqIds, seqLens):
                readLoader = ReadLoader(seqLen, bAllReads, minAlignPer, maxEditDistPer)
                bamfile.fetch(seqId, 0, seqLen, callback=readLoader)

                start = 0
                end = windowSize
                windowCoverages = []
                while(end < seqLen):
                    windowCoverages.append(sum(readLoader.coverage[start:end]) / windowSize)

                    start = end
                    try:
                        end += windowSize
                    except:
                        print '*****************'
                        print end
                        print windowSize
                        print '******************'

                coverage = float(sum(readLoader.coverage)) / seqLen

                queueOut.put((seqId, seqLen, coverage, windowCoverages, readLoader.numReads,
                                readLoader.numDuplicates, readLoader.numSecondary, readLoader.numFailedQC,
                                readLoader.numFailedAlignLen, readLoader.numFailedEditDist,
                                readLoader.numFailedProperPair, readLoader.numMappedReads))

            bamfile.close()

    def __writerThread(self, coverageInfo, numRefSeqs, writerQueue):
        """Store or write results of worker threads in a single thread."""
        totalReads = 0
        totalDuplicates = 0
        totalSecondary = 0
        totalFailedQC = 0
        totalFailedAlignLen = 0
        totalFailedEditDist = 0
        totalFailedProperPair = 0
        totalMappedReads = 0

        processedRefSeqs = 0
        while True:
            seqId, seqLen, coverage, windowCoverages, numReads, numDuplicates, numSecondary, numFailedQC, numFailedAlignLen, numFailedEditDist, numFailedProperPair, numMappedReads = writerQueue.get(block=True, timeout=None)
            if seqId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                processedRefSeqs += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) reference sequences.' % (processedRefSeqs, numRefSeqs, float(processedRefSeqs) * 100 / numRefSeqs)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

                totalReads += numReads
                totalDuplicates += numDuplicates
                totalSecondary += numSecondary
                totalFailedQC += numFailedQC
                totalFailedAlignLen += numFailedAlignLen
                totalFailedEditDist += numFailedEditDist
                totalFailedProperPair += numFailedProperPair
                totalMappedReads += numMappedReads

            coverageInfo[seqId] = [coverage, windowCoverages]

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

            print ''
            print '    # total reads: %d' % totalReads
            print '      # properly mapped reads: %d (%.1f%%)' % (totalMappedReads, float(totalMappedReads) * 100 / totalReads)
            print '      # duplicate reads: %d (%.1f%%)' % (totalDuplicates, float(totalDuplicates) * 100 / totalReads)
            print '      # secondary reads: %d (%.1f%%)' % (totalSecondary, float(totalSecondary) * 100 / totalReads)
            print '      # reads failing QC: %d (%.1f%%)' % (totalFailedQC, float(totalFailedQC) * 100 / totalReads)
            print '      # reads failing alignment length: %d (%.1f%%)' % (totalFailedAlignLen, float(totalFailedAlignLen) * 100 / totalReads)
            print '      # reads failing edit distance: %d (%.1f%%)' % (totalFailedEditDist, float(totalFailedEditDist) * 100 / totalReads)
            print '      # reads not properly paired: %d (%.1f%%)' % (totalFailedProperPair, float(totalFailedProperPair) * 100 / totalReads)
            print ''
