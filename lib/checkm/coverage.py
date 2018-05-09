###############################################################################
#
# coverage.py - calculate coverage of all sequences
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
import ntpath
import traceback
from collections import defaultdict

import pysam

from checkm.defaultValues import DefaultValues
from checkm.common import reassignStdOut, restoreStdOut, binIdFromFilename
from checkm.util.seqUtils import readFasta

from numpy import mean, sqrt


class CoverageStruct():
    def __init__(self, seqLen, mappedReads, coverage):
        self.seqLen = seqLen
        self.mappedReads = mappedReads
        self.coverage = coverage


class Coverage():
    """Calculate coverage of all sequences."""
    def __init__(self, threads):
        self.logger = logging.getLogger()

        self.totalThreads = threads

    def run(self, binFiles, bamFiles, outFile, bAllReads, minAlignPer, maxEditDistPer, minQC):
        """Calculate coverage of sequences for each BAM file."""

        # determine bin assignment of each sequence
        self.logger.info('  Determining bin assignment of each sequence.')

        seqIdToBinId = {}
        seqIdToSeqLen = {}
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)

            seqs = readFasta(binFile)
            for seqId, seq in seqs.iteritems():
                seqIdToBinId[seqId] = binId
                seqIdToSeqLen[seqId] = len(seq)

        # process each fasta file
        self.logger.info("  Processing %d file(s) with %d threads.\n" % (len(bamFiles), self.totalThreads))

        # make sure all BAM files are sorted
        self.numFiles = len(bamFiles)
        for bamFile in bamFiles:
            if not os.path.exists(bamFile + '.bai'):
                self.logger.error('  [Error] BAM file is either unsorted or not indexed: ' + bamFile + '\n')
                sys.exit()

        # calculate coverage of each BAM file
        coverageInfo = {}
        numFilesStarted = 0
        for bamFile in bamFiles:
            numFilesStarted += 1
            self.logger.info('  Processing %s (%d of %d):' % (ntpath.basename(bamFile), numFilesStarted, len(bamFiles)))

            coverageInfo[bamFile] = mp.Manager().dict()
            coverageInfo[bamFile] = self.__processBam(bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, coverageInfo[bamFile])

        # redirect output
        self.logger.info('  Writing coverage information to file.')
        oldStdOut = reassignStdOut(outFile)

        header = 'Sequence Id\tBin Id\tSequence length (bp)'
        for bamFile in bamFiles:
            header += '\tBam Id\tCoverage\tMapped reads'

        print(header)

        # get length of all seqs
        for bamFile, seqIds in coverageInfo.iteritems():
            for seqId in seqIds.keys():
                seqIdToSeqLen[seqId] = seqIds[seqId].seqLen

        # write coverage stats for all scaffolds to file
        for seqId, seqLen in seqIdToSeqLen.iteritems():
            rowStr = seqId + '\t' + seqIdToBinId.get(seqId, DefaultValues.UNBINNED) + '\t' + str(seqLen)
            for bamFile in bamFiles:
                bamId = binIdFromFilename(bamFile)

                if seqId in coverageInfo[bamFile]:
                    rowStr += '\t%s\t%f\t%d' % (bamId, coverageInfo[bamFile][seqId].coverage, coverageInfo[bamFile][seqId].mappedReads)
                else:
                    rowStr += '\t%s\t%f\t%d' % (bamId, 0, 0)

            print(rowStr)

        # restore stdout
        restoreStdOut(outFile, oldStdOut)

    def __processBam(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, coverageInfo):
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
            workerProc = [mp.Process(target=self.__workerThread, args=(bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__writerThread, args=(coverageInfo, len(refSeqIds), writerQueue))

            writeProc.start()

            for p in workerProc:
                p.start()

            for p in workerProc:
                p.join()

            writerQueue.put((None, None, None, None, None, None, None, None, None, None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            print traceback.format_exc()
            for p in workerProc:
                p.terminate()

            writeProc.terminate()

        return coverageInfo

    def __workerThread(self, bamFile, bAllReads, minAlignPer, maxEditDistPer, minQC, queueIn, queueOut):
        """Process each data item in parallel."""
        while True:
            seqIds, seqLens = queueIn.get(block=True, timeout=None)
            if seqIds == None:
                break

            bamfile = pysam.Samfile(bamFile, 'rb')

            for seqId, seqLen in zip(seqIds, seqLens):
                numReads = 0
                numMappedReads = 0
                numDuplicates = 0
                numSecondary = 0
                numFailedQC = 0
                numFailedAlignLen = 0
                numFailedEditDist = 0
                numFailedProperPair = 0
                coverage = 0

                for read in bamfile.fetch(seqId, 0, seqLen):
                    numReads += 1

                    if read.is_unmapped:
                        pass
                    elif read.is_duplicate:
                        numDuplicates += 1
                    elif read.is_secondary or read.is_supplementary:
                        numSecondary += 1
                    elif read.is_qcfail or read.mapping_quality < minQC:
                        numFailedQC += 1
                    elif read.query_alignment_length < minAlignPer * read.query_length:
                        numFailedAlignLen += 1
                    elif read.get_tag('NM') > maxEditDistPer * read.query_length:
                        numFailedEditDist += 1
                    elif not bAllReads and not read.is_proper_pair:
                        numFailedProperPair += 1
                    else:
                        numMappedReads += 1

                        # Note: the alignment length (query_alignment_length) is used instead of the
                        # read length (query_length) as this bring the calculated coverage
                        # in line with 'samtools depth' (at least when the min
                        # alignment length and edit distance thresholds are zero).
                        coverage += read.query_alignment_length

                coverage = float(coverage) / seqLen

                queueOut.put((seqId, seqLen, coverage, numReads,
                                numDuplicates, numSecondary, numFailedQC,
                                numFailedAlignLen, numFailedEditDist,
                                numFailedProperPair, numMappedReads))

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
            seqId, seqLen, coverage, numReads, numDuplicates, numSecondary, numFailedQC, numFailedAlignLen, numFailedEditDist, numFailedProperPair, numMappedReads = writerQueue.get(block=True, timeout=None)
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

            coverageInfo[seqId] = CoverageStruct(seqLen=seqLen, mappedReads=numMappedReads, coverage=coverage)

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

    def parseCoverage(self, coverageFile):
        """Read coverage information from file."""
        coverageStats = {}
        bHeader = True
        for line in open(coverageFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            binId = lineSplit[1]

            if binId not in coverageStats:
                coverageStats[binId] = {}

            if seqId not in coverageStats[binId]:
                coverageStats[binId][seqId] = {}

            for i in xrange(3, len(lineSplit), 3):
                bamId = lineSplit[i]
                coverage = float(lineSplit[i + 1])
                coverageStats[binId][seqId][bamId] = coverage

        return coverageStats

    def binProfiles(self, coverageFile):
        """Read coverage information for each bin."""
        binCoverages = defaultdict(lambda: defaultdict(list))
        binStats = defaultdict(dict)

        bHeader = True
        for line in open(coverageFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            binId = lineSplit[1]
            seqLen = int(lineSplit[2])

            # calculate mean coverage (weighted by scaffold length)
            # for each bin under each BAM file
            for i in xrange(3, len(lineSplit), 3):
                bamId = lineSplit[i]
                coverage = float(lineSplit[i + 1])
                binCoverages[binId][bamId].append(coverage)

                if bamId not in binStats[binId]:
                    binStats[binId][bamId] = [0, 0]

                binLength = binStats[binId][bamId][0] + seqLen
                weight = float(seqLen) / binLength
                meanBinCoverage = coverage * weight + binStats[binId][bamId][1] * (1 - weight)

                binStats[binId][bamId] = [binLength, meanBinCoverage]

        profiles = defaultdict(dict)
        for binId in binStats:
            for bamId, stats in binStats[binId].iteritems():
                binLength, meanBinCoverage = stats
                coverages = binCoverages[binId][bamId]

                varCoverage = 0
                if len(coverages) > 1:
                    varCoverage = mean(map(lambda x: (x - meanBinCoverage) ** 2, coverages))

                profiles[binId][bamId] = [meanBinCoverage, sqrt(varCoverage)]

        return profiles
