###############################################################################
#
# binComparer.py - compare two sets of bins
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

import logging
from collections import defaultdict

from common import binIdFromFilename
from checkm.util.seqUtils import readFasta, readFastaSeqIds


class BinComparer(object):
    def __init__(self):
        self.logger = logging.getLogger()

    def __readBins(self, binFiles):
        bins = {}
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)
            bins[binId] = set(readFastaSeqIds(binFile))

        return bins

    def __binningStats(self, bins, seqLens):
        totalUniqueBinnedSeqs = 0
        totalUniqueBinnedBases = 0

        binStats = {}
        processedSeqs = set()
        repeats = set()
        for binId, seqs in bins.iteritems():
            numBinnedBases = 0
            for seqId in seqs:
                numBinnedBases += seqLens[seqId]
                if seqId not in processedSeqs:
                    processedSeqs.add(seqId)
                    totalUniqueBinnedBases += seqLens[seqId]
                    totalUniqueBinnedSeqs += 1
                else:
                    repeats.add(seqId)

            binStats[binId] = [len(seqs), numBinnedBases]

        return binStats, totalUniqueBinnedSeqs, totalUniqueBinnedBases, len(repeats)

    def report(self, binFiles1, binFiles2, seqFile, outputFile):
        # determine total number of sequences
        self.logger.info('  Reading sequences.')
        seqs = readFasta(seqFile)

        seqLens = {}
        totalBases = 0
        numSeq1K = 0
        totalBases1K = 0
        numSeq5K = 0
        totalBases5K = 0
        for seqId, seq in seqs.iteritems():
            seqLen = len(seq)
            seqLens[seqId] = seqLen
            totalBases += seqLen
            if seqLen >= 1000:
                numSeq1K += 1
                totalBases1K += seqLen
            if seqLen >= 5000:
                numSeq5K += 1
                totalBases5K += seqLen

        # determine sequences in each bin
        bins1 = self.__readBins(binFiles1)
        bins2 = self.__readBins(binFiles2)

        # determine bin stats
        binStats1, totalUniqueBinnedSeqs1, totalUniqueBinnedBases1, numRepeats1 = self.__binningStats(bins1, seqLens)
        binStats2, totalUniqueBinnedSeqs2, totalUniqueBinnedBases2, numRepeats2 = self.__binningStats(bins2, seqLens)

        # sort bins by size
        binStats1 = sorted(binStats1.iteritems(), key=lambda x: x[1][1], reverse=True)
        binStats2 = sorted(binStats2.iteritems(), key=lambda x: x[1][1], reverse=True)

        # report summary results
        self.logger.info('    Total seqs = %d (%.2f Mbp)' % (len(seqs), float(totalBases) / 1e6))
        self.logger.info('      # seqs > 1 kbp = %d (%.2f Mbp)' % (numSeq1K, float(totalBases1K) / 1e6))
        self.logger.info('      # seqs > 5 kbp = %d (%.2f Mbp)' % (numSeq5K, float(totalBases5K) / 1e6))
        self.logger.info('')
        self.logger.info('  Binned seqs statistics:')
        self.logger.info('    1) # bins: %s, # binned seqs: %d (%.2f%%), # binned bases: %.2f Mbp (%.2f%%), # seqs in multiple bins: %d'
                                % (len(bins1),
                                   totalUniqueBinnedSeqs1,
                                   float(totalUniqueBinnedSeqs1) * 100 / len(seqs),
                                   float(totalUniqueBinnedBases1) / 1e6,
                                   float(totalUniqueBinnedBases1) * 100 / totalBases,
                                   numRepeats1))
        self.logger.info('    2) # bins: %s, # binned seqs: %d (%.2f%%), # binned bases: %.2f Mbp (%.2f%%), # seqs in multiple bins: %d'
                                % (len(bins2),
                                   totalUniqueBinnedSeqs2,
                                   float(totalUniqueBinnedSeqs2) * 100 / len(seqs),
                                   float(totalUniqueBinnedBases2) / 1e6,
                                   float(totalUniqueBinnedBases2) * 100 / totalBases,
                                   numRepeats2))

        # output report
        fout = open(outputFile, 'w')
        for data in binStats2:
            fout.write('\t' + data[0])
        fout.write('\tunbinned\t# seqs\t# bases (Mbp)\tBest match\t% bases in common\t% seqs in common\n')

        totalSeqsInCommon2 = defaultdict(int)
        maxBasesInCommon2 = defaultdict(int)
        maxSeqsInCommon2 = defaultdict(int)
        bestMatchingBin2 = {}
        binnedSeqs2 = defaultdict(set)
        for data1 in binStats1:
            binId1 = data1[0]
            fout.write(binId1)

            seqs1 = bins1[binId1]

            maxBasesInCommon = 0
            maxSeqsInCommon = 0
            bestMatchingBin = 'n/a'
            binnedSeqs = set()
            for data2 in binStats2:
                binId2 = data2[0]
                seqs2 = bins2[binId2]

                seqsInCommon = seqs1.intersection(seqs2)
                binnedSeqs.update(seqsInCommon)
                numSeqsInCommon = len(seqsInCommon)
                fout.write('\t' + str(numSeqsInCommon))

                basesInCommon = 0
                for seqId in seqsInCommon:
                    basesInCommon += seqLens[seqId]

                if basesInCommon > maxBasesInCommon:
                    maxBasesInCommon = basesInCommon
                    maxSeqsInCommon = numSeqsInCommon
                    bestMatchingBin = binId2

                if basesInCommon > maxBasesInCommon2[binId2]:
                    maxBasesInCommon2[binId2] = basesInCommon
                    maxSeqsInCommon2[binId2] = numSeqsInCommon
                    bestMatchingBin2[binId2] = binId1

                binnedSeqs2[binId2].update(seqsInCommon)
            fout.write('\t%d\t%d\t%.2f\t%s\t%.2f\t%.2f\n' % (len(seqs1) - len(binnedSeqs),
                                                             data1[1][0],
                                                             float(data1[1][1]) / 1e6,
                                                             bestMatchingBin,
                                                             float(maxBasesInCommon) * 100 / data1[1][1],
                                                             float(maxSeqsInCommon) * 100 / data1[1][0],
                                                             ))

        fout.write('unbinned')
        for data in binStats2:
            binId = data[0]
            fout.write('\t%d' % (len(bins2[binId]) - len(binnedSeqs2[binId])))
        fout.write('\n')

        fout.write('# seqs')
        for data in binStats2:
            fout.write('\t%d' % data[1][0])
        fout.write('\n')

        fout.write('# bases (Mbp)')
        for data in binStats2:
            fout.write('\t%.2f' % (float(data[1][1]) / 1e6))
        fout.write('\n')

        fout.write('Best match')
        for data in binStats2:
            binId = data[0]
            fout.write('\t%s' % bestMatchingBin2.get(binId, 'n/a'))
        fout.write('\n')

        fout.write('% bases in common')
        for data in binStats2:
            binId = data[0]
            fout.write('\t%.2f' % (float(maxBasesInCommon2[binId]) * 100 / data[1][1]))
        fout.write('\n')

        fout.write('% seqs in common')
        for data in binStats2:
            binId = data[0]
            fout.write('\t%.2f' % (float(maxSeqsInCommon2[binId]) * 100 / data[1][0]))
        fout.write('\n')

        fout.close()
