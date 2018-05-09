###############################################################################
#
# unbinned.py - identify unbinned sequences
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

from checkm.common import checkFileExists

from checkm.util.seqUtils import baseCount, readFasta


class Unbinned():
    def __init__(self):
        self.logger = logging.getLogger()

    def run(self, binFiles, seqFile, outSeqFile, outStatsFile, minSeqLen):
        checkFileExists(seqFile)

        # get list of sequences in bins
        self.logger.info('  Reading binned sequences.')

        binnedSeqs = {}
        totalBinnedBases = 0
        for binFile in binFiles:
            seqs = readFasta(binFile)
            binnedSeqs.update(seqs)
            for seq in seqs.values():
                totalBinnedBases += len(seq)

        self.logger.info('    Read %d (%.2f Mbp) binned sequences.' % (len(binnedSeqs), float(totalBinnedBases) / 1e6))

        # get list of all sequences
        self.logger.info('  Reading all sequences.')
        allSeqs = readFasta(seqFile)
        totalBases = 0
        for seq in allSeqs.values():
            totalBases += len(seq)
        self.logger.info('    Read %d (%.2f Mbp) sequences.' % (len(allSeqs), float(totalBases) / 1e6))

        # write all unbinned sequences
        self.logger.info('  Identifying unbinned sequences >= %d bp.' % minSeqLen)
        seqOut = open(outSeqFile, 'w')

        statsOut = open(outStatsFile, 'w')
        statsOut.write('Sequence Id\tLength\tGC\n')

        unbinnedCount = 0
        unbinnedBases = 0
        for seqId, seq in allSeqs.iteritems():
            if seqId not in binnedSeqs:
                if len(seq) >= minSeqLen:
                    unbinnedCount += 1
                    seqOut.write('>' + seqId + '\n')
                    seqOut.write(seq + '\n')

                    unbinnedBases += len(seq)

                    a, c, g, t = baseCount(seq)

                    statsOut.write('%s\t%d\t%.2f\n' % (seqId, len(seq), float(g + c) * 100 / (a + c + g + t)))

        seqOut.close()
        statsOut.close()

        self.logger.info('    Identified %d (%.2f Mbp) unbinned sequences.' % (unbinnedCount, float(unbinnedBases) / 1e6))

        self.logger.info('')
        self.logger.info('  Percentage of unbinned sequences: %.2f%%' % (unbinnedCount * 100.0 / len(allSeqs)))
        self.logger.info('  Percentage of unbinned bases: %.2f%%' % (unbinnedBases * 100.0 / totalBases))
