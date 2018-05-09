###############################################################################
#
# aminoAcidIdentity.py - calculate AAI between aligned marker genes
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

import os
import sys
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.common import getBinIdsFromOutDir
from checkm.util.seqUtils import readFasta


class AminoAcidIdentity():
    """Calculate AAI between sequences aligned to an HMM."""
    def __init__(self):
        self.logger = logging.getLogger()
        self.aaiRawScores = defaultdict(dict)
        self.aaiHetero = defaultdict(dict)
        self.aaiMeanBinHetero = {}

    def run(self, aaiStrainThreshold, outDir, alignmentOutputFile):
        """Calculate AAI between input alignments."""

        self.logger.info('  Calculating AAI between multi-copy marker genes.')

        if alignmentOutputFile:
            fout = open(alignmentOutputFile, 'w')

        # calculate AAI for duplicate marker genes
        binIds = getBinIdsFromOutDir(outDir)
        aaiOutputDir = os.path.join(outDir, 'storage', 'aai_qa')
        for binId in binIds:
            binPath = os.path.join(aaiOutputDir, binId)
            if not os.path.exists(binPath):
                continue

            for f in os.listdir(binPath):
                if not f.endswith('.masked.faa'):
                    continue

                markerId = f[0:f.find('.')]

                seqs = readFasta(os.path.join(binPath, f))

                # calculate AAI between all pairs of seqs
                for i in xrange(0, len(seqs)):
                    seqIdI = seqs.keys()[i]
                    binIdI = seqIdI[0:seqIdI.find(DefaultValues.SEQ_CONCAT_CHAR)]

                    seqI = seqs[seqIdI]

                    for j in xrange(i + 1, len(seqs)):
                        seqIdJ = seqs.keys()[j]
                        binIdJ = seqIdJ[0:seqIdJ.find(DefaultValues.SEQ_CONCAT_CHAR)]

                        seqJ = seqs[seqIdJ]

                        if binIdI == binIdJ:
                            aai = self.aai(seqI, seqJ)

                            if alignmentOutputFile:
                                fout.write(binId + ',' + markerId + '\n')
                                fout.write(seqIdI + '\t' + seqI + '\n')
                                fout.write(seqIdJ + '\t' + seqJ + '\n')
                                fout.write('AAI: %.3f\n' % aai)
                                fout.write('\n')

                            if binIdI not in self.aaiRawScores:
                                self.aaiRawScores[binIdI] = defaultdict(list)
                            self.aaiRawScores[binIdI][markerId].append(aai)
                        else:
                            # something is wrong as the bin Ids should always be the same
                            self.logger.error('  [Error] Bin ids do not match.')
                            sys.exit()

        if alignmentOutputFile:
            fout.close()

        # calculate strain heterogeneity for each marker gene in each bin
        self.aaiHetero, self.aaiMeanBinHetero = self.strainHetero(self.aaiRawScores, aaiStrainThreshold)

    def strainHetero(self, aaiScores, aaiStrainThreshold):
        """Calculate strain heterogeneity."""
        aaiHetero = defaultdict(dict)
        aaiMeanBinHetero = {}

        for binId, markerIds in aaiScores.iteritems():
            strainCount = 0
            multiCopyPairs = 0

            aaiHetero[binId] = {}

            for markerId, aaiScores in markerIds.iteritems():
                localStrainCount = 0
                for aaiScore in aaiScores:
                    multiCopyPairs += 1
                    if aaiScore > aaiStrainThreshold:
                        strainCount += 1
                        localStrainCount += 1

                strainHetero = float(localStrainCount) / len(aaiScores)
                aaiHetero[binId][markerId] = strainHetero

            aaiMeanBinHetero[binId] = 100 * float(strainCount) / multiCopyPairs

        return aaiHetero, aaiMeanBinHetero

    def aai(self, seq1, seq2):
        """Calculate amino acid identity between sequences."""
        assert len(seq1) == len(seq2)

        # calculation of AAI should ignore missing data at
        # the start of end of each sequence
        startIndex = 0
        for i in xrange(0, len(seq1)):
            if seq1[i] == '-' or seq2[i] == '-':
                startIndex = i + 1
            else:
                break

        endIndex = len(seq1)
        for i in xrange(len(seq1) - 1, 0, -1):
            if seq1[i] == '-' or seq2[i] == '-':
                endIndex = i
            else:
                break

        mismatches = 0
        seqLen = 0
        for i in xrange(startIndex, endIndex):
            if seq1[i] != seq2[i]:
                mismatches += 1
                seqLen += 1
            elif seq1[i] == '-' and seq2[i] == '-':
                pass
            else:
                seqLen += 1

        if seqLen == 0:
            return 0.0

        return 1.0 - (float(mismatches) / seqLen)
