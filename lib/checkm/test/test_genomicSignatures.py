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


import unittest

from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.binStatistics import BinStatistics


class VerifyBinStatistics(unittest.TestCase):
    def testGC(self):
        """Verify computation of GC."""
        binStats = BinStatistics(threads=1)

        seqs = {'S1': 'ACgt', 'S2': 'GGgg', 'S3': 'TTtt', 'S4': 'NNNN'}

        seqStats = defaultdict(dict)
        meanGC, _ = binStats.calculateGC(seqs, seqStats)

        self.assertAlmostEqual(seqStats['S1']['GC'], 0.5)
        self.assertAlmostEqual(seqStats['S2']['GC'], 1.0)
        self.assertAlmostEqual(seqStats['S3']['GC'], 0.0)
        self.assertAlmostEqual(seqStats['S4']['GC'], 0.0)

        self.assertAlmostEqual(meanGC, 6.0 / 12.0)

    def testScaffoldLengthStats(self):
        """Verify computation of scaffold length statistics."""
        binStats = BinStatistics(threads=1)

        scaffolds = {'S1': 'ACGT' + DefaultValues.CONTIG_BREAK + 'ACGT', 'S2': 'ACGTACGT', 'S3': 'TTtt'}

        scaffoldStats = defaultdict(dict)
        maxScaffoldLen, maxContigLen, totalScaffoldBps, _, _, numContigs = binStats.calculateScaffoldLengthStats(scaffolds, scaffoldStats)

        self.assertAlmostEqual(scaffoldStats['S1']['Length'], len(DefaultValues.CONTIG_BREAK) + 8)
        self.assertAlmostEqual(scaffoldStats['S1']['Total contig length'], 8)
        self.assertAlmostEqual(scaffoldStats['S1']['# contigs'], 2)

        self.assertAlmostEqual(maxScaffoldLen, len(DefaultValues.CONTIG_BREAK) + 8)
        self.assertAlmostEqual(maxContigLen, 8)
        self.assertAlmostEqual(totalScaffoldBps, len(DefaultValues.CONTIG_BREAK) + 8 + 8 + 4)
        self.assertAlmostEqual(numContigs, 4)

    def testCodingBases(self):
        """Verify computation of coding bases."""

        binStats = BinStatistics(threads=1)

        aaGenes = {'S1_C1': 'ACGTACGT', 'S1_C2': 'ACGTACGT', 'S3_C1': 'TTtt'}

        seqStats = defaultdict(dict)
        codingBasePairs = binStats._BinStatistics__calculateCodingBases(aaGenes, seqStats)

        self.assertAlmostEqual(seqStats['S1']['# ORFs'], 2)
        self.assertAlmostEqual(seqStats['S1']['Coding bases'], len(aaGenes['S1_C1']) * 3 + len(aaGenes['S1_C2']) * 3)

        self.assertAlmostEqual(codingBasePairs, len(aaGenes['S1_C1']) * 3 + len(aaGenes['S1_C2']) * 3 + len(aaGenes['S3_C1']) * 3)

if __name__ == "__main__":
    unittest.main()
