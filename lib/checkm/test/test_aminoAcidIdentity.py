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

from checkm.aminoAcidIdentity import AminoAcidIdentity


class VerifyAminoAcidIdentity(unittest.TestCase):
    def testIdentiticalSeqs(self):
        """Verify computation of AAI on identical sequences."""
        aai = AminoAcidIdentity()

        score = aai.aai('ACGT', 'ACGT')
        self.assertAlmostEqual(score, 1.0)

    def testDiffSeqs(self):
        """Verify computation of AAI on two completely different sequences."""
        aai = AminoAcidIdentity()

        score = aai.aai('ACGT', 'TGCA')
        self.assertAlmostEqual(score, 0.0)

    def testNonOverlappingSeqs(self):
        """Verify computation of AAI on two non-overlapping sequences."""
        aai = AminoAcidIdentity()

        score = aai.aai('ACGT----', '----TGCA')
        self.assertAlmostEqual(score, 0.0)

    def testPartialOverlappingSeqs(self):
        """Verify computation of AAI on two partially overlapping sequences."""
        aai = AminoAcidIdentity()

        score = aai.aai('ACGT--', '--GTAC')
        self.assertAlmostEqual(score, 1.0)

    def testIncompleteSeqs(self):
        """Verify computation of AAI on incomplete sequence."""
        aai = AminoAcidIdentity()

        score = aai.aai('AAAACGTTTT', '---ACGG---')
        self.assertAlmostEqual(score, 3.0 / 4.0)

    def testSimilarSeqs(self):
        """Verify computation of AAI on two similar sequences."""
        aai = AminoAcidIdentity()

        score = aai.aai('ACGT', 'ACGG')
        self.assertAlmostEqual(score, 3.0 / 4.0)

    def testIdenticalInsertionSeqs(self):
        """Verify computation of AAI on two identical sequences with multiple insertions."""
        aai = AminoAcidIdentity()

        score = aai.aai('A-C-G-T', 'A-C-G-T')
        self.assertAlmostEqual(score, 1.0)

    def testDiffInsertionSeqs(self):
        """Verify computation of AAI on two different sequences with multiple insertions."""
        aai = AminoAcidIdentity()

        score = aai.aai('A-C-G-T', 'AACCGGT')
        self.assertAlmostEqual(score, 4.0 / 7.0)


class VerifyStrainHeterogeneity(unittest.TestCase):
    def testNoStrainHetero(self):
        """Verify computation of strain heterogeneity score on divergent sequences."""
        aai = AminoAcidIdentity()

        aaiScores = defaultdict(dict)
        aaiScores['b1'] = {'g1': [0.1], 'g2': [0.1], 'g3': [0.1]}

        aaiHetero, aaiMeanBinHetero = aai.strainHetero(aaiScores, 0.9)

        self.assertAlmostEqual(aaiHetero['b1']['g1'], 0.0)
        self.assertAlmostEqual(aaiHetero['b1']['g2'], 0.0)
        self.assertAlmostEqual(aaiHetero['b1']['g3'], 0.0)

        self.assertAlmostEqual(aaiMeanBinHetero['b1'], 0.0)

    def testAllStrainHetero(self):
        """Verify computation of strain heterogeneity score on highly similar sequences."""
        aai = AminoAcidIdentity()

        aaiScores = defaultdict(dict)
        aaiScores['b1'] = {'g1': [0.95], 'g2': [0.95], 'g3': [0.95]}

        aaiHetero, aaiMeanBinHetero = aai.strainHetero(aaiScores, 0.9)

        self.assertAlmostEqual(aaiHetero['b1']['g1'], 1.0)
        self.assertAlmostEqual(aaiHetero['b1']['g2'], 1.0)
        self.assertAlmostEqual(aaiHetero['b1']['g3'], 1.0)

        self.assertAlmostEqual(aaiMeanBinHetero['b1'], 100.0)

    def testMixedStrainHetero(self):
        """Verify computation of strain heterogeneity score on sequences with variable similarity."""
        aai = AminoAcidIdentity()

        aaiScores = defaultdict(dict)
        aaiScores['b1'] = {'g1': [0.95], 'g2': [0.1], 'g3': [0.1]}

        aaiHetero, aaiMeanBinHetero = aai.strainHetero(aaiScores, 0.9)

        self.assertAlmostEqual(aaiHetero['b1']['g1'], 1.0)
        self.assertAlmostEqual(aaiHetero['b1']['g2'], 0.0)
        self.assertAlmostEqual(aaiHetero['b1']['g3'], 0.0)

        self.assertAlmostEqual(aaiMeanBinHetero['b1'], 1.0 * 100 / 3.0)

    def testMultiCopyStrainHetero(self):
        """Verify computation of strain heterogeneity score when there are multiple copies of a sequence."""
        aai = AminoAcidIdentity()

        aaiScores = defaultdict(dict)
        aaiScores['b1'] = {'g1': [0.95, 0.95, 0.95], 'g2': [0.1, 0.1, 0.1], 'g3': [0.95, 0.1, 0.1]}

        aaiHetero, aaiMeanBinHetero = aai.strainHetero(aaiScores, 0.9)

        self.assertAlmostEqual(aaiHetero['b1']['g1'], 1.0)
        self.assertAlmostEqual(aaiHetero['b1']['g2'], 0.0)
        self.assertAlmostEqual(aaiHetero['b1']['g3'], 1.0 / 3.0)

        self.assertAlmostEqual(aaiMeanBinHetero['b1'], 4.0 * 100 / 9.0)


if __name__ == "__main__":
    unittest.main()
