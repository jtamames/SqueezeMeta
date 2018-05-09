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

from checkm.genomicSignatures import GenomicSignatures


class VerifyGenomicSignatures(unittest.TestCase):
    def testGenomicSignature(self):
        """Verify computation of genomic signature."""
        gs = GenomicSignatures(K=2, threads=1)

        sig = gs.seqSignature('AACC')
        kmerOrder = gs.canonicalKmerOrder()
        aaIndex = kmerOrder.index('AA')
        acIndex = kmerOrder.index('AC')
        ccIndex = kmerOrder.index('AC')
        atIndex = kmerOrder.index('AT')

        self.assertEqual(sig[aaIndex], 1.0 / 3.0)
        self.assertEqual(sig[acIndex], 1.0 / 3.0)
        self.assertEqual(sig[ccIndex], 1.0 / 3.0)
        self.assertEqual(sig[atIndex], 0)

    def testDistanceZero(self):
        """Verify computation of distances between genomic signatures."""
        gs = GenomicSignatures(K=2, threads=1)

        sig1 = gs.seqSignature('AACC')
        sig2 = gs.seqSignature('AACC')

        dist = gs.distance(sig1, sig2)

        self.assertEqual(dist, 0)

    def testDistanceMax(self):
        """Verify computation of distances between genomic signatures."""
        gs = GenomicSignatures(K=2, threads=1)

        sig1 = gs.seqSignature('AAAA')
        sig2 = gs.seqSignature('GGGG')

        dist = gs.distance(sig1, sig2)

        self.assertAlmostEqual(dist, 2.0)

if __name__ == "__main__":
    unittest.main()
