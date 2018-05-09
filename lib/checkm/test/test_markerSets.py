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

from checkm.markerSets import MarkerSet, BinMarkerSets


class VerifyMarkerSet(unittest.TestCase):
    def testMarkerSet(self):
        """Verify marker set data structure."""

        markers = [set(['a', 'b']), set(['c'])]
        ms = MarkerSet(0, 'k__Bacteria', 100, markers)

        markerGenes, markerSets = ms.size()
        self.assertEqual(markerGenes, 3)
        self.assertEqual(markerSets, 2)

        self.assertEqual(ms.numMarkers(), 3)
        self.assertEqual(ms.numSets(), 2)

        self.assertEqual(ms.getMarkerGenes(), set(['a', 'b', 'c']))


class VerifyBinMarkerSets(unittest.TestCase):
    def testBinMarkerSets(self):
        """Verify bin marker set data structure."""

        bms = BinMarkerSets(0, BinMarkerSets.TAXONOMIC_MARKER_SET)

        ms1 = MarkerSet(1, 'k__Bacteria', 100, [set(['a', 'b']), set(['c'])])
        bms.addMarkerSet(ms1)

        ms2 = MarkerSet(2, 'k__Bacteria', 100, [set(['d', 'e']), set(['f'])])
        bms.addMarkerSet(ms2)

        self.assertEqual(bms.getMarkerGenes(), set(['a', 'b', 'c', 'd', 'e', 'f']))
        self.assertEqual(bms.mostSpecificMarkerSet(), ms1)
        self.assertEqual(bms.selectedMarkerSet(), ms1)

if __name__ == "__main__":
    unittest.main()
