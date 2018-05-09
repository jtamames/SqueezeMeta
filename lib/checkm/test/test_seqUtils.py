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

from checkm.util.seqUtils import baseCount, calculateN50


class VerifySeqUtils(unittest.TestCase):
    def testBaseCount(self):
        """Verify computation of base count on mixed-case sequence."""
        a, c, g, t = baseCount('ACGTacgtNnUu')
        self.assertEqual(a, 2)
        self.assertEqual(c, 2)
        self.assertEqual(g, 2)
        self.assertEqual(t, 4)

    def testScaffoldLengthStats(self):
        """Verify computation of N50."""
        n50 = calculateN50([1, 1, 2, 2, 2, 2, 10])
        self.assertEqual(n50, 10)

if __name__ == "__main__":
    unittest.main()
