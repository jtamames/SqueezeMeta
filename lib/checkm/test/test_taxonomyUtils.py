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

from checkm.util.taxonomyUtils import appendTaxonomyRanks, LCA


class VerifyTaxonomyUtils(unittest.TestCase):
    def testAppendTaxonomyRanks(self):
        """Verify computation of base count on mixed-case sequence."""
        r = appendTaxonomyRanks(['k', 'p', 'c'], ranks=3)

        self.assertEqual(r, ['k__k', 'p__p', 'c__c'])

    def testLCA(self):
        """Verify computation of lowest-common ancestor."""
        lca = LCA(['a', 'b', 'c', 'd', 'e', 'f'], ['a', 'b', 'c', 'x', 'y', 'z'])
        self.assertEqual(lca, ['a', 'b', 'c', 'o__unclassified', 'f__unclassified', 'g__unclassified'])

if __name__ == "__main__":
    unittest.main()
