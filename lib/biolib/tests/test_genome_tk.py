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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
from unittest import TestCase

from nose.tools import assert_equals

import biolib.genome_tk as genome_tk


class GenomeTkTests(TestCase):
    @classmethod
    def setup_class(self):
        """Setup class variables before any tests."""
        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(self.test_dir, 'test_data')

    def test_unique(self):
        """Verify GenomeTk.unique()"""
        unique_test_data_dir = os.path.join(self.test_data_dir, 'unique')

        genome_files = [os.path.join(unique_test_data_dir, f) for f in os.listdir(unique_test_data_dir)]

        duplicates = genome_tk.unique(genome_files)

        gt = {'genome2': {'genome1': set(['c_dup'])}, 'genome1': {'genome2': set(['c_dup']), 'genome1': ['b_dup']}}
        assert_equals(duplicates, gt)
