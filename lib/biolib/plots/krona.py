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
import shutil
import logging
import tempfile

import biolib.external.execute as execute
from biolib.common import alphanumeric_sort


class Krona():
    """Wrapper for creating Krona plots with KronaTools."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

        execute.check_on_path('ktImportText')

    def create(self, profiles, output_file):
        """Create Krona plot.

        Profiles for multiple items (e.g., genome, metagenome) can
        be specified. The complete hierarchy for each unique element
        should be specified as a semicolon separated string, e.g.,

            k__Bacteria;c__Firmicutes;...;s__

        The number of hits to each unique element is specified in
        the profiles dictionary, e.g.,

            d[unique_id][element_str] = 10

        Parameters
        ----------
        profiles: d[unique_id][element_str] -> count
            Number of hits to specific elements for each item.
        output_file : str
            Name of output file.
        """

        # create temporary files for each item
        cmd = 'ktImportText -o %s' % output_file
        tmp_dir = tempfile.mkdtemp()
        for unique_id in alphanumeric_sort(list(profiles.keys())):
            tmp_file = os.path.join(tmp_dir, unique_id)
            fout = open(tmp_file, 'w')
            for element_str, num_hits in profiles[unique_id].items():
                elements = [x.strip() for x in element_str.split(';')]
                fout.write(str(num_hits) + '\t' + '\t'.join(elements) + '\n')
            fout.close()

            cmd += ' %s,%s' % (tmp_file, unique_id)

        # create krona plot
        execute.run(cmd)

        # clean up temporary files
        shutil.rmtree(tmp_dir)
