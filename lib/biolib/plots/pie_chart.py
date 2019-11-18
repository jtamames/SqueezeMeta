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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2015'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

from .abstract_plot import AbstractPlot
from .color_brewer import ColorBrewer


class PieChart(AbstractPlot):
    """Generate pie chart."""

    def __init__(self, options=None):
        """Initialization."""
        AbstractPlot.__init__(self, options)

    def plot(self):  # sizes, labels, min_size=1, other_label='Other'):
        # set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axes = self.fig.add_subplot(111, aspect='equal')

        # create plot
        labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
        sizes = [15, 30, 45, 10]
        colours = list(ColorBrewer().maps['qualPastel1'].values())[0:len(sizes)]

        axes.pie(sizes,
                 labels=labels,
                 colors=colours,
                 autopct='%1.1f%%',
                 shadow=False,
                 startangle=90)

        self.draw()
