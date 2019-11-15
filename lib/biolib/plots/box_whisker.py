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

from matplotlib.ticker import NullFormatter
from matplotlib import collections

from numpy import (sum as np_sum,
                    array as np_array)


class BoxWhisker(AbstractPlot):
    """Box and whisker plot."""

    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, data, data_labels, xlabel, ylabel, ylim=None):
        # set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)

        bp = axis.boxplot(data, sym='r.', whis=[5, 95])
        for box in bp['boxes']:
            box.set( color='black', linewidth=1)
        for whisker in bp['whiskers']:
            whisker.set(color='black', linewidth=1)
        for median in bp['medians']:
            median.set(color='red', linewidth=1)
        for flier in bp['fliers']:
            flier.set(marker='+', color='red')

        axis.set_ylabel(ylabel)
        axis.set_xlabel(xlabel)
        axis.set_xticklabels(data_labels)
        
        axis.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
        axis.set_axisbelow(True)
        
        if ylim:
            axis.set_ylim(ylim)
        
        self.prettify(axis)
        if False:
            for line in axis.yaxis.get_ticklines(): 
                line.set_color(self.axes_colour)
                    
            for line in axis.xaxis.get_ticklines(): 
                line.set_color(self.axes_colour)
                
            for loc, spine in axis.spines.items():
                spine.set_color(self.axes_colour)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()
