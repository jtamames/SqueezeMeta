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

from matplotlib.ticker import NullFormatter, FuncFormatter

from numpy import (zeros_like as np_zeros_like,
                    ones_like as np_ones_like,
                    array as np_array)

def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    return s

class Histogram(AbstractPlot):
    """Histogram plot."""

    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, x, 
                    xlabel, 
                    ylabel, 
                    bins, 
                    color, 
                    yticks=None, 
                    xticks=None, 
                    xlim=None, 
                    invert_xaxis=False,
                    stacked=False,
                    show_legend=False):
        # set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axis = self.fig.add_subplot(111)

        # plot histogram as percentages
        num_entries = sum([len(c) for c in x])
        weights = []
        for c in range(0, len(x)):
            w = [100.0/num_entries for _ in range(0, len(x[c]))]
            weights.append(w)
        
        #weights = []
        #for c in xrange(0, len(x)):
        #    num_entries = len(x[c])
        #    w = [100.0/num_entries for _ in xrange(0, len(x[c]))]
        #    weights.append(w)

        pdf, bins, patches = axis.hist(x, bins=bins, 
                                            #normed=normed,
                                            weights=weights,
                                            color=color, 
                                            lw=0.5, 
                                            histtype='bar', 
                                            stacked=stacked)
        
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        
        #formatter = FuncFormatter(to_percent)
        #axis.yaxis.set_major_formatter(formatter)
        
        if show_legend:
            lgnd = axis.legend(
                               loc='upper left',
                               ncol=1,
                               fontsize=self.options.tick_font_size,
                               handletextpad=0.5,
                               markerscale=2,
                               frameon=False)
        
        # *** Prettify histogram plot
        if yticks:
            axis.set_yticks(yticks)

        if xlim:
            axis.set_xlim(xlim)
        
        if xticks:
            axis.set_xticks(xticks)
            
        if invert_xaxis:
            axis.invert_xaxis()

        self.prettify(axis)

        self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()
