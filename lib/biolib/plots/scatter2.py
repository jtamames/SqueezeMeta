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
                    min as np_min,
                    max as np_max,
                    polyfit as np_polyfit,
                    poly1d as np_poly1d)
from scipy.stats import pearsonr


class Scatter2(AbstractPlot):
    """Scatter plot with optional histograms."""

    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, x, y, xlabel, ylabel, show_histograms=False, num_bins=0, xlim=None, ylim=None):
        # set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        # setup histograms
        self.histogram_size = 0.5
        if show_histograms:
                histogram_sizeX = self.histogram_size / self.options.width
                histogram_sizeY = self.histogram_size / self.options.height
        else:
                histogram_sizeX = 0.0
                histogram_sizeY = 0.0

        padding = 0.1   # inches
        xoffset_figspace = (0.4 + padding)/self.options.width
        yoffset_figspace = (0.3 + padding)/self.options.height
        axes_scatter = self.fig.add_axes([xoffset_figspace, 
                                            yoffset_figspace,
                                            1.0 - xoffset_figspace - histogram_sizeX - (2*padding)/self.options.width, 
                                            1.0 - yoffset_figspace - histogram_sizeY - (2*padding)/self.options.height])

        axes_top_histogram = self.fig.add_axes([xoffset_figspace, 
                                                    1.0 - histogram_sizeY - padding/self.options.height,
                                                    1.0 - xoffset_figspace - histogram_sizeX - (2*padding)/self.options.width, 
                                                    histogram_sizeY])

        axes_right_histogram = self.fig.add_axes([1.0 - histogram_sizeX - padding/self.options.width, 
                                                    yoffset_figspace,
                                                    histogram_sizeX, 
                                                    1.0 - yoffset_figspace - histogram_sizeY - (2*padding)/self.options.height])

        # scatter plot
        axes_scatter.scatter(x, y, c='grey', s=10, lw=0.5)
       
        axes_scatter.set_xlabel(xlabel)
        axes_scatter.set_ylabel(ylabel)
        
        if xlim:
            axes_scatter.set_xlim(xlim)
        if ylim:
            axes_scatter.set_ylim(ylim)
            
        # linear regression line
        fit = np_polyfit(x,y,1)
        fit_fn = np_poly1d(fit)
        r, p = pearsonr(x, y)
        axes_scatter.plot(axes_scatter.get_xlim(), fit_fn(axes_scatter.get_xlim()), 'r--', alpha=0.5, zorder=5)
        axes_scatter.text(min(axes_scatter.get_xlim())+0.01, 99.8, "y=%.1fx + %.1f" % (fit[0], fit[1]), fontsize=8)
        axes_scatter.text(min(axes_scatter.get_xlim())+0.01, 99.6, "r=%.2f" % r, fontsize=8)
        
        # plot y=x
        lims = [
            np_min([axes_scatter.get_xlim(), axes_scatter.get_ylim()]),  # min of both axes
            np_max([axes_scatter.get_xlim(), axes_scatter.get_ylim()]),  # max of both axes
        ]

        # now plot both limits against eachother
        print('y=x')
        axes_scatter.plot(lims, lims, 'r-', alpha=0.5, zorder=5)
        
        # *** Prettify scatter plot
        for line in axes_scatter.yaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
                
        for line in axes_scatter.xaxis.get_ticklines(): 
            line.set_color(self.axes_colour)
            
        for loc, spine in axes_scatter.spines.items():
            spine.set_color(self.axes_colour)

        # plot histograms
        if not show_histograms:
            for a in axes_scatter.yaxis.majorTicks:
                    a.tick1On=True
                    a.tick2On=False
                
            for a in axes_scatter.xaxis.majorTicks:
                    a.tick1On=True
                    a.tick2On=False
                    
            for line in axes_scatter.yaxis.get_ticklines(): 
                line.set_color(self.axes_colour)
            
            for line in axes_scatter.xaxis.get_ticklines(): 
                line.set_color(self.axes_colour)

            for loc, spine in axes_scatter.spines.items():
                    if loc in ['right','top']:
                            spine.set_color('none')
                    else:
                        spine.set_color(self.axes_colour)
            
        else: # show histograms 
            # get bin weights for percentage histogram plot
            weights = [100.0/len(x) for _ in range(0, len(x))]
        
            # plot top histogram
            axes_top_histogram.xaxis.set_major_formatter(NullFormatter())
            pdf, bins, patches = axes_top_histogram.hist(x, 
                                                            bins=num_bins, 
                                                            weights=weights,
                                                            color='grey',
                                                            lw=0.5, 
                                                            histtype='bar', 
                                                            stacked=True)
                                                                                           
            max_y = max(pdf)
            axes_top_histogram.yaxis.set_label_position('right')
            axes_top_histogram.set_xlim(axes_scatter.get_xlim())
            axes_top_histogram.set_yticks([int(max_y)])
            axes_top_histogram.set_yticklabels(['%d%%' % int(max_y)])
            axes_top_histogram.set_ylim([0, max_y*1.05])
            
            axes_top_histogram.yaxis.tick_right()
            

            # plot right histogram
            axes_right_histogram.yaxis.set_major_formatter(NullFormatter())
            pdf, bins, patches = axes_right_histogram.hist(y, 
                                                            bins=num_bins, 
                                                            orientation='horizontal',
                                                            weights=weights,
                                                            color='grey',
                                                            lw=0.5, 
                                                            histtype='bar', 
                                                            stacked=True)
                                                            
            max_x = max(pdf)
            axes_right_histogram.set_ylim(axes_scatter.get_ylim())
            axes_right_histogram.set_xticks([int(max_x)])
            axes_right_histogram.set_xticklabels(['%d%%' % int(max_x)])
            axes_right_histogram.set_xlim([0, max_x*1.05])

            # *** Prettify histogram plot
            for a in axes_top_histogram.yaxis.majorTicks:
                    a.tick1On=False
                    a.tick2On=True
                
            for a in axes_top_histogram.xaxis.majorTicks:
                    a.tick1On=True
                    a.tick2On=False
                    
            for line in axes_top_histogram.yaxis.get_ticklines(): 
                line.set_color(self.axes_colour)
            
            for line in axes_top_histogram.xaxis.get_ticklines(): 
                line.set_color(self.axes_colour)

            for loc, spine in axes_top_histogram.spines.items():
                    if loc in ['left','top']:
                            spine.set_color('none')
                    else:
                        spine.set_color(self.axes_colour)

            for a in axes_right_histogram.yaxis.majorTicks:
                    a.tick1On=True
                    a.tick2On=False
                
            for a in axes_right_histogram.xaxis.majorTicks:
                    a.tick1On=True
                    a.tick2On=False
                    
            for line in axes_right_histogram.yaxis.get_ticklines(): 
                line.set_color(self.axes_colour)
            
            for line in axes_right_histogram.xaxis.get_ticklines(): 
                line.set_color(self.axes_colour)

            for loc, spine in axes_right_histogram.spines.items():
                    if loc in ['right','top']:
                            spine.set_color('none') 
                    else:
                        spine.set_color(self.axes_colour)

        #self.fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        self.draw()
