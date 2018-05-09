###############################################################################
#
# parallelCoordPlot.py - Create a parallel coordinate plot.
#
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

from matplotlib.collections import LineCollection
import matplotlib.ticker as ticker
from matplotlib.cm import get_cmap

from AbstractPlot import AbstractPlot


class ParallelCoordPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def createColorMapGC(self):
        return get_cmap('RdYlBu')

    def plot(self, binIdToHighlight, seqStats, coverageStats):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        # create data points for each sequence
        data = []
        xlabels = ['GC']
        bHighlight = []
        for binId in seqStats:
            for seqId in seqStats[binId]:
                data.append([seqStats[binId][seqId]['GC']] + coverageStats[binId][seqId].values())

                if binId == binIdToHighlight:
                    bHighlight.append(True)
                else:
                    bHighlight.append(False)

                if len(xlabels) == 1:
                    for bamId in coverageStats[binId][seqId]:
                        xlabels.append(bamId)

        dims = len(data[0])
        x = range(dims)
        axes = []
        for i in xrange(1, dims):
            a = self.fig.add_subplot(1, dims - 1, i)
            axes.append(a)

        # Calculate the limits on the data
        min_max_range = list()
        for m in zip(*data):
            mn = min(m)
            mx = max(m)
            if mn == mx:
                mn -= 0.5
                mx = mn + 1.
            r = float(mx - mn)
            min_max_range.append((mn, mx, r))

        # Normalize the data sets
        norm_data_sets = list()
        for ds in data:
            nds = [(value - min_max_range[dimension][0]) /
                    min_max_range[dimension][2]
                    for dimension, value in enumerate(ds)]
            norm_data_sets.append(nds)
        data = norm_data_sets

        for ax in axes:
            ax.set_ylim(0, 1)

        # Plot the datasets on all the subplots
        colourMapGC = self.createColorMapGC()
        for i, ax in enumerate(axes):
            colourList = []
            lines = []

            highlightColourList = []
            highlightLines = []
            for dsi, d in enumerate(data):
                if bHighlight[dsi]:
                    highlightLines.append(((x[0], d[0]), (x[1], d[1])))
                    highlightColourList.append(colourMapGC(d[0]))
                else:
                    lines.append(((x[0], d[0]), (x[1], d[1])))
                    colourList.append((0.9, 0.9, 0.9, 0.1))

            ax.add_collection(LineCollection(lines, colors=colourList, zorder=0))
            ax.add_collection(LineCollection(highlightLines, colors=highlightColourList, zorder=1))

            ax.set_xlim([x[i], x[i + 1]])

        # Set the y axis labels
        for dimension, (axx, xx) in enumerate(zip(axes, x[:-1])):
            axx.xaxis.set_major_locator(ticker.FixedLocator([xx]))
            ticks = len(axx.get_yticklabels())
            labels = list()
            step = min_max_range[dimension][2] / (ticks - 1)
            mn = min_max_range[dimension][0]
            for i in xrange(ticks):
                v = mn + i * step
                labels.append('%4.2f' % v)

            axx.set_yticklabels(labels)
            axx.set_xticklabels(xlabels[dimension:dimension + 1])
            for label in axx.get_xticklabels():
                label.set_rotation(90)

            for loc, spine in axx.spines.iteritems():
                if loc in ['bottom', 'top']:
                    spine.set_color('none')
                else:
                    spine.set_color(self.axesColour)

        # Move the final axis' ticks to the right-hand side
        axx = axes[-1].twinx()
        dimension += 1
        axx.xaxis.set_major_locator(ticker.FixedLocator([x[-2], x[-1]]))
        ticks = len(axx.get_yticklabels())
        step = min_max_range[dimension][2] / (ticks - 1)
        mn = min_max_range[dimension][0]
        labels = ['%4.2f' % (mn + i * step) for i in xrange(ticks)]

        axx.set_yticklabels(labels)
        for tick in axx.yaxis.get_major_ticks():
            tick.set_pad(-6)
            tick.label2.set_horizontalalignment('right')

        axx.set_xticklabels(xlabels[-2:])

        # Stack the subplots
        self.fig.tight_layout(pad=1)
        self.fig.subplots_adjust(wspace=0)

        self.draw()
