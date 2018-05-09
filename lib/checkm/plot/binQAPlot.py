###############################################################################
#
# binQAPlot.py - Bar plot of bin completeness, contamination,
#                  and strain heterogeneity.
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

import sys
import logging
from operator import itemgetter

import numpy as np
import matplotlib as mpl

from AbstractPlot import AbstractPlot

from checkm.common import binIdFromFilename


class BinQAPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)
        self.logger = logging.getLogger()

    def __sortBinsByCompleteness(self, binFiles, binStatsExt):
        sortedBinIds = []
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)

            sortedBinIds.append([binId, binStatsExt[binId]['Completeness']])

        sortedBinIds.sort(key=itemgetter(1, 0))

        return [x[0] for x in sortedBinIds]

    def plot(self, binFiles, binStatsExt, bIgnoreHetero, aaiHetero):
        self.fig.clear()

        rowPadding = 1.1
        height = rowPadding * self.options.row_height * len(binFiles)  # make room for each row with 5% border

        legendHeight = 0.5 * self.options.row_height + 0.5
        height += legendHeight  # make room for legends
        self.fig.set_size_inches(self.options.width, height)

        if height * self.options.dpi >= 32768:
            errorMsg = '[Error] There are too many bins to plot.'
            errorMsg += '\nThe resulting plot would be %d pixels in height and the maximum allowed size is 32768.' % (height * self.options.dpi)
            errorMsg += '\nPlease reduce the number of bins to be plotted or decrease the DPI (--dpi).'
            self.logger.error(errorMsg)
            return False

        fracRowHeight = self.options.row_height / height

        self.logger.info('  Plotting bin completeness, contamination, and strain heterogeneity.')

        sortedBinIds = self.__sortBinsByCompleteness(binFiles, binStatsExt)
        xLabelBounds = self.xLabelExtents(sortedBinIds, self.options.font_size)
        labelWidth = 1.05 * xLabelBounds.width

        for i, binId in enumerate(sortedBinIds):
            if self.logger.getEffectiveLevel() <= logging.INFO:
                statusStr = '    Plotting bin %d of %d (%.2f%%) bins.' % (i + 1, len(binFiles), float(i + 1) * 100 / len(binFiles))
                sys.stdout.write('%s\r' % statusStr)
                sys.stdout.flush()

            axes = self.fig.add_axes([labelWidth, 1.0 - rowPadding * (i + 1) * fracRowHeight, 1.0 - labelWidth, fracRowHeight])
            # axes = self.fig.add_subplot(len(binFiles), 1, i+1)

            # create a vector for each column indicating:
            #  0) single-copy
            #  1-4) increasing levels of strain heterogeneity
            #  5-8) increasing levels of contamination (2 to 5+ copies)
            #  9) missing
            data = []
            for i in range(binStatsExt[binId]['1']):
                data.append([1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

            if bIgnoreHetero:
                for i in range(binStatsExt[binId]['2']):
                    data.append([0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0])
                for i in range(binStatsExt[binId]['3']):
                    data.append([0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0])
                for i in range(binStatsExt[binId]['4']):
                    data.append([0, 0, 0, 0, 0, 0, 0, 1.0, 0, 0])
                for i in range(binStatsExt[binId]['5+']):
                    data.append([0, 0, 0, 0, 0, 0, 0, 0, 1.0, 0])
            else:
                for geneCountNumber, index in zip(['GCN2', 'GCN3', 'GCN4', 'GCN5+'], [1, 2, 3, 4]):
                    strainHetero = []
                    for markerId in binStatsExt[binId][geneCountNumber]:
                        strainHetero.append(aaiHetero[binId].get(markerId, 0))

                    for sh in sorted(strainHetero, reverse=True):
                        d = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                        d[index] = sh
                        d[index + 4] = 1.0 - sh
                        data.append(d)

            for i in range(binStatsExt[binId]['0']):
                data.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 1.0])

            blues = ['#e0f3f8', '#abd9e9', '#74add1', '#4575b4']  # blue gradient
            reds = ['#fee090', '#fdae61', '#f46d43', '#d73027']  # red gradient

            colors = ['#b2df8a']  # green
            colors += blues
            colors += reds
            colors += ['#777777']  # grey

            self.skylinePlot(binId, axes, data, colors)

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stdout.write('\n')

        # legends
        discreteColourMap = mpl.colors.ListedColormap([colors[0]])
        axisColourMap = self.fig.add_axes([0.05, 0.25 / height, 0.1875, 0.5 * fracRowHeight])

        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=discreteColourMap, norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='horizontal', drawedges=True)
        colourBar.set_ticks([0.5])
        colourBar.set_ticklabels(['1'])
        colourBar.outline.set_linewidth(0.5)
        colourBar.dividers.set_linewidth(0.5)
        axisColourMap.set_title('Single-copy')

        discreteColourMap = mpl.colors.ListedColormap([colors[-1]])
        # axisColourMap = self.fig.add_axes([0.5 - 4.5*(self.options.row_height/self.options.width), self.options.row_height/height, 4*(self.options.row_height/self.options.width), 0.5*self.options.row_height/height])
        axisColourMap = self.fig.add_axes([0.2875, 0.25 / height, 0.1875, 0.5 * fracRowHeight])
        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=discreteColourMap, norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='horizontal', drawedges=True)
        colourBar.set_ticks([0.5])
        colourBar.set_ticklabels(['0'])
        colourBar.outline.set_linewidth(0.5)
        colourBar.dividers.set_linewidth(0.5)
        axisColourMap.set_title('Missing')

        discreteColourMap = mpl.colors.ListedColormap(blues)
        # axisColourMap = self.fig.add_axes([0.5 + 0.5*(self.options.row_height/self.options.width), self.options.row_height/height, 4*(self.options.row_height/self.options.width), 0.5*self.options.row_height/height])
        axisColourMap = self.fig.add_axes([0.525, 0.25 / height, 0.1875, 0.5 * fracRowHeight])
        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=discreteColourMap, norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='horizontal', drawedges=True)
        colourBar.set_ticks([0.125, 0.375, 0.625, 0.875])
        colourBar.set_ticklabels(['2', '3', '4', '5+'])
        colourBar.outline.set_linewidth(0.5)
        colourBar.dividers.set_linewidth(0.5)
        axisColourMap.set_title('Heterogeneity')

        discreteColourMap = mpl.colors.ListedColormap(reds)
        # axisColourMap = self.fig.add_axes([0.5 + 5.5*(self.options.row_height/self.options.width), self.options.row_height/height, 4*(self.options.row_height/self.options.width), 0.5*self.options.row_height/height])
        axisColourMap = self.fig.add_axes([0.7625, 0.25 / height, 0.1875, 0.5 * fracRowHeight])
        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=discreteColourMap, norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='horizontal', drawedges=True)
        colourBar.set_ticks([0.125, 0.375, 0.625, 0.875])
        colourBar.set_ticklabels(['2', '3', '4', '5+'])
        colourBar.outline.set_linewidth(0.5)
        colourBar.dividers.set_linewidth(0.5)
        axisColourMap.set_title('Contamination')

        self.draw()

        return True

    def skylinePlot(self,
                       binId,
                       ax,  # axes to plot onto
                       data,  # data to plot
                       colors  # colors for each level
                       ):

        # make sure this makes sense
        data_copy = np.copy(data).transpose()
        data_shape = np.shape(data_copy)

        # determine the number of bars and corresponding levels from the shape of the data
        num_bars = data_shape[1]
        levels = data_shape[0]
        bar_width = 1.0
        x = np.arange(num_bars) + 0.5 * (1.0 - bar_width)

        # stack the data --
        # replace the value in each level by the cumulative sum of all preceding levels
        data_stack = np.reshape([float(i) for i in np.ravel(np.cumsum(data_copy, axis=0))], data_shape)

        # bars
        ax.bar(x,
               data_stack[0],
               color=colors[0],
               width=bar_width,
               linewidth=0.5,
               align='center'
               )

        for i in np.arange(1, levels):
            ax.bar(x,
                   data_copy[i],
                   bottom=data_stack[i - 1],
                   color=colors[i],
                   width=bar_width,
                   linewidth=0.5,
                   align='center'
                   )

        # limits
        ax.set_xlim(-0.5, num_bars)
        ax.set_ylim(0, np.max(data_stack))

        # labels
        ax.set_yticks(np.arange(0.5, np.max(data_stack)))
        ax.set_yticklabels([binId])

        # ticks
        ax.set_xticks([])
        for a in ax.yaxis.majorTicks:
            a.tick1On = False
            a.tick2On = False

        # borders
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
