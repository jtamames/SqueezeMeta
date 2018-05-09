###############################################################################
#
# markerGenePosPlot.py - Create a plot showing the position of marker genes.
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

import os
import math
import operator

import numpy as np

from AbstractPlot import AbstractPlot

from checkm.defaultValues import DefaultValues
from checkm.util.seqUtils import readFasta
from checkm.prodigal import ProdigalFastaParser
from checkm.common import binIdFromFilename

from matplotlib.patches import Rectangle

import matplotlib as mpl


class MarkerGenePosPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def getMarkerGenesPerSeq(self, markerGeneStats):
        markerGenePos = {}
        markerGeneNum = {}

        for geneId in markerGeneStats:
            scaffoldId = geneId[0:geneId.rfind('_')]

            for markerId, hitList in markerGeneStats[geneId].iteritems():
                for hit in hitList:
                    start = hit[0]
                    end = hit[1]
                markerGenePos[scaffoldId] = markerGenePos.get(scaffoldId, []) + [[geneId, markerId, start, end]]
                markerGeneNum[markerId] = len(hitList)

        return markerGenePos, markerGeneNum

    def roundUpToNearest100(self, x):
        return int(math.ceil(x / 100.0)) * 100

    def plot(self, binFile, markerGeneStats, binStats):
        binId = binIdFromFilename(binFile)

        markerGenesPerSeq, _markerGeneNum = self.getMarkerGenesPerSeq(markerGeneStats)

        if len(markerGenesPerSeq) == 0:
            return False

        # Get length of sequences with one or more marker genes
        seqs = readFasta(binFile)
        seqLens = {}
        longestSeq = 0
        binSize = 0
        for seqId, seq in seqs.iteritems():
            seqLen = len(seq)
            binSize += seqLen

            if seqId not in markerGenesPerSeq:
                continue

            seqLens[seqId] = seqLen
            if seqLen > longestSeq:
                longestSeq = seqLen

        sortedSeqLens = sorted(seqLens.iteritems(), key=operator.itemgetter(1), reverse=True)

        MAX_BINS = 100
        plotBinSize = self.roundUpToNearest100(float(longestSeq) / MAX_BINS)
        yLabels = [x[0] for x in sortedSeqLens]

        # get position of genes in bin
        prodigalFastaParser = ProdigalFastaParser()
        geneFile = os.path.join(self.options.out_folder, 'bins', binId, DefaultValues.PRODIGAL_AA)
        genePos = prodigalFastaParser.genePositions(geneFile)

        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        yLabelBounds = self.yLabelExtents(yLabels, self.options.font_size)

        heightBottomLabels = 0.4 + self.options.fig_padding  # inches
        widthSideLabel = yLabelBounds.width * self.options.width + self.options.fig_padding  # inches

        widthPerBin = (self.options.width - widthSideLabel - self.options.fig_padding) / MAX_BINS

        titleHeight = 0.2
        HEIGHT_PER_ROW = 0.2
        height = HEIGHT_PER_ROW * len(sortedSeqLens) + heightBottomLabels + self.options.fig_padding + titleHeight
        rowBinHeight = widthPerBin / HEIGHT_PER_ROW

        self.fig.set_size_inches(self.options.width, height)
        axes = self.fig.add_axes([widthSideLabel / self.options.width, heightBottomLabels / height, \
                                                                        1.0 - (widthSideLabel + self.options.fig_padding) / self.options.width, \
                                                                        1.0 - (heightBottomLabels + self.options.fig_padding + titleHeight) / height])

        # set plot axis
        axes.set_xlim([0, MAX_BINS + 0.1])
        axes.set_xlabel('Position (' + str(plotBinSize) + ' bp/bin)')

        axes.set_ylim([0, len(sortedSeqLens)])
        axes.set_yticks(np.arange(0.5, len(sortedSeqLens) + 0.5, 1.0))

        axes.set_yticklabels(yLabels)

        # legend
        colours = [(1.0, 1.0, 1.0), (127 / 255.0, 201 / 255.0, 127 / 255.0), (255 / 255.0, 192 / 255.0, 134 / 255.0), (190 / 255.0, 174 / 255.0, 212 / 255.0), (0.0, 0.0, 0.0)]
        discreteColourMap = mpl.colors.ListedColormap(colours)
        axisColourMap = self.fig.add_axes([self.options.fig_padding / self.options.width, self.options.fig_padding / height, 0.15, 0.03 * (self.options.width / height)])
        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=discreteColourMap, norm=mpl.colors.Normalize(vmin=0, vmax=1), orientation='horizontal', drawedges=True)
        colourBar.set_ticks([0.1, 0.3, 0.5, 0.7, 0.9])
        colourBar.set_ticklabels(['0', '1', '2', '3', '4+'])
        # colourBar.outline.set_color(self.axesColour)
        colourBar.outline.set_linewidth(0.5)
        # colourBar.dividers.set_color(self.axesColour)
        colourBar.dividers.set_linewidth(0.5)

        for a in axisColourMap.xaxis.majorTicks:
            a.tick1On = False
            a.tick2On = False

        # plot each bin
        binPosX = 0.5
        for seqId, seqLen in sortedSeqLens:
            markerCount = [0] * int(math.ceil(float(seqLen) / plotBinSize))
            for geneId, _markerGeneId, geneStartPos, _geneEndPos in markerGenesPerSeq[seqId]:
                binPos = int(float(genePos[geneId][0] + geneStartPos) / plotBinSize)
                markerCount[binPos] += 1

            for i in xrange(0, len(markerCount)):
                if markerCount[i] < len(colours):
                    axes.add_patch(Rectangle((i + 0.1, binPosX - 0.4 * rowBinHeight), 0.8, 0.8 * rowBinHeight, facecolor=colours[markerCount[i]], lw=0.2))
                else:
                    axes.add_patch(Rectangle((i + 0.1, binPosX - 0.4 * rowBinHeight), 0.8, 0.8 * rowBinHeight, facecolor=colours[-1], lw=0.2))

            binPosX += 1.0

        # set plot title
        titleStr = binId + '\n'
        titleStr += '(%.2f Mbp, %d seqs, %.2f%% complete, %.2f%% contamination)' % (float(binSize) / 1e6, len(seqs), binStats['Completeness'], binStats['Contamination'])
        axes.set_title(titleStr)

        # Prettify plot
        for a in axes.yaxis.majorTicks:
            a.tick1On = False
            a.tick2On = False

        for a in axes.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axes.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axes.xaxis.get_ticklines():
            line.set_color(self.axesColour)
            line.set_ms(2)

        for loc, spine in axes.spines.iteritems():
            if loc in ['left', 'right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)

        self.draw()

        return True
