###############################################################################
#
# codingDensityPlots.py - Create a CD histogram and a delta-CD plot.
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
import sys

import numpy as np

from AbstractPlot import AbstractPlot

from checkm.prodigal import ProdigalGeneFeatureParser
from checkm.common import readDistribution, findNearest, binIdFromFilename
from checkm.binTools import BinTools
from checkm.util.seqUtils import readFasta, baseCount
from checkm.defaultValues import DefaultValues


class CodingDensityPlots(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, fastaFile, distributionsToPlot):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesHist = self.fig.add_subplot(121)
        axesDeltaCD = self.fig.add_subplot(122)

        self.plotOnAxes(fastaFile, distributionsToPlot, axesHist, axesDeltaCD)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plotOnAxes(self, fastaFile, distributionsToPlot, axesHist, axesDeltaCD):
        # parse Prodigal output
        gffFile = os.path.join(self.options.out_folder, 'bins', binIdFromFilename(fastaFile), DefaultValues.PRODIGAL_GFF)
        if not os.path.exists(gffFile):
            print 'Missing gene feature file (%s). This plot if not compatible with the --genes option.' % DefaultValues.PRODIGAL_GFF
            sys.exit()

        prodigalParser = ProdigalGeneFeatureParser(gffFile)

        # Read reference distributions from file
        dist = readDistribution('cd_dist')

        # get coding density for windows
        seqs = readFasta(fastaFile)

        data = []
        seqLens = []
        for seqId, seq in seqs.iteritems():
            start = 0
            end = self.options.cd_window_size

            seqLen = len(seq)
            seqLens.append(seqLen)

            while(end < seqLen):
                codingBases = prodigalParser.codingBases(seqId, start, end)

                a, c, g, t = baseCount(seq[start:end])
                data.append(float(codingBases) / (a + c + g + t))

                start = end
                end += self.options.cd_window_size

        if len(data) == 0:
            axesHist.set_xlabel('[Error] No seqs >= %d, the specified window size' % self.options.cd_window_size)
            return

        # Histogram plot
        bins = [0.0]
        binWidth = self.options.cd_bin_width
        binEnd = binWidth
        while binEnd <= 1.0:
            bins.append(binEnd)
            binEnd += binWidth

        axesHist.hist(data, bins=bins, normed=True, color=(0.5, 0.5, 0.5))
        axesHist.set_xlabel('% coding density')
        axesHist.set_ylabel('% windows (' + str(self.options.cd_window_size) + ' bp)')

        # Prettify plot
        for a in axesHist.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axesHist.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axesHist.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axesHist.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axesHist.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)

        # get CD bin statistics
        binTools = BinTools()
        meanCD, deltaCDs, _ = binTools.codingDensityDist(seqs, prodigalParser)

        # Delta-CD vs sequence length plot
        axesDeltaCD.scatter(deltaCDs, seqLens, c=abs(deltaCDs), s=10, lw=0.5, cmap='gray_r')
        axesDeltaCD.set_xlabel(r'$\Delta$ CD (mean coding density = %.1f%%)' % (meanCD * 100))
        axesDeltaCD.set_ylabel('Sequence length (kbp)')

        _, yMaxSeqs = axesDeltaCD.get_ylim()
        xMinSeqs, xMaxSeqs = axesDeltaCD.get_xlim()

        # plot reference distributions
        for distToPlot in distributionsToPlot:
            closestCD = findNearest(np.array(dist.keys()), meanCD)

            # find closest distribution values
            sampleSeqLen = dist[closestCD].keys()[0]
            d = dist[closestCD][sampleSeqLen]
            cdLowerBoundKey = findNearest(d.keys(), (100 - distToPlot) / 2.0)
            cdUpperBoundKey = findNearest(d.keys(), (100 + distToPlot) / 2.0)

            xL = []
            xU = []
            y = []
            for windowSize in dist[closestCD]:
                xL.append(dist[closestCD][windowSize][cdLowerBoundKey])
                xU.append(dist[closestCD][windowSize][cdUpperBoundKey])
                y.append(windowSize)

            # sort by y-values
            sortIndexY = np.argsort(y)
            xL = np.array(xL)[sortIndexY]
            xU = np.array(xU)[sortIndexY]
            y = np.array(y)[sortIndexY]
            axesDeltaCD.plot(xL, y, 'r--', lw=0.5, zorder=0)
            axesDeltaCD.plot(xU, y, 'r--', lw=0.5, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axesDeltaCD.set_ylim([0, yMaxSeqs])

        # ensure x-axis is set appropriately for sequences
        axesDeltaCD.set_xlim([xMinSeqs, xMaxSeqs])

        # draw vertical line at x=0
        axesDeltaCD.vlines(0, 0, yMaxSeqs, linestyle='dashed', color=self.axesColour, zorder=0)

        # Change sequence lengths from bp to kbp
        yticks = axesDeltaCD.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.1f' % (float(seqLen) / 1000)
            label = label.replace('.0', '')  # remove trailing zero
            kbpLabels.append(label)
        axesDeltaCD.set_yticklabels(kbpLabels)

        # Prettify plot
        for a in axesDeltaCD.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axesDeltaCD.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axesDeltaCD.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axesDeltaCD.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axesDeltaCD.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)
