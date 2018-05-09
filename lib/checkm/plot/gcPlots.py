###############################################################################
#
# gcPlots.py - Create a GC histogram and delta-GC plot.
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

import numpy as np

from AbstractPlot import AbstractPlot

from checkm.binTools import BinTools
from checkm.util.seqUtils import readFasta, baseCount
from checkm.common import findNearest, readDistribution


class GcPlots(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, fastaFile, distributionsToPlot):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesHist = self.fig.add_subplot(121)
        axesDeltaGC = self.fig.add_subplot(122)

        self.plotOnAxes(fastaFile, distributionsToPlot, axesHist, axesDeltaGC)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plotOnAxes(self, fastaFile, distributionsToPlot, axesHist, axesDeltaGC):
        # Read reference distributions from file
        dist = readDistribution('gc_dist')

        # get GC for windows
        seqs = readFasta(fastaFile)

        data = []
        seqLens = []
        for _, seq in seqs.iteritems():
            start = 0
            end = self.options.gc_window_size

            seqLen = len(seq)
            seqLens.append(seqLen)

            while(end < seqLen):
                a, c, g, t = baseCount(seq[start:end])
                try:
                    data.append(float(g + c) / (a + c + g + t))
                except:
                    # it is possible to reach a long stretch of
                    # N's that causes a division by zero error

                    pass

                start = end
                end += self.options.gc_window_size

        if len(data) == 0:
            axesHist.set_xlabel('[Error] No seqs >= %d, the specified window size' % self.options.gc_window_size)
            return

        # Histogram plot
        bins = [0.0]
        binWidth = self.options.gc_bin_width
        binEnd = binWidth
        while binEnd <= 1.0:
            bins.append(binEnd)
            binEnd += binWidth

        axesHist.hist(data, bins=bins, normed=True, color=(0.5, 0.5, 0.5))
        axesHist.set_xlabel('% GC')
        axesHist.set_ylabel('% windows (' + str(self.options.gc_window_size) + ' bp)')

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

        # get GC bin statistics
        binTools = BinTools()
        meanGC, deltaGCs, _ = binTools.gcDist(seqs)

        # Delta-GC vs Sequence length plot
        axesDeltaGC.scatter(deltaGCs, seqLens, c=abs(deltaGCs), s=10, lw=0.5, cmap='gray_r')
        axesDeltaGC.set_xlabel(r'$\Delta$ GC (mean GC = %.1f%%)' % (meanGC * 100))
        axesDeltaGC.set_ylabel('Sequence length (kbp)')

        _, yMaxSeqs = axesDeltaGC.get_ylim()
        xMinSeqs, xMaxSeqs = axesDeltaGC.get_xlim()

        # plot reference distributions
        for distToPlot in distributionsToPlot:
            closestGC = findNearest(np.array(dist.keys()), meanGC)

            # find closest distribution values
            sampleSeqLen = dist[closestGC].keys()[0]
            d = dist[closestGC][sampleSeqLen]
            gcLowerBoundKey = findNearest(d.keys(), (100 - distToPlot) / 2.0)
            gcUpperBoundKey = findNearest(d.keys(), (100 + distToPlot) / 2.0)

            xL = []
            xU = []
            y = []
            for windowSize in dist[closestGC]:
                xL.append(dist[closestGC][windowSize][gcLowerBoundKey])
                xU.append(dist[closestGC][windowSize][gcUpperBoundKey])
                y.append(windowSize)

            # sort by y-values
            sortIndexY = np.argsort(y)
            xL = np.array(xL)[sortIndexY]
            xU = np.array(xU)[sortIndexY]
            y = np.array(y)[sortIndexY]
            axesDeltaGC.plot(xL, y, 'r--', lw=0.5, zorder=0)
            axesDeltaGC.plot(xU, y, 'r--', lw=0.5, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axesDeltaGC.set_ylim([0, yMaxSeqs])

        # ensure x-axis is set appropriately for sequences
        axesDeltaGC.set_xlim([xMinSeqs, xMaxSeqs])

        # draw vertical line at x=0
        axesDeltaGC.vlines(0, 0, yMaxSeqs, linestyle='dashed', color=self.axesColour, zorder=0)

        # Change sequence lengths from bp to kbp
        yticks = axesDeltaGC.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.1f' % (float(seqLen) / 1000)
            label = label.replace('.0', '')  # remove trailing zero
            kbpLabels.append(label)
        axesDeltaGC.set_yticklabels(kbpLabels)

        # Prettify plot
        for a in axesDeltaGC.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axesDeltaGC.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axesDeltaGC.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axesDeltaGC.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axesDeltaGC.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)
