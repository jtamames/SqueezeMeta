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

from checkm.util.seqUtils import readFasta
from checkm.common import readDistribution, findNearest
from checkm.genomicSignatures import GenomicSignatures
from checkm.binTools import BinTools


class TetraDistPlots(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, fastaFile, tetraSigs, distributionsToPlot):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesHist = self.fig.add_subplot(121)
        axesDeltaTD = self.fig.add_subplot(122)

        self.plotOnAxes(fastaFile, tetraSigs, distributionsToPlot, axesHist, axesDeltaTD)

        self.fig.tight_layout(pad=1, w_pad=1)
        self.draw()

    def plotOnAxes(self, fastaFile, tetraSigs, distributionsToPlot, axesHist, axesDeltaTD):
        # Read reference distributions from file
        dist = readDistribution('td_dist')

        # get tetranucleotide signature for bin
        seqs = readFasta(fastaFile)

        binTools = BinTools()
        binSig = binTools.binTetraSig(seqs, tetraSigs)

        # get tetranucleotide distances for windows
        genomicSig = GenomicSignatures(K=4, threads=1)

        data = []
        seqLens = []
        deltaTDs = []
        for seqId, seq in seqs.iteritems():
            start = 0
            end = self.options.td_window_size

            seqLen = len(seq)
            seqLens.append(seqLen)
            deltaTDs.append(genomicSig.distance(tetraSigs[seqId], binSig))

            while(end < seqLen):
                windowSig = genomicSig.seqSignature(seq[start:end])
                data.append(genomicSig.distance(windowSig, binSig))

                start = end
                end += self.options.td_window_size

        if len(data) == 0:
            axesHist.set_xlabel('[Error] No seqs >= %d, the specified window size' % self.options.td_window_size)
            return

        deltaTDs = np.array(deltaTDs)

        # Histogram plot
        bins = [0.0]
        binWidth = self.options.td_bin_width
        binEnd = binWidth
        while binEnd <= 1.0:
            bins.append(binEnd)
            binEnd += binWidth

        axesHist.hist(data, bins=bins, normed=True, color=(0.5, 0.5, 0.5))
        axesHist.set_xlabel(r'$\Delta$ TD')
        axesHist.set_ylabel('% windows (' + str(self.options.td_window_size) + ' bp)')

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
        meanTD, deltaTDs = binTools.tetraDiffDist(seqs, genomicSig, tetraSigs, binSig)

        # Delta-TD vs Sequence length plot
        axesDeltaTD.scatter(deltaTDs, seqLens, c=abs(deltaTDs), s=10, lw=0.5, cmap='gray_r')
        axesDeltaTD.set_xlabel(r'$\Delta$ TD (mean TD = %.2f)' % meanTD)
        axesDeltaTD.set_ylabel('Sequence length (kbp)')

        _, yMaxSeqs = axesDeltaTD.get_ylim()
        xMinSeqs, xMaxSeqs = axesDeltaTD.get_xlim()

        # plot reference distributions
        for distToPlot in distributionsToPlot:
            boundKey = findNearest(dist[dist.keys()[0]].keys(), distToPlot)

            x = []
            y = []
            for windowSize in dist:
                x.append(dist[windowSize][boundKey])
                y.append(windowSize)

            # sort by y-values
            sortIndexY = np.argsort(y)
            x = np.array(x)[sortIndexY]
            y = np.array(y)[sortIndexY]

            # make sure x-values are strictly decreasing as y increases
            # as this is conservative and visually satisfying
            for i in xrange(0, len(x) - 1):
                for j in xrange(i + 1, len(x)):
                    if x[j] > x[i]:
                        if j == len(x) - 1:
                            x[j] = x[i]
                        else:
                            x[j] = (x[j - 1] + x[j + 1]) / 2  # interpolate values from neighbours

                        if x[j] > x[i]:
                            x[j] = x[i]

            axesDeltaTD.plot(x, y, 'r--', lw=0.5, zorder=0)

        # ensure y-axis include zero and covers all sequences
        axesDeltaTD.set_ylim([0, yMaxSeqs])

        # ensure x-axis is set appropriately for sequences
        axesDeltaTD.set_xlim([xMinSeqs, xMaxSeqs])

        # draw vertical line at x=0
        axesDeltaTD.vlines(0, 0, yMaxSeqs, linestyle='dashed', color=self.axesColour, zorder=0)

        # Change sequence lengths from bp to kbp
        yticks = axesDeltaTD.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.1f' % (float(seqLen) / 1000)
            label = label.replace('.0', '')  # remove trailing zero
            kbpLabels.append(label)
        axesDeltaTD.set_yticklabels(kbpLabels)

        # Prettify plot
        for a in axesDeltaTD.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axesDeltaTD.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axesDeltaTD.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axesDeltaTD.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axesDeltaTD.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)
