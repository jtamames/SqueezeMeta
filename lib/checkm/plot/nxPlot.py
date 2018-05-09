###############################################################################
#
# nxPlot.py - Create a Nx-plot.
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


class NxPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def calculateNx(self, x, seqs):
        seqLens = []
        for _, seq in seqs.iteritems():
            seqLens.append(len(seq))

        sumSeqLens = sum(seqLens)

        seqLens.sort(reverse=True)
        x.sort()

        testSum = 0
        xIndex = 0
        nxThreshold = x[xIndex] * sumSeqLens
        nx = []
        for seqLen in seqLens:
            testSum += seqLen

            while(testSum >= nxThreshold):
                nx.append(seqLen)
                if len(nx) == len(x):
                    break

                xIndex += 1
                nxThreshold = x[xIndex] * sumSeqLens

        return nx

    def plot(self, fastaFile):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axes = self.fig.add_subplot(111)

        # calculate Nx
        seqs = readFasta(fastaFile)
        x = np.arange(0, 1.0 + 0.5 * self.options.step_size, self.options.step_size)
        nx = self.calculateNx(x, seqs)

        # Create plot
        axes.plot(x, nx, 'ko-', ms=4)
        axes.set_xlabel('Nx')
        axes.set_ylabel('Sequence length (kbp)')

        # Change sequence lengths from bp to kbp
        yticks = axes.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.1f' % (float(seqLen) / 1000)
            label = label.replace('.0', '')  # remove trailing zero
            kbpLabels.append(label)
        axes.set_yticklabels(kbpLabels)

        # Prettify plot
        for a in axes.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in axes.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in axes.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in axes.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in axes.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)

        self.fig.tight_layout(pad=1)
        self.draw()
