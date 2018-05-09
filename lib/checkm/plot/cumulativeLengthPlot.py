###############################################################################
#
# cumulativeLengthPlot.py - Create a cumulative sequence length plot.
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


class CumulativeLengthPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, fastaFile):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axes = self.fig.add_subplot(111)

        # calculate cumulative sequence length
        seqs = readFasta(fastaFile)

        seqLens = []
        for seq in seqs.values():
            seqLens.append(len(seq))

        seqLens.sort(reverse=True)
        x = np.arange(0, len(seqLens))

        y = []
        cumLen = 0
        for seqLen in seqLens:
            cumLen += seqLen
            y.append(cumLen)

        # Create plot
        axes.plot(x, y, 'k-',)
        axes.set_xlabel('Sequence index')
        axes.set_ylabel('Cumulative sequence length (Mbp)')

        # ensure y-axis include zero
        _, end = axes.get_ylim()
        axes.set_ylim([0, end])

        # Change sequence lengths from bp to kbp
        yticks = axes.get_yticks()
        kbpLabels = []
        for seqLen in yticks:
            label = '%.2f' % (float(seqLen) / 1e6)
            label = label.replace('.00', '')  # remove trailing zeros
            if label[-1] == '0':
                label = label[0:-1]
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
