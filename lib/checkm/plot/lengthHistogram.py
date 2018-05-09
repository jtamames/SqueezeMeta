###############################################################################
#
# seqLenPlot.py - Create a sequence length distribution histogram.
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

from matplotlib.ticker import MaxNLocator

from checkm.util.seqUtils import readFasta


class LengthHistogram(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, fastaFile):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)
        axes = self.fig.add_subplot(111)

        # calculate sequence lengths (in kb)
        seqs = readFasta(fastaFile)

        seqLens = []
        for seq in seqs.values():
            seqLens.append(float(len(seq)) / 1e3)

        # set unequal bin sizes (in kb)
        bins = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1e12]
        counts, _edges = np.histogram(seqLens, bins=bins)

        # create histogram
        axes.bar(left=np.arange(0.1, len(counts)), height=counts, width=0.8, color=(0.5, 0.5, 0.5))
        axes.set_xlabel('Sequence length (kbp)')
        axes.set_ylabel('Number sequences (out of %d)' % len(seqs))

        # ensure y-axis include zero
        _, end = axes.get_ylim()
        axes.set_ylim([0, end])
        axes.get_yaxis().set_major_locator(MaxNLocator(integer=True))

        # Change sequence lengths from bp to kbp
        axes.set_xlim([0, len(counts)])
        axes.set_xticks(np.arange(0.5, len(counts)))
        axes.set_xticklabels(['<1', '1-2', '2-5', '5-10', '10-20', '20-50', '50-100', '100-200', '200-500', '>500'])

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
