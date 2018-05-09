###############################################################################
#
# codingDensityPlots.py - Create a GC histogram and a delta-CD plot.
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

from AbstractPlot import AbstractPlot

from checkm.plot.gcPlots import GcPlots
from checkm.plot.codingDensityPlots import CodingDensityPlots
from checkm.plot.tetraDistPlots import TetraDistPlots


class DistributionPlots(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)
        self.options = options

    def plot(self, fastaFile, tetraSigs, distributionsToPlot):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        axesHistGC = self.fig.add_subplot(321)
        axesDeltaGC = self.fig.add_subplot(322)
        axesHistTD = self.fig.add_subplot(323)
        axesDeltaTD = self.fig.add_subplot(324)
        axesHistCD = self.fig.add_subplot(325)
        axesDeltaCD = self.fig.add_subplot(326)

        gcPlots = GcPlots(self.options)
        gcPlots.plotOnAxes(fastaFile, distributionsToPlot, axesHistGC, axesDeltaGC)

        tetraDistPlots = TetraDistPlots(self.options)
        tetraDistPlots.plotOnAxes(fastaFile, tetraSigs, distributionsToPlot, axesHistTD, axesDeltaTD)

        codingDensityPlots = CodingDensityPlots(self.options)
        codingDensityPlots.plotOnAxes(fastaFile, distributionsToPlot, axesHistCD, axesDeltaCD)

        self.fig.tight_layout(pad=1, w_pad=2, h_pad=2)
        self.draw()
