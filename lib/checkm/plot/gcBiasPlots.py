###############################################################################
#
# gcBiasPlots.py - Create a pair of plots used to explore GC bias on coverage.
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

from checkm.util.seqUtils import readFasta, baseCount

from numpy import mean, array, log, poly1d, polyfit


class GcBiasPlot(AbstractPlot):
    def __init__(self, options):
        AbstractPlot.__init__(self, options)

    def plot(self, binFile, coverageProfile):
        # Set size of figure
        self.fig.clear()
        self.fig.set_size_inches(self.options.width, self.options.height)

        windowAxes = self.fig.add_subplot(121)
        seqAxes = self.fig.add_subplot(122)

        self.plotOnAxes(binFile, coverageProfile, windowAxes, seqAxes)

        self.fig.tight_layout(pad=1)
        self.draw()

    def plotOnAxes(self, binFile, coverageProfile, windowAxes, seqAxes):

        # get GC for windows
        seqs = readFasta(binFile)

        gcProfile = {}
        for seqId, seq in seqs.iteritems():
            start = 0
            end = self.options.window_size

            windowGCs = []
            while(end < len(seq)):
                a, c, g, t = baseCount(seq[start:end])
                windowGCs.append(float(g + c) / (a + c + g + t))

                start = end
                end += self.options.window_size

            a, c, g, t = baseCount(seq)
            seqGC = float(g + c) / (a + c + g + t)
            gcProfile[seqId] = [seqGC, windowGCs]

        # plot GC vs coverage for windows
        gc = []
        coverage = []
        for seqId, gcInfo in gcProfile.iteritems():
            gc += gcInfo[1]
            coverage += coverageProfile[seqId][1]

        windowAxes.scatter(gc, coverage, c=abs(array(coverage)), s=10, lw=0.5, cmap='gray_r')
        windowAxes.set_xlabel('GC (mean = %.1f%%)' % (mean(gc) * 100))
        windowAxes.set_ylabel('Coverage (mean = %.1f)' % mean(coverage))

        # plot linear regression line
        if len(gc) > 1:
            slope, inter = polyfit(gc, coverage, 1)
            fit_fn = poly1d([slope, inter])  # fit_fn is now a function which takes in x and returns an estimate for y
            windowAxes.plot([min(gc), max(gc)], fit_fn([min(gc), max(gc)]), '--r', lw=0.5)
            windowAxes.set_title('GC vs. Coverage\n(window size = %d bp, slope = %.2f)' % (self.options.window_size, slope))
        else:
            # not possible to calculate best fit line
            windowAxes.set_title('GC vs. Coverage\n(window size = %d bp, no best fit line)' % self.options.window_size)

        # Prettify plot
        for a in windowAxes.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in windowAxes.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in windowAxes.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in windowAxes.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in windowAxes.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)

        # plot GC vs coverage for entire sequences
        gc = []
        coverage = []
        seqLen = []
        for seqId, gcInfo in gcProfile.iteritems():
            gc.append(gcInfo[0])
            coverage.append(coverageProfile[seqId][0])
            seqLen.append(len(seqs[seqId]))

        # set marker size proportional to sequence length
        markerSize = log(array(seqLen))  # log-scale
        markerSize = (markerSize - min(markerSize)) / max(markerSize)  # normalize between 0 and 1
        markerSize = markerSize * 200 + 10  # normalize between 10 and 200

        seqAxes.scatter(gc, coverage, c=abs(array(coverage)), s=markerSize, lw=0.5, cmap='gray_r')
        seqAxes.set_xlabel('GC (mean = %.1f%%)' % (mean(gc) * 100))
        seqAxes.set_ylabel('Coverage (mean = %.1f)' % mean(coverage))
        seqAxes.set_title('GC vs. Coverage\nIndividual Sequences')

        # Prettify plot
        for a in seqAxes.yaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for a in seqAxes.xaxis.majorTicks:
            a.tick1On = True
            a.tick2On = False

        for line in seqAxes.yaxis.get_ticklines():
            line.set_color(self.axesColour)

        for line in seqAxes.xaxis.get_ticklines():
            line.set_color(self.axesColour)

        for loc, spine in seqAxes.spines.iteritems():
            if loc in ['right', 'top']:
                spine.set_color('none')
            else:
                spine.set_color(self.axesColour)
