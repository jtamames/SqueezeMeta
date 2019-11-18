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



import re
from itertools import permutations

from numpy import (array as np_array,
                    zeros as np_zeros)

import scipy.cluster.hierarchy as cluster
import scipy.spatial.distance as dist

import matplotlib.pyplot as pylab
from matplotlib.colors import ListedColormap

from biolib.plots.abstract_plot import AbstractPlot
from biolib.common import alphanumeric_sort

# import mpld3
# from comparem.plots.mpld3_plugins import Tooltip


class Heatmap(AbstractPlot):
    def __init__(self, infile, outfile):
        AbstractPlot.__init__(self, None)
        
        self.outfile = outfile
        self.genomes = None
        self._parse_data(infile)
        
        self.colormap = pylab.cm.bwr

        self.discreteColourMap = ListedColormap([(141/255.0, 211/255.0, 199/255.0),(255/255.0, 255/255.0, 179/255.0),
                                                    (190/255.0, 186/255.0, 218/255.0),(251/255.0, 128/255.0, 114/255.0),
                                                    (128/255.0, 177/255.0, 211/255.0),(253/255.0, 180/255.0, 98/255.0),
                                                    (179/255.0, 222/255.0, 105/255.0),(252/255.0, 205/255.0, 229/255.0),
                                                    (217/255.0, 217/255.0, 217/255.0), (188/255.0, 128/255.0, 189/255.0),
                                                    (204/255.0, 235/255.0, 197/255.0),(255/255.0, 237/255.0, 111/255.0)])
        
    def _parse_data(self, infile):
        data = {}
        with open(infile) as fp:
            fp.readline()
            genomes = set()
            for line in fp:
                fields = line.rstrip().split('\t')
                fields[0] = re.sub(r'_genes$', "", fields[0])
                fields[2] = re.sub(r'_genes$', "", fields[2])
                genomes.add(fields[0])
                genomes.add(fields[2])
                try:
                    data[fields[0]][fields[2]] = [float(fields[5]), float(fields[7])]
                except KeyError:
                    data[fields[0]] = {}
                    data[fields[0]][fields[2]] = [float(fields[5]), float(fields[7])]
                except IndexError as e:
                    print(fields)
                    raise e

        self.perc_ids = np_zeros([len(genomes), len(genomes)])
        self.perc_aln = np_zeros([len(genomes), len(genomes)])
        genome_to_index = {}
        self.genomes = [None] * len(genomes)
        for n, g in enumerate(alphanumeric_sort(genomes)):
            genome_to_index[g] = n
            self.genomes[n] = g

        self.genomes = np_array(self.genomes)
        for g1, g2 in permutations(genomes, 2):
            try:
                self.perc_ids[genome_to_index[g1]][genome_to_index[g2]] = 100.0 - data[g1][g2][0]
                self.perc_aln[genome_to_index[g1], genome_to_index[g2]] = data[g1][g2][1]
            except:
                self.perc_ids[genome_to_index[g1]][genome_to_index[g2]] = 100.0 - data[g2][g1][0]
                self.perc_aln[genome_to_index[g1], genome_to_index[g2]] = data[g2][g1][1]
        
    def plotDendrogram(self, matrix, axis, clusteringThreshold, orientation):

        d = dist.pdist(matrix)
        linkage = cluster.linkage(dist.squareform(d), method='average', metric='cityblock')
        dendrogram = cluster.dendrogram(linkage, orientation=orientation, link_color_func=lambda k: 'k')
        index = cluster.fcluster(linkage, clusteringThreshold * max(linkage[:,2]), 'distance')
        axis.set_xticks([])
        axis.set_yticks([])

        return index, dendrogram['leaves']

    def plot(self, cluster=False, method='average', metric='euclidean'):
        matrix = np_array(self.perc_ids)
        
        rowHeaders = self.genomes
        clusteringThreshold = 40

        # setup figure
        self.fig.clear()
        self.fig.set_size_inches(0.4 * len(self.genomes) + 2, 0.4 * len(self.genomes) + 2)

        # position all figure elements
        colourBarWidth = 0.015
        heatmapMargin = 0.005

        rowDendrogramX = 0.05
        rowDendrogramY = 0.2
        rowDendrogramW = 0.15
        rowDendrogramH = 0.6

        rowClusterBarX = rowDendrogramX + rowDendrogramW + heatmapMargin
        rowClusterBarY = rowDendrogramY
        rowClusterBarW = colourBarWidth
        rowClusterBarH = rowDendrogramH

        colDendrogramX = rowClusterBarX + rowClusterBarW + heatmapMargin
        colDendrogramY = rowDendrogramY + rowDendrogramH + heatmapMargin + colourBarWidth + heatmapMargin
        colDendrogramW = rowDendrogramH
        colDendrogramH = rowDendrogramW

        colClusterBarX = colDendrogramX
        colClusterBarY = rowDendrogramY + rowDendrogramH + heatmapMargin
        colClusterBarW = colDendrogramW
        colClusterBarH = colourBarWidth

        heatmapX = rowClusterBarX + rowClusterBarW + heatmapMargin
        heatmapY = rowDendrogramY
        heatmapW = colDendrogramW
        heatmapH = rowDendrogramH

        legendX = rowDendrogramX
        legendY = 0.9
        legendW = rowDendrogramW + heatmapMargin + colourBarWidth
        legendH = 0.05

        # plot dendrograms
        axisRowDendrogram = self.fig.add_axes([rowDendrogramX, rowDendrogramY, rowDendrogramW, rowDendrogramH], frame_on=False)
        ind1, leafIndex1 = self.plotDendrogram(matrix, axisRowDendrogram, clusteringThreshold, 'right')

        axisColDendrogram = self.fig.add_axes([colDendrogramX, colDendrogramY, colDendrogramW, colDendrogramH], frame_on=False)
        ind2, leafIndex2 = self.plotDendrogram(matrix.T, axisColDendrogram, clusteringThreshold, 'top')

        # plot column clustering bars
        matrix = matrix[:,leafIndex2]
        ind2 = ind2[:,leafIndex2]

        axc = self.fig.add_axes([colClusterBarX, colClusterBarY, colClusterBarW, colClusterBarH])  # axes for column side colorbar
        dc = numpy.array(ind2, dtype=int)
        dc.shape = (1,len(ind2))
        im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=self.discreteColourMap)
        axc.set_xticks([])
        axc.set_yticks([])

        # plot row clustering bars
        matrix = matrix[leafIndex1,:]
        ind1 = ind1[leafIndex1,:]

        axr = self.fig.add_axes([rowClusterBarX, rowClusterBarY, rowClusterBarW, rowClusterBarH])
        dr = numpy.array(ind1, dtype=int)
        dr.shape = (len(ind1),1)
        im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=self.discreteColourMap)
        axr.set_xticks([])
        axr.set_yticks([])

        # determine scale for colour map
        minValue = 1e6
        maxValue = 0
        for row in matrix:
            minValue = min(minValue, min(row))
            maxValue = max(maxValue, max(row))
        norm = mpl.colors.Normalize(minValue, maxValue)

        # plot heatmap
        axisHeatmap = self.fig.add_axes([heatmapX, heatmapY, heatmapW, heatmapH])
        axisHeatmap.matshow(matrix, aspect='auto', origin='lower', cmap = self.colormap, norm = norm)
        axisHeatmap.set_xticks([])
        axisHeatmap.set_yticks([])

        # row and column labels
        for i in range(0, len(rowHeaders)):
            axisHeatmap.text(matrix.shape[1] - 0.5, i, '  ' + rowHeaders[leafIndex1[i]], horizontalalignment="left")

        for i in range(0, len(colHeaders)):
            axisHeatmap.text(i, -0.5, '  ' + colHeaders[leafIndex2[i]], rotation = 270, verticalalignment="top")

        # plot colour map legend
        axisColourMap = self.fig.add_axes([legendX, legendY, legendW, legendH], frame_on=False)  # axes for colorbar
        colourBar = mpl.colorbar.ColorbarBase(axisColourMap, cmap=self.colormap, norm=norm, orientation='horizontal')
        #axisColourMap.set_title("Relative Abundance")
        colourBar.set_ticks([minValue, 0.5*(maxValue-minValue) + minValue, maxValue])
        colourBar.set_ticklabels(['%.1f' % (minValue*100.0) + '%', '%.1f' % ((0.5*(maxValue-minValue) + minValue)*100.0) + '%', '%.1f' % (maxValue*100.0) + '%'])

        for i in range(0, len(rowHeaders)):
            axisHeatmap.plot([-0.5, len(colHeaders)-0.5], [i-0.5,i-0.5], color='white', linestyle='-', linewidth=1)

        for i in range(0, len(colHeaders)):
            axisHeatmap.plot([i-0.5, i-0.5], [-0.5, len(rowHeaders)-0.5], color='white', linestyle='-', linewidth=1)

        self.fig.tight_layout()
        self.fig.savefig(self.outfile)
