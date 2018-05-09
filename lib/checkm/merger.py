###############################################################################
#
# merger.py - identify bins with complementary sets of marker genes
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
import logging

from checkm.common import checkDirExists
from checkm.resultsParser import ResultsParser


class Merger():
    def __init__(self):
        self.logger = logging.getLogger()

    def run(self, binFiles, outDir, hmmTableFile,
                binIdToModels, binIdToBinMarkerSets,
                minDeltaComp, maxDeltaCont,
                minMergedComp, maxMergedCont):
        checkDirExists(outDir)

        self.logger.info('  Comparing marker sets between all pairs of bins.')

        # ensure all bins are using the same marker set
        markerGenesI = binIdToBinMarkerSets[binIdToBinMarkerSets.keys()[0]].mostSpecificMarkerSet().getMarkerGenes()
        for binIdJ in binIdToBinMarkerSets:
            if markerGenesI != binIdToBinMarkerSets[binIdJ].mostSpecificMarkerSet().getMarkerGenes():
                self.logger.error('  [Error] All bins must use the same marker set to assess potential mergers.')
                sys.exit(0)

        # parse HMM information
        resultsParser = ResultsParser(binIdToModels)

        # get HMM hits to each bin
        resultsParser.parseBinHits(outDir, hmmTableFile)

        # determine union and intersection of marker sets for each pair of bins
        outputFile = os.path.join(outDir, "merger.tsv")
        fout = open(outputFile, 'w')
        fout.write('Bin Id 1\tBin Id 2')
        fout.write('\tBin 1 completeness\tBin 1 contamination')
        fout.write('\tBin 2 completeness\tBin 2 contamination')
        fout.write('\tDelta completeness\tDelta contamination\tMerger delta')
        fout.write('\tMerged completeness\tMerged contamination\n')

        binMarkerHits = resultsParser.results
        binIds = sorted(binMarkerHits.keys())
        for i in xrange(0, len(binMarkerHits)):
            binIdI = binIds[i]

            geneCountsI = binMarkerHits[binIdI].geneCounts(binIdToBinMarkerSets[binIdI].mostSpecificMarkerSet(), binMarkerHits[binIdI].markerHits, True)
            completenessI, contaminationI = geneCountsI[6:8]

            for j in xrange(i + 1, len(binMarkerHits)):
                binIdJ = binIds[j]

                geneCountsJ = binMarkerHits[binIdJ].geneCounts(binIdToBinMarkerSets[binIdJ].mostSpecificMarkerSet(), binMarkerHits[binIdJ].markerHits, True)
                completenessJ, contaminationJ = geneCountsJ[6:8]

                # merge together hits from both bins and calculate completeness and contamination
                mergedHits = {}
                for markerId, hits in binMarkerHits[binIdI].markerHits.iteritems():
                    mergedHits[markerId] = list(hits)

                for markerId, hits in binMarkerHits[binIdJ].markerHits.iteritems():
                    if markerId in mergedHits:
                        mergedHits[markerId].extend(hits)
                    else:
                        mergedHits[markerId] = hits

                geneCountsMerged = binMarkerHits[binIdI].geneCounts(binIdToBinMarkerSets[binIdJ].mostSpecificMarkerSet(), mergedHits, True)
                completenessMerged, contaminationMerged = geneCountsMerged[6:8]

                if not (completenessMerged >= minMergedComp and contaminationMerged < maxMergedCont):
                    continue

                # calculate merged statistics
                deltaComp = completenessMerged - max(completenessI, completenessJ)
                deltaCont = contaminationMerged - max(contaminationI, contaminationJ)
                delta = deltaComp - deltaCont

                if deltaComp >= minDeltaComp and deltaCont < maxDeltaCont:
                    fout.write('%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' %
                                                                        (binIdI, binIdJ,
                                                                         completenessI, contaminationI,
                                                                         completenessJ, contaminationJ,
                                                                         deltaComp, deltaCont, delta,
                                                                         completenessMerged, contaminationMerged))

        fout.close()

        return outputFile
