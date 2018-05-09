###############################################################################
#
# binComparer.py - compare two sets of bins
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

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2014"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = ""

import logging
import math
import csv

from common import binIdFromFilename
from checkm.util.seqUtils import readFasta


class UnionBin:
    def __init__(self, binningIndex, completeness, contamination, binFile):
        self.binningIndex = binningIndex
        self.completeness = completeness
        self.contamination = contamination
        self.binId = binIdFromFilename(binFile)
        self.seqs = readFasta(binFile)
        self.binFile = binFile

    def numBasesOverlapping(self, anotherUnionBin):
        commonContigs = set(self.seqs.keys()).intersection(anotherUnionBin.seqs.keys())
        return sum([len(self.seqs[contig]) for contig in commonContigs])

    def numBases(self):
        return sum([len(seq) for seq in self.seqs.values()])

    def compContSquaredScored(self):
        # The divide by 100s are not necessary for direct comparison, but have the
        # advantage of a perfect bin having a score of 1 and terrible bins 0.
        return self.completeness / 100 * math.pow(1 - (self.contamination / 100), 2)


class UnionCheckmQaTsv:
    def __init__(self, tsvFile):
        csv.register_dialect('tabffs', delimiter='\t', quoting=csv.QUOTE_NONE)
        self.binIdToCompleteness = {}
        self.binIdToContamination = {}
        with open(tsvFile) as f:
            for row in csv.DictReader(f, dialect='tabffs'):
                binId = row['Bin Id']
                self.binIdToCompleteness[binId] = float(row['Completeness'])
                self.binIdToContamination[binId] = float(row['Contamination'])

    def completeness(self, binId):
        return self.binIdToCompleteness[binId]

    def contamination(self, binId):
        return self.binIdToContamination[binId]


class BinUnion(object):
    def __init__(self):
        self.logger = logging.getLogger()

    def report(self, binFolders, binFileSets, checkmQaTsvs, unionBinOutputFile, multiplyBinnedContigsFile, minCompleteness=0, maxContamination=0):
        # Read QA files
        qas = [UnionCheckmQaTsv(f) for f in checkmQaTsvs]

        bestCandidates = self.getBestCandidates(binFileSets, qas, minCompleteness, maxContamination)

        numBinsBinningWise = [0] * len(binFolders)
        for candidate in bestCandidates:
            numBinsBinningWise[candidate.binningIndex] += 1
        self.logger.info("")
        for binIndex, numBins in enumerate(numBinsBinningWise):
            self.logger.info("   Kept %i out of %i bins from %s" % (numBins,
                                                                len(binFileSets[binIndex]),
                                                                binFolders[binIndex]
                                                                ))

        self.logger.info("")
        with open(multiplyBinnedContigsFile, 'w') as multiplyBinnedOutput:
            self.printMultiplyBinnedContigs(bestCandidates, multiplyBinnedOutput)

        with open(unionBinOutputFile, 'w') as out:
            for candidate in bestCandidates:
                out.write(candidate.binFile)
                out.write("\n")

        self.logger.info("")
        self.logger.info("   Wrote %i bins to %s" % (len(bestCandidates),
                                                  unionBinOutputFile,
                                                  ))

    def getBestCandidates(self, binFileSets, qas, minCompleteness, maxContamination):

        # Take the first set of bins as the best set yet
        bestCandidates = []
        for f in binFileSets[0]:
            if qas[0].completeness(binIdFromFilename(f)) >= minCompleteness and qas[0].contamination(binIdFromFilename(f)) <= maxContamination:
                bestCandidates.append(UnionBin(0,
                                       qas[0].completeness(binIdFromFilename(f)),
                                       qas[0].contamination(binIdFromFilename(f)),
                                       f
                                       ))

        # For each bin in the second or after set,
        for binningIndex, binFileSet in enumerate(binFileSets):
            if binningIndex == 0:
                continue

            currentRoundCandidatesToAdd = []
            for binFile in binFileSet:
                # Is it >50% (by sequence) aligned with any of the bins in the best set?
                binId = binIdFromFilename(binFile)

                if qas[binningIndex].completeness(binId) >= minCompleteness and qas[binningIndex].contamination(binId) <= maxContamination:

                    current = UnionBin(binningIndex,
                                       qas[binningIndex].completeness(binId),
                                       qas[binningIndex].contamination(binId),
                                       binFile)

                    fiftyPercent = 0.5 * current.numBases()
                    accountedFor = False
                    for i, bestBin in enumerate(bestCandidates):
                        overlap = current.numBasesOverlapping(bestBin)
                        fiftyPercentBest = 0.5 * bestBin.numBases()

                        if overlap > fiftyPercent or overlap > fiftyPercentBest:
                            self.logger.debug("Comparing best bin %s and current bin %s, overlap is %i" % (bestBin.binId, current.binId, overlap))

                        if overlap > fiftyPercent and overlap > fiftyPercentBest:
                            accountedFor = True
                            # Choose the best one
                            if current.compContSquaredScored() > bestBin.compContSquaredScored():
                                self.logger.debug("The newly found bin is better, going with that")
                                # Found a better one, replace the best bin with that
                                bestCandidates[i] = current
                                # There's a bug here, but is sufficiently rare and hard to fix that meh. If a multiple bins have
                                # the same contig, then it is possible that a bin can overlap > 50% with more
                                # than one bin. So by breaking out of this for loop we may not be replacing the 'optimal'
                                # bin. But then should it be a 1:1 swap anyway? meh.
                                break
                        elif overlap > fiftyPercent or overlap > fiftyPercentBest:
                            self.logger.warn("Bins %s and %s with sizes %i and %i overlap by %i bases and so have unusual overlap ratios, proceeding as if they are distinct bins" % (
                                                                                                                                                          bestBin.binId,
                                                                                                                                                          current.binId,
                                                                                                                                                          bestBin.numBases(),
                                                                                                                                                          current.numBases(),
                                                                                                                                                          overlap
                                                                                                                                                          ))
                            # Bins don't overlap, continue to go through the loop again

                    if not accountedFor:
                        currentRoundCandidatesToAdd.append(current)

            # Add all the bins that hit no other bins to the bestCandidates list
            # Do this after so that bins are not compared to themselves (saves some time?)
            for b in currentRoundCandidatesToAdd:
                self.logger.debug("Adding unmatched bin %s from %s" % (b.binId, b.binningIndex))
                bestCandidates.append(b)

        return bestCandidates

    def printMultiplyBinnedContigs(self, bestCandidates, multiplyBinnedOutput):
        contigToBin = {}
        for binn in bestCandidates:
            for contigName in binn.seqs.keys():
                if contigName in contigToBin:
                    contigToBin[contigName].append(binn)
                else:
                    contigToBin[contigName] = [binn]

        # IPython.embed()
        numMultiplyBinnedContigs = 0
        multiplyBinnedContigsLength = 0
        for contigName, binList in contigToBin.iteritems():
            if len(binList) > 1:
                numMultiplyBinnedContigs += 1
                multiplyBinnedContigsLength += len(binList[0].seqs[contigName])

                for binn in binList:
                    multiplyBinnedOutput.write("\t".join((contigName,
                                                           str(len(binList[0].seqs[contigName])),
                                                           binn.binFile
                                                           )))
                    multiplyBinnedOutput.write("\n")

        if numMultiplyBinnedContigs != 0:
            self.logger.warn("   [Warning] Found %i contigs totaling %i bp that are shared between multiple bins" % (numMultiplyBinnedContigs,
                                                                                                          multiplyBinnedContigsLength))
            self.logger.warn("   [Warning] Wrote summary to %s" % multiplyBinnedOutput.name)
