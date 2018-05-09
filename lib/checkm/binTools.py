###############################################################################
#
# binTools.py - functions for exploring and modifying bins
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
import gzip

import numpy as np

from common import binIdFromFilename, checkFileExists, readDistribution, findNearest
from checkm.util.seqUtils import readFasta, writeFasta, baseCount
from checkm.genomicSignatures import GenomicSignatures
from checkm.prodigal import ProdigalGeneFeatureParser
from checkm.defaultValues import DefaultValues


class BinTools():
    """Functions for exploring and modifying bins."""
    def __init__(self, threads=1):
        self.logger = logging.getLogger()

    def __removeSeqs(self, seqs, seqsToRemove):
        """Remove sequences. """
        missingSeqIds = set(seqsToRemove).difference(set(seqs.keys()))
        if len(missingSeqIds) > 0:
            self.logger.error('  [Error] Missing sequence(s) specified for removal: ' + ', '.join(missingSeqIds) + '\n')
            sys.exit()

        for seqId in seqsToRemove:
            seqs.pop(seqId)

    def __addSeqs(self, seqs, refSeqs, seqsToAdd):
        """Add sequences. """
        missingSeqIds = set(seqsToAdd).difference(set(refSeqs.keys()))
        if len(missingSeqIds) > 0:
            self.logger.error('  [Error] Missing sequence(s) specified for addition: ' + ', '.join(missingSeqIds) + '\n')
            sys.exit()

        for seqId in seqsToAdd:
            seqs[seqId] = refSeqs[seqId]

    def modify(self, binFile, seqFile, seqsToAdd, seqsToRemove, outputFile):
        """Add and remove sequences from a file."""
        binSeqs = readFasta(binFile)

        # add sequences to bin
        if seqsToAdd != None:
            refSeqs = readFasta(seqFile)
            self.__addSeqs(binSeqs, refSeqs, seqsToAdd)

        # remove sequences from bin
        if seqsToRemove != None:
            self.__removeSeqs(binSeqs, seqsToRemove)

        # save modified bin
        writeFasta(binSeqs, outputFile)

    def removeOutliers(self, binFile, outlierFile, outputFile):
        """Remove sequences specified as outliers in the provided file."""

        binSeqs = readFasta(binFile)
        binIdToModify = binIdFromFilename(binFile)

        # get files to remove
        checkFileExists(outlierFile)
        seqsToRemove = []
        bHeader = True
        for line in open(outlierFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            binId = lineSplit[0]

            if binId == binIdToModify:
                seqId = lineSplit[1]
                seqsToRemove.append(seqId)

        # remove sequences from bin
        if len(seqsToRemove) > 0:
            self.__removeSeqs(binSeqs, seqsToRemove)

        # save modified bin
        writeFasta(binSeqs, outputFile)

    def unique(self, binFiles):
        """Check if sequences are assigned to multiple bins."""

        # read sequence IDs from all bins,
        # while checking for duplicate sequences within a bin
        binSeqs = {}
        for f in binFiles:
            binId = binIdFromFilename(f)

            if f.endswith('.gz'):
                openFile = gzip.open
            else:
                openFile = open

            seqIds = set()
            for line in openFile(f):
                if line[0] == '>':
                    seqId = line[1:].split(None, 1)[0]

                    if seqId in seqIds:
                        print '  [Warning] Sequence %s found multiple times in bin %s.' % (seqId, binId)
                    seqIds.add(seqId)

            binSeqs[binId] = seqIds

        # check for sequences assigned to multiple bins
        bDuplicates = False
        binIds = binSeqs.keys()
        for i in xrange(0, len(binIds)):
            for j in xrange(i + 1, len(binIds)):
                seqInter = set(binSeqs[binIds[i]]).intersection(set(binSeqs[binIds[j]]))

                if len(seqInter) > 0:
                    bDuplicates = True
                    print '  Sequences shared between %s and %s: ' % (binIds[i], binIds[j])
                    for seqId in seqInter:
                        print '    ' + seqId
                    print ''

        if not bDuplicates:
            print '  No sequences assigned to multiple bins.'

    def gcDist(self, seqs):
        """GC statistics for bin."""
        GCs = []
        gcTotal = 0
        basesTotal = 0
        for _, seq in seqs.iteritems():
            a, c, g, t = baseCount(seq)
            gc = g + c
            bases = a + c + g + t

            GCs.append(float(gc) / (bases))

            gcTotal += gc
            basesTotal += bases

        meanGC = float(gcTotal) / basesTotal
        deltaGCs = np.array(GCs) - meanGC

        return meanGC, deltaGCs, GCs

    def codingDensityDist(self, seqs, prodigalParser):
        """Coding density statistics for bin."""
        CDs = []

        codingBasesTotal = 0
        basesTotal = 0
        for seqId, seq in seqs.iteritems():
            codingBases = prodigalParser.codingBases(seqId)

            CDs.append(float(codingBases) / len(seq))
            codingBasesTotal += codingBases
            basesTotal += len(seq)

        meanCD = float(codingBasesTotal) / basesTotal
        deltaCDs = np.array(CDs) - meanCD

        return meanCD, deltaCDs, CDs

    def binTetraSig(self, seqs, tetraSigs):
        """Tetranucleotide signature for bin. """
        binSize = 0
        for _, seq in seqs.iteritems():
            binSize += len(seq)

        bInit = True
        for seqId, seq in seqs.iteritems():
            weightedTetraSig = tetraSigs[seqId] * (float(len(seq)) / binSize)
            if bInit:
                binSig = weightedTetraSig
                bInit = False
            else:
                binSig += weightedTetraSig

        return binSig

    def tetraDiffDist(self, seqs, genomicSig, tetraSigs, binSig):
        """TD statistics for bin."""
        deltaTDs = np.zeros(len(seqs))
        for i, seqId in enumerate(seqs.keys()):
            deltaTDs[i] = genomicSig.distance(tetraSigs[seqId], binSig)

        return np.mean(deltaTDs), deltaTDs

    def identifyOutliers(self, outDir, binFiles, tetraProfileFile, distribution, reportType, outputFile):
        """Identify sequences that are outliers."""

        self.logger.info('  Reading reference distributions.')
        gcBounds = readDistribution('gc_dist')
        cdBounds = readDistribution('cd_dist')
        tdBounds = readDistribution('td_dist')

        fout = open(outputFile, 'w')
        fout.write('Bin Id\tSequence Id\tSequence length\tOutlying distributions')
        fout.write('\tSequence GC\tMean bin GC\tLower GC bound (%s%%)\tUpper GC bound (%s%%)' % (distribution, distribution))
        fout.write('\tSequence CD\tMean bin CD\tLower CD bound (%s%%)' % distribution)
        fout.write('\tSequence TD\tMean bin TD\tUpper TD bound (%s%%)\n' % distribution)

        self.logger.info('')
        processedBins = 0
        for binFile in binFiles:
            binId = binIdFromFilename(binFile)

            processedBins += 1
            self.logger.info('  Finding outliers in %s (%d of %d).' % (binId, processedBins, len(binFiles)))

            seqs = readFasta(binFile)

            meanGC, deltaGCs, seqGC = self.gcDist(seqs)

            genomicSig = GenomicSignatures(K=4, threads=1)
            tetraSigs = genomicSig.read(tetraProfileFile)
            binSig = self.binTetraSig(seqs, tetraSigs)
            meanTD, deltaTDs = self.tetraDiffDist(seqs, genomicSig, tetraSigs, binSig)

            gffFile = os.path.join(outDir, 'bins', binId, DefaultValues.PRODIGAL_GFF)
            if not os.path.exists(gffFile):
                self.logger.error('  [Error] Missing gene feature file (%s). This plot if not compatible with the --genes option.\n' % DefaultValues.PRODIGAL_GFF)
                sys.exit()

            prodigalParser = ProdigalGeneFeatureParser(gffFile)
            meanCD, deltaCDs, CDs = self.codingDensityDist(seqs, prodigalParser)

            # find keys into GC and CD distributions
            closestGC = findNearest(np.array(gcBounds.keys()), meanGC)
            sampleSeqLen = gcBounds[closestGC].keys()[0]
            d = gcBounds[closestGC][sampleSeqLen]
            gcLowerBoundKey = findNearest(d.keys(), (100 - distribution) / 2.0)
            gcUpperBoundKey = findNearest(d.keys(), (100 + distribution) / 2.0)

            closestCD = findNearest(np.array(cdBounds.keys()), meanCD)
            sampleSeqLen = cdBounds[closestCD].keys()[0]
            d = cdBounds[closestCD][sampleSeqLen]
            cdLowerBoundKey = findNearest(d.keys(), (100 - distribution) / 2.0)

            tdBoundKey = findNearest(tdBounds[tdBounds.keys()[0]].keys(), distribution)

            index = 0
            for seqId, seq in seqs.iteritems():
                seqLen = len(seq)

                # find GC, CD, and TD bounds
                closestSeqLen = findNearest(gcBounds[closestGC].keys(), seqLen)
                gcLowerBound = gcBounds[closestGC][closestSeqLen][gcLowerBoundKey]
                gcUpperBound = gcBounds[closestGC][closestSeqLen][gcUpperBoundKey]

                closestSeqLen = findNearest(cdBounds[closestCD].keys(), seqLen)
                cdLowerBound = cdBounds[closestCD][closestSeqLen][cdLowerBoundKey]

                closestSeqLen = findNearest(tdBounds.keys(), seqLen)
                tdBound = tdBounds[closestSeqLen][tdBoundKey]

                outlyingDists = []
                if deltaGCs[index] < gcLowerBound or deltaGCs[index] > gcUpperBound:
                    outlyingDists.append('GC')

                if deltaCDs[index] < cdLowerBound:
                    outlyingDists.append('CD')

                if deltaTDs[index] > tdBound:
                    outlyingDists.append('TD')

                if (reportType == 'any' and len(outlyingDists) >= 1) or (reportType == 'all' and len(outlyingDists) == 3):
                    fout.write(binId + '\t' + seqId + '\t%d' % len(seq) + '\t' + ','.join(outlyingDists))
                    fout.write('\t%.1f\t%.1f\t%.1f\t%.1f' % (seqGC[index] * 100, meanGC * 100, (meanGC + gcLowerBound) * 100, (meanGC + gcUpperBound) * 100))
                    fout.write('\t%.1f\t%.1f\t%.1f' % (CDs[index] * 100, meanCD * 100, (meanCD + cdLowerBound) * 100))
                    fout.write('\t%.3f\t%.3f\t%.3f' % (deltaTDs[index], meanTD, tdBound) + '\n')

                index += 1

        fout.close()
