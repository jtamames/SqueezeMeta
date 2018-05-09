###############################################################################
#
# binStatistics.py - calculate statistics for each putative genome bin
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
import multiprocessing as mp
import math
import logging

from checkm.defaultValues import DefaultValues
from checkm.util.seqUtils import readFasta, baseCount, calculateN50
from checkm.common import binIdFromFilename, makeSurePathExists
from checkm.prodigal import ProdigalGeneFeatureParser

from numpy import mean


class BinStatistics():
    """Statistics for genomes.

    This class calculates statistics for genomes comprised
    of one or more scaffolds. Mean statistics are weighted by
    scaffold length. The following statistics are calculated:
     - bin assignment
     - mean GC
     - mean and median scaffold length
     - N50
     - mean coverage
     - mean tetranucleotide signature

    """

    def __init__(self, threads=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger()
        self.totalThreads = threads

    def calculate(self, binFiles, outDir, binStatsFile):
        """Calculate statistics for each putative genome bin."""

        # process each bin
        self.logger.info("  Calculating genome statistics for %d bins with %d threads:" % (len(binFiles), self.totalThreads))

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binFile in binFiles:
            workerQueue.put(binFile)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        try:
            calcProc = [mp.Process(target=self.__processBin, args=(outDir, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__reportProgress, args=(outDir, binStatsFile, len(binFiles), writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put((None, None))
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

    def __processBin(self, outDir, queueIn, queueOut):
        """Thread safe bin processing."""
        while True:
            binFile = queueIn.get(block=True, timeout=None)
            if binFile == None:
                break

            binStats = {}

            binId = binIdFromFilename(binFile)
            binDir = os.path.join(outDir, 'bins', binId)
            makeSurePathExists(binDir)

            # read scaffolds
            scaffolds = readFasta(binFile)

            # calculate GC statistics
            GC, stdGC = self.calculateGC(scaffolds)
            binStats['GC'] = GC
            binStats['GC std'] = stdGC

            # calculate statistics related to contigs and scaffolds
            maxScaffoldLen, maxContigLen, genomeSize, scaffold_N50, contig_N50, scaffoldAvgLen, contigAvgLen, numContigs, numAmbiguousBases = self.calculateSeqStats(scaffolds)
            binStats['Genome size'] = genomeSize
            binStats['# ambiguous bases'] = numAmbiguousBases
            binStats['# scaffolds'] = len(scaffolds)
            binStats['# contigs'] = numContigs
            binStats['Longest scaffold'] = maxScaffoldLen
            binStats['Longest contig'] = maxContigLen
            binStats['N50 (scaffolds)'] = scaffold_N50
            binStats['N50 (contigs)'] = contig_N50
            binStats['Mean scaffold length'] = scaffoldAvgLen
            binStats['Mean contig length'] = contigAvgLen

            # calculate coding density statistics
            codingDensity, translationTable, numORFs = self.calculateCodingDensity(binDir, scaffolds, genomeSize)
            binStats['Coding density'] = codingDensity
            binStats['Translation table'] = translationTable
            binStats['# predicted genes'] = numORFs

            queueOut.put((binId, binStats))

    def __reportProgress(self, outDir, binStatsFile, numBins, queueIn):
        """Report number of processed bins and write statistics to file."""

        storagePath = os.path.join(outDir, 'storage')
        fout = open(os.path.join(storagePath, binStatsFile), 'w')

        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()

        while True:
            binId, curBinStats = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            fout.write(binId + '\t' + str(curBinStats) + '\n')
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

        fout.close()

    def calculateGC(self, seqs, seqStats=None):
        """Calculate fraction of nucleotides that are G or C."""
        totalGC = 0
        totalAT = 0
        gcPerSeq = []
        for seqId, seq in seqs.iteritems():
            a, c, g, t = baseCount(seq)

            gc = g + c
            at = a + t

            totalGC += gc
            totalAT += at

            if (gc + at) > 0:
                gcContent = float(gc) / (gc + at)
            else:
                gcContent = 0.0

            if seqStats:
                seqStats[seqId]['GC'] = gcContent

            if len(seq) > DefaultValues.MIN_SEQ_LEN_GC_STD:
                gcPerSeq.append(gcContent)

        if (totalGC + totalAT) > 0:
            GC = float(totalGC) / (totalGC + totalAT)
        else:
            GC = 0.0

        varGC = 0
        if len(gcPerSeq) > 1:
            varGC = mean(map(lambda x: (x - GC) ** 2, gcPerSeq))

        return GC, math.sqrt(varGC)

    def calculateSeqStats(self, scaffolds, seqStats=None):
        """Calculate scaffold length statistics (min length, max length, total length, N50, # contigs)."""
        scaffoldLens = []
        contigLens = []
        numAmbiguousBases = 0
        for scaffoldId, scaffold in scaffolds.iteritems():
            scaffoldLen = len(scaffold)
            scaffoldLens.append(scaffoldLen)

            splitScaffold = scaffold.split(DefaultValues.CONTIG_BREAK)
            lenContigsInScaffold = []
            for contig in splitScaffold:
                contigLen = len(contig.replace('N', ''))
                if contigLen > 0:
                    lenContigsInScaffold.append(contigLen)

            contigLens += lenContigsInScaffold

            if seqStats:
                seqStats[scaffoldId]['Length'] = scaffoldLen
                seqStats[scaffoldId]['Total contig length'] = sum(lenContigsInScaffold)
                seqStats[scaffoldId]['# contigs'] = len(lenContigsInScaffold)

            numAmbiguousBases += scaffold.count('N') + scaffold.count('n')

        scaffold_N50 = calculateN50(scaffoldLens)
        contig_N50 = calculateN50(contigLens)

        return max(scaffoldLens), max(contigLens), sum(scaffoldLens), scaffold_N50, contig_N50, mean(scaffoldLens), mean(contigLens), len(contigLens), numAmbiguousBases

    def calculateCodingDensity(self, outDir, scaffolds, genomeSize):
        """Calculate coding density of putative genome bin."""
        gffFile = os.path.join(outDir, DefaultValues.PRODIGAL_GFF)
        if os.path.exists(gffFile):
            prodigalParserGFF = ProdigalGeneFeatureParser(gffFile)

            aaFile = os.path.join(outDir, DefaultValues.PRODIGAL_AA)  # use AA file as nucleotide file is optional
            aaGenes = readFasta(aaFile)

            codingBasePairs = 0  # self.__calculateCodingBases(aaGenes)
            for scaffold_id in scaffolds.keys():
                codingBasePairs += prodigalParserGFF.codingBases(scaffold_id)

            return float(codingBasePairs) / genomeSize, prodigalParserGFF.translationTable, len(aaGenes)
        else:
            # there is no gene feature file (perhaps the user specified pre-calculated genes)
            # so calculating the coding density is not possible
            return -1, -1, -1

    def __calculateCodingBases(self, aaGenes):
        """Calculate number of coding bases in a set of genes."""
        codingBasePairs = 0
        for _geneId, gene in aaGenes.iteritems():
            codingBasePairs += len(gene) * 3

        return codingBasePairs

    def sequenceStats(self, outDir, binFile):
        """Calculate statistics for all sequences within a bin."""

        # read scaffolds
        seqs = readFasta(binFile)

        seqStats = {}
        for seqId in seqs:
            seqStats[seqId] = {}

        self.calculateGC(seqs, seqStats)
        self.calculateSeqStats(seqs, seqStats)

        binId = binIdFromFilename(binFile)
        aaFile = os.path.join(outDir, 'bins', binId, DefaultValues.PRODIGAL_AA)
        if os.path.exists(aaFile):
            aaGenes = readFasta(aaFile)
            for geneId, gene in aaGenes.iteritems():
                seqId = geneId[0:geneId.rfind('_')]
                seqStats[seqId]['# ORFs'] = seqStats[seqId].get('# ORFs', 0) + 1
                seqStats[seqId]['Coding bases'] = seqStats[seqId].get('Coding bases', 0) + len(gene) * 3
        else:
            # missing amino acid file likely indicates users used a pre-called gene file, so
            # just set some defaults
            seqStats[seqId]['# ORFs'] = seqStats[seqId].get('# ORFs', 0) + 1
            seqStats[seqId]['Coding bases'] = seqStats[seqId].get('Coding bases', 0) + len(gene) * 3

        return seqStats
