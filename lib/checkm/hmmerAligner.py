###############################################################################
#
# hmmerAlign.py - runs HMMER and provides functions for parsing output
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
import uuid
import logging
import tempfile
import multiprocessing as mp
from collections import defaultdict

from checkm.defaultValues import DefaultValues
from checkm.common import makeSurePathExists
from checkm.util.seqUtils import readFasta
from checkm.hmmer import HMMERRunner
from checkm.resultsParser import ResultsParser


class HmmerAligner:
    def __init__(self, threads):
        self.logger = logging.getLogger()
        self.totalThreads = threads

        self.outputFormat = 'Pfam'

    def makeAlignmentTopHit(self,
                               outDir,
                               hmmModelFile,
                               hmmTableFile,
                               binIdToModels,
                               bIgnoreThresholds,
                               evalueThreshold,
                               lengthThreshold,
                               bReportHitStats,
                               alignOutputDir,
                               bKeepUnmaskedAlign=False
                               ):
        """Align top hits in each bin. Assumes all bins are using the same marker genes."""

        self.logger.info("  Extracting marker genes to align.")

        # parse HMM information
        resultsParser = ResultsParser(binIdToModels)

        # get HMM hits to each bin
        resultsParser.parseBinHits(outDir, hmmTableFile, False, bIgnoreThresholds, evalueThreshold, lengthThreshold)

        # extract the ORFs to align
        markerSeqs, markerStats = self.__extractMarkerSeqsTopHits(outDir, resultsParser)

        # generate individual HMMs required to create multiple sequence alignments
        binId = binIdToModels.keys()[0]
        hmmModelFiles = {}
        self.__makeAlignmentModels(hmmModelFile, binIdToModels[binId], hmmModelFiles)

        # align each of the marker genes
        makeSurePathExists(alignOutputDir)
        self.__alignMarkerGenes(markerSeqs, markerStats, bReportHitStats, hmmModelFiles, alignOutputDir, bKeepUnmaskedAlign)

        # remove the temporary HMM files
        for fileName in hmmModelFiles:
            os.remove(hmmModelFiles[fileName])

        return resultsParser

    def makeAlignmentToPhyloMarkers(self,
                                       outDir,
                                       hmmModelFile,
                                       hmmTableFile,
                                       binIdToModels,
                                       bIgnoreThresholds,
                                       evalueThreshold,
                                       lengthThreshold,
                                       bReportHitStats,
                                       alignOutputDir,
                                       bKeepUnmaskedAlign=False
                                       ):
        """Align hits to a set of common marker genes."""

        self.logger.info("  Extracting marker genes to align.")

        # parse HMM information
        resultsParser = ResultsParser(binIdToModels)

        # get HMM hits to each bin
        resultsParser.parseBinHits(outDir, hmmTableFile, False, bIgnoreThresholds, evalueThreshold, lengthThreshold)

        # extract the ORFs to align
        markerSeqs, markerStats = self.__extractMarkerSeqsUnique(outDir, resultsParser)

        # generate individual HMMs required to create multiple sequence alignments
        binId = binIdToModels.keys()[0]
        hmmModelFiles = {}
        self.__makeAlignmentModels(hmmModelFile, binIdToModels[binId], hmmModelFiles)

        # align each of the marker genes
        makeSurePathExists(alignOutputDir)
        self.__alignMarkerGenes(markerSeqs, markerStats, bReportHitStats, hmmModelFiles, alignOutputDir, bKeepUnmaskedAlign)

        # remove the temporary HMM files
        for fileName in hmmModelFiles:
            os.remove(hmmModelFiles[fileName])

        return resultsParser

    def makeAlignmentsOfMultipleHits(self,
                                       outDir,
                                       markerFile,
                                       hmmTableFile,
                                       binIdToModels,
                                       binIdToBinMarkerSets,
                                       bIgnoreThresholds,
                                       evalueThreshold,
                                       lengthThreshold,
                                       alignOutputDir,
                                       ):
        """Align markers with multiple hits within a bin."""

        makeSurePathExists(alignOutputDir)

        # parse HMM information
        resultsParser = ResultsParser(binIdToModels)

        # get HMM hits to each bin
        resultsParser.parseBinHits(outDir, hmmTableFile, False, bIgnoreThresholds, evalueThreshold, lengthThreshold)

        # align any markers with multiple hits in a bin
        self.logger.info('  Aligning marker genes with multiple hits in a single bin:')

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binId in binIdToModels:
            workerQueue.put(binId)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        try:
            calcProc = [mp.Process(target=self.__createMSA, args=(resultsParser, binIdToBinMarkerSets, markerFile, outDir, alignOutputDir, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__reportBinProgress, args=(len(binIdToModels), writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

    def __createMSA(self, resultsParser, binIdToBinMarkerSets, hmmModelFile, outDir, alignOutputDir, queueIn, queueOut):
        """Create multiple sequence alignment for markers with multiple hits in a bin."""

        HF = HMMERRunner(mode='fetch')

        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            markersWithMultipleHits = self.__extractMarkersWithMultipleHits(outDir, binId, resultsParser, binIdToBinMarkerSets[binId])

            if len(markersWithMultipleHits) != 0:
                # create multiple sequence alignments for markers with multiple hits
                binAlignOutputDir = os.path.join(alignOutputDir, binId)
                makeSurePathExists(binAlignOutputDir)
                for markerId in markersWithMultipleHits:
                    tempModelFile = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
                    HF.fetch(hmmModelFile, markerId, tempModelFile)

                    self.__alignMarker(markerId, markersWithMultipleHits[markerId], None, False, binAlignOutputDir, tempModelFile, bKeepUnmaskedAlign=False)

                    os.remove(tempModelFile)

            queueOut.put(binId)

    def __reportBinProgress(self, numBins, queueIn):
        """Report number of processed bins."""

        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()

        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished processing %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def __alignMarkerGenes(self, markerSeqs, markerStats, bReportHitStats, hmmModelFiles, alignOutputDir, bKeepUnmaskedAlign=False, bReportProgress=True):
        """Align marker genes with HMMs in parallel."""

        if bReportProgress:
            self.logger.info("  Aligning %d marker genes with %d threads:" % (len(hmmModelFiles), self.totalThreads))

        # process each bin in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for markerId in hmmModelFiles:
            workerQueue.put(markerId)

        for _ in range(self.totalThreads):
            workerQueue.put(None)

        try:
            calcProc = [mp.Process(target=self.__alignMarkerParallel, args=(markerSeqs, markerStats, bReportHitStats, alignOutputDir, hmmModelFiles, bKeepUnmaskedAlign, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__reportAlignmentProgress, args=(len(hmmModelFiles), bReportProgress, writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

    def __alignMarkerParallel(self, markerSeqs, markerStats, bReportHitStats, alignOutputDir, hmmModelFiles, bKeepUnmaskedAlign, queueIn, queueOut):
        while True:
            markerId = queueIn.get(block=True, timeout=None)
            if markerId == None:
                break

            self.__alignMarker(markerId, markerSeqs[markerId], markerStats[markerId], bReportHitStats, alignOutputDir, hmmModelFiles[markerId], bKeepUnmaskedAlign)

            queueOut.put(markerId)

    def __alignMarker(self, markerId, binSeqs, binStats, bReportHitStats, alignOutputDir, hmmModelFile, bKeepUnmaskedAlign):
        unalignSeqFile = os.path.join(alignOutputDir, markerId + '.unaligned.faa')
        fout = open(unalignSeqFile, 'w')
        numSeqs = 0
        for binId, seqs in binSeqs.iteritems():
            for seqId, seq in seqs.iteritems():
                header = '>' + binId + DefaultValues.SEQ_CONCAT_CHAR + seqId
                if bReportHitStats:
                    header += ' [e-value=%.4g,score=%.1f]' % (binStats[binId][seqId][0], binStats[binId][seqId][1])

                fout.write(header + '\n')
                fout.write(seq + '\n')
                numSeqs += 1
        fout.close()

        if numSeqs > 0:
            alignSeqFile = os.path.join(alignOutputDir, markerId + '.aligned.faa')
            HA = HMMERRunner(mode='align')
            HA.align(hmmModelFile, unalignSeqFile, alignSeqFile, writeMode='>', outputFormat=self.outputFormat, trim=False)

            makedSeqFile = os.path.join(alignOutputDir, markerId + '.masked.faa')
            self.__maskAlignment(alignSeqFile, makedSeqFile)

            if not bKeepUnmaskedAlign:
                os.remove(alignSeqFile)

        os.remove(unalignSeqFile)

    def __reportAlignmentProgress(self, numMarkers, bReportProgress, queueIn):
        """Report number of processed markers."""

        numProcessedGenes = 0
        if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished aligning %d of %d (%.2f%%) marker genes.' % (numProcessedGenes, numMarkers, float(numProcessedGenes) * 100 / numMarkers)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()

        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedGenes += 1
                statusStr = '    Finished aligning %d of %d (%.2f%%) marker genes.' % (numProcessedGenes, numMarkers, float(numProcessedGenes) * 100 / numMarkers)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def __maskAlignment(self, inputFile, outputFile):
        """Read HMMER alignment in STOCKHOLM format and output masked alignment in FASTA format."""
        # read STOCKHOLM alignment
        seqs = {}
        seqStats = {}
        for line in open(inputFile):
            line = line.rstrip()
            if line == '' or line[0] == '#' or line == '//':
                if 'GC RF' in line:
                    mask = line.split('GC RF')[1].strip()
                elif '=GS' in line:
                    # read additional sequence informations
                    lineSplit = line.split()
                    seqId = lineSplit[1]
                    stats = lineSplit[3].strip()
                    seqStats[seqId] = stats
                continue
            else:
                lineSplit = line.split()
                seqs[lineSplit[0]] = lineSplit[1].upper().replace('.', '-').strip()

        # output masked sequences in FASTA format
        fout = open(outputFile, 'w')
        for seqId, seq in seqs.iteritems():
            if seqStats:
                fout.write('>%s %s\n' % (seqId, seqStats[seqId]))
            else:
                fout.write('>' + seqId + '\n')

            maskedSeq = ''.join([seq[i] for i in xrange(0, len(seq)) if mask[i] == 'x'])
            fout.write(maskedSeq + '\n')
        fout.close()

    def __extractMarkerSeqsTopHits(self, outDir, resultsParser):
        """Extract marker sequences from top hits (assume all bins use the same HMM file)."""

        markerSeqs = defaultdict(dict)
        markerStats = defaultdict(dict)
        for binId in resultsParser.results:
            # read ORFs for bin
            aaGeneFile = os.path.join(outDir, 'bins', binId, DefaultValues.PRODIGAL_AA)
            binORFs = readFasta(aaGeneFile)

            # extract ORFs hitting a marker
            for markerId, hits in resultsParser.results[binId].markerHits.iteritems():
                markerSeqs[markerId][binId] = {}
                markerStats[markerId][binId] = {}

                # sort hits from highest to lowest e-value in order to ensure only the best hit
                # to a given target is retained
                hits.sort(key=lambda x: x.full_e_value, reverse=True)
                topHit = hits[0]
                markerSeqs[markerId][binId][topHit.target_name] = self.__extractSeq(topHit.target_name, binORFs)
                markerStats[markerId][binId][topHit.target_name] = [topHit.full_e_value, topHit.full_score]

        return markerSeqs, markerStats

    def __extractMarkerSeqsUnique(self, outDir, resultsParser):
        """Extract marker sequences with a single unique hit."""

        markerSeqs = defaultdict(dict)
        markerStats = defaultdict(dict)
        for binId in resultsParser.results:
            # read ORFs for bin
            aaGeneFile = os.path.join(outDir, 'bins', binId, DefaultValues.PRODIGAL_AA)
            binORFs = readFasta(aaGeneFile)

            # extract ORFs hitting a marker
            for markerId, hits in resultsParser.results[binId].markerHits.iteritems():
                markerSeqs[markerId][binId] = {}
                markerStats[markerId][binId] = {}

                # only record hits which are unique
                if len(hits) == 1:
                    hit = hits[0]
                    markerSeqs[markerId][binId][hit.target_name] = self.__extractSeq(hit.target_name, binORFs)
                    markerStats[markerId][binId][hit.target_name] = [hit.full_e_value, hit.full_score]

        return markerSeqs, markerStats

    def __extractSeq(self, seqId, seqs):
        """Extract sequence data."""

        if DefaultValues.SEQ_CONCAT_CHAR in seqId:
            seqIds = seqId.split(DefaultValues.SEQ_CONCAT_CHAR)

            seq = ''
            for seqId in seqIds:
                tempSeq = seqs[seqId]
                if tempSeq[-1] == '*':
                    tempSeq = tempSeq[0:-1]  # remove final '*' inserted by prodigal

                seq += tempSeq

            rtnSeq = seq
        else:
            rtnSeq = seqs[seqId]

            if rtnSeq[-1] == '*':
                rtnSeq = rtnSeq[0:-1]  # remove final '*' inserted by prodigal

        return rtnSeq

    def __extractMarkersWithMultipleHits(self, outDir, binId, resultsParser, binMarkerSet):
        """Extract markers with multiple hits within a single bin."""

        markersWithMultipleHits = defaultdict(dict)

        aaGeneFile = os.path.join(outDir, 'bins', binId, DefaultValues.PRODIGAL_AA)
        binORFs = readFasta(aaGeneFile)

        markerGenes = binMarkerSet.selectedMarkerSet().getMarkerGenes()
        for markerId, hits in resultsParser.results[binId].markerHits.iteritems():
            if markerId not in markerGenes or len(hits) < 2:
                continue

            # sort hits from highest to lowest e-value in order to ensure only the best hit
            # to a given target is retained
            hits.sort(key=lambda x: x.full_e_value, reverse=True)

            # Note: this data structure is used to mimic that used by __extractMarkerSeqsTopHits()
            markersWithMultipleHits[markerId][binId] = {}
            for hit in hits:
                markersWithMultipleHits[markerId][binId][hit.target_name] = self.__extractSeq(hit.target_name, binORFs)

        return markersWithMultipleHits

    def __makeAlignmentModels(self, hmmModelFile, modelKeys, hmmModelFiles, bReportProgress=True):
        """Make temporary HMM files used to create HMM alignments."""

        if bReportProgress:
            self.logger.info("  Extracting %d HMMs with %d threads:" % (len(modelKeys), self.totalThreads))

        # process each marker in parallel
        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for modelId in modelKeys:
            fetchFilename = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
            hmmModelFiles[modelId] = fetchFilename
            workerQueue.put((modelId, fetchFilename))

        for _ in range(self.totalThreads):
            workerQueue.put((None, None))

        try:
            calcProc = [mp.Process(target=self.__extractModel, args=(hmmModelFile, workerQueue, writerQueue)) for _ in range(self.totalThreads)]
            writeProc = mp.Process(target=self.__reportModelExtractionProgress, args=(len(modelKeys), bReportProgress, writerQueue))

            writeProc.start()

            for p in calcProc:
                p.start()

            for p in calcProc:
                p.join()

            writerQueue.put(None)
            writeProc.join()
        except:
            # make sure all processes are terminated
            for p in calcProc:
                p.terminate()

            writeProc.terminate()

    def __extractModel(self, hmmModelFile, queueIn, queueOut):
        """Extract HMM."""
        HF = HMMERRunner(mode='fetch')

        while True:
            modelId, fetchFilename = queueIn.get(block=True, timeout=None)
            if modelId == None:
                break

            HF.fetch(hmmModelFile, modelId, fetchFilename)

            queueOut.put(modelId)

    def __reportModelExtractionProgress(self, numMarkers, bReportProgress, queueIn):
        """Report number of extracted HMMs."""

        numModelsExtracted = 0
        if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
            statusStr = '    Finished extracting %d of %d (%.2f%%) HMMs.' % (numModelsExtracted, numMarkers, float(numModelsExtracted) * 100 / numMarkers)
            sys.stderr.write('%s\r' % statusStr)
            sys.stderr.flush()

        while True:
            modelId = queueIn.get(block=True, timeout=None)
            if modelId == None:
                break

            if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
                numModelsExtracted += 1
                statusStr = '    Finished extracting %d of %d (%.2f%%) HMMs.' % (numModelsExtracted, numMarkers, float(numModelsExtracted) * 100 / numMarkers)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if bReportProgress and self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')
