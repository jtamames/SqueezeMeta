###############################################################################
#
# markerSet.py - Calculate and process marker sets.
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
import uuid
import tempfile
import shutil
import multiprocessing as mp
import cPickle as pickle
import gzip

from checkm.defaultValues import DefaultValues
from checkm.hmmer import HMMERRunner
from checkm.hmmerModelParser import HmmModelParser

from checkm.util.pfam import PFAM


class BinMarkerSets():
    """A collection of one or more marker sets associated with a bin."""

    # type of marker set
    TAXONOMIC_MARKER_SET = 1
    TREE_MARKER_SET = 2
    HMM_MODELS_SET = 3

    def __init__(self, binId, markerSetType):
        self.logger = logging.getLogger()
        self.markerSets = []
        self.binId = binId
        self.markerSetType = markerSetType
        self.selectedLinageSpecificMarkerSet = None

    def numMarkerSets(self):
        """Number of marker sets associated with bin."""
        return len(self.markerSets)

    def addMarkerSet(self, markerSet):
        """Add marker set to bin."""
        self.markerSets.append(markerSet)

    def markerSetIter(self):
        """Generator function for iterating over marker sets."""
        for markerSet in self.markerSets:
            yield markerSet

    def getMarkerGenes(self):
        """Get marker genes from all marker sets."""
        markerGenes = set()
        for ms in self.markerSets:
            markerGenes.update(ms.getMarkerGenes())

        return markerGenes

    def mostSpecificMarkerSet(self):
        return self.markerSets[0]

    def treeMarkerSet(self):
        pass

    def selectedMarkerSet(self):
        """Return the 'selected' marker set for this bin."""
        if self.markerSetType == self.TAXONOMIC_MARKER_SET:
            return self.mostSpecificMarkerSet()
        elif self.markerSetType == self.TREE_MARKER_SET:
            return self.selectedLinageSpecificMarkerSet
        else:
            # there should be a single marker set associate with this bin
            if len(self.markerSets) == 1:
                return self.markerSets[0]

            self.logger.error('  [Error] Expect a single marker set to be associated with each bin.')
            sys.exit()

    def setLineageSpecificSelectedMarkerSet(self, selectedMarkerSetMap):
        uid = self.mostSpecificMarkerSet().UID

        selectedId = selectedMarkerSetMap[uid]
        self.selectedLinageSpecificMarkerSet = None
        
        while not self.selectedLinageSpecificMarkerSet:
            for ms in self.markerSets:
                if ms.UID == selectedId:
                    self.selectedLinageSpecificMarkerSet = ms
                    break
                
            if not self.selectedLinageSpecificMarkerSet:
                # This is a hack for the reduced tree. Since not all
                # marker sets are in the reduced tree it is possible the
                # selected marker set might not be avaliable. In this case,
                # we should move to the next suitable marker set. Ideally,
                # this could be avoided by just forcing in the selected
                # marker set.
                selectedId = selectedMarkerSetMap[selectedId]
            else:
                break

        if self.selectedLinageSpecificMarkerSet == None:
            # something has gone wrong
            self.logger.error('  [Error] Failed to set a selected lineage-specific marker set.')
            sys.exit()

    def removeMarkers(self, markersToRemove):
        """Remove specified markers from all marker sets."""
        for markerSet in self.markerSets:
            curMarkerSet = markerSet.markerSet

            newMarkerSet = []
            for ms in curMarkerSet:
                newMS = ms - markersToRemove

                if len(newMS) != 0:
                    newMarkerSet.append(newMS)

            markerSet.markerSet = newMarkerSet

    def write(self, fout):
        """Write marker set to file."""
        fout.write(self.binId)
        fout.write('\t' + str(len(self.markerSets)))
        for ms in self.markerSets:
            fout.write('\t' + str(ms))
        fout.write('\n')

    def read(self, line):
        """Construct bin marker set data from line."""
        lineSplit = line.split('\t')
        numMarkerSets = int(lineSplit[1])
        for i in xrange(0, numMarkerSets):
            uid = lineSplit[i * 4 + 2]
            lineageStr = lineSplit[i * 4 + 3]
            numGenomes = int(lineSplit[i * 4 + 4])
            markerSet = eval(lineSplit[i * 4 + 5])
            self.markerSets.append(MarkerSet(uid, lineageStr, numGenomes, markerSet))


class MarkerSet():
    """A collection of marker genes organized into co-located sets."""
    def __init__(self, UID, lineageStr, numGenomes, markerSet):
        self.logger = logging.getLogger()

        self.UID = UID  # unique ID of marker set
        self.lineageStr = lineageStr  # taxonomic string associated with marker set
        self.numGenomes = numGenomes  # number of genomes used to calculate marker set
        self.markerSet = markerSet  # marker genes organized into co-located sets

    def __repr__(self):
        return str(self.UID) + '\t' + self.lineageStr + '\t' + str(self.numGenomes) + '\t' + str(self.markerSet)

    def size(self):
        """Number of marker genes and marker gene sets."""
        numMarkerGenes = 0
        for m in self.markerSet:
            numMarkerGenes += len(m)

        return numMarkerGenes, len(self.markerSet)

    def numMarkers(self):
        """Number of marker genes."""
        return self.size()[0]

    def numSets(self):
        """Number of marker sets."""
        return len(self.markerSet)

    def getMarkerGenes(self):
        """Get marker genes within marker set."""
        markerGenes = set()
        for m in self.markerSet:
            for marker in m:
                markerGenes.add(marker)

        return markerGenes
        
    def removeMarkers(self, markersToRemove):
        """Remove specified markers from marker sets."""
        newMarkerSet = []
        for ms in self.markerSet:
            newMS = ms - markersToRemove

            if len(newMS) != 0:
                newMarkerSet.append(newMS)

        self.markerSet = newMarkerSet

    def genomeCheck(self, hits, bIndividualMarkers):
        """Calculate genome completeness and contamination."""
        if bIndividualMarkers:
            present = 0
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                if marker in hits:
                    present += 1
                    multiCopyCount += (len(hits[marker]) - 1)

            percComp = 100 * float(present) / self.numMarkers()
            percCont = 100 * float(multiCopyCount) / self.numMarkers()
        else:
            comp = 0.0
            cont = 0.0
            for ms in self.markerSet:
                present = 0
                multiCopy = 0
                for marker in ms:
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multiCopy += (count - 1)

                comp += float(present) / len(ms)
                cont += float(multiCopy) / len(ms)

            percComp = 100 * comp / len(self.markerSet)
            percCont = 100 * cont / len(self.markerSet)

        return percComp, percCont


class MarkerSetParser():
    """Parse marker set file."""

    def __init__(self, threads=1):
        self.logger = logging.getLogger()
        self.numThreads = threads

    def getMarkerSets(self, outDir, binIds, markerFile, excludeMarkersFile=None):
        """Determine marker set for each bin."""

        # determine type of marker set file
        markerFileType = self.markerFileType(markerFile)

        # get marker set for each bin
        binIdToBinMarkerSets = {}

        if markerFileType == BinMarkerSets.TAXONOMIC_MARKER_SET:
            binMarkerSets = self.parseTaxonomicMarkerSetFile(markerFile)

            for binId in binIds:
                binIdToBinMarkerSets[binId] = binMarkerSets
        elif markerFileType == BinMarkerSets.TREE_MARKER_SET:
            binIdToBinMarkerSets = self.parseLineageMarkerSetFile(markerFile)
        else:
            markers = [set()]
            modelParser = HmmModelParser(markerFile)
            for model in modelParser.parse():
                markers[0].add(model.acc)
            markerSet = MarkerSet(0, "N/A", -1, markers)

            for binId in binIds:
                binMarkerSets = BinMarkerSets(binId, BinMarkerSets.HMM_MODELS_SET)
                binMarkerSets.addMarkerSet(markerSet)
                binIdToBinMarkerSets[binId] = binMarkerSets

        # remove marker genes specified by user or marker for exclusion
        markersToExclude = set()
        if excludeMarkersFile:
            markersToExclude = self.readExcludeMarkersFile(excludeMarkersFile)

        markersToExclude.update(DefaultValues.MARKERS_TO_EXCLUDE)
        for binId, binMarkerSet in binIdToBinMarkerSets.iteritems():
            binMarkerSet.removeMarkers(markersToExclude)

        return binIdToBinMarkerSets

    def readExcludeMarkersFile(self, excludeMarkersFile):
        """Parse file specifying markers to exclude."""
        markersToExclude = set()
        for line in open(excludeMarkersFile):
            if line[0] == '#':
                continue

            marker = line.strip()
            markersToExclude.add(marker)

        return markersToExclude

    def createHmmModels(self, outDir, binIds, markerFile):
        """Create HMM model for each bins marker set."""

        # determine type of marker set file
        markerFileType = self.markerFileType(markerFile)

        # get HMM file for each bin
        binIdToModels = {}
        if markerFileType == BinMarkerSets.TAXONOMIC_MARKER_SET:
            hmmModelFile = self.createHmmModelFile(binIds.keys()[0], markerFile)

            modelParser = HmmModelParser(hmmModelFile)
            models = modelParser.models()
            for binId in binIds:
                binIdToModels[binId] = models

            os.remove(hmmModelFile)
        elif markerFileType == BinMarkerSets.TREE_MARKER_SET:
            binIdToModels = self.__createLineageHmmModels(binIds, markerFile)
        else:
            modelParser = HmmModelParser(markerFile)
            models = modelParser.models()
            for binId in binIds:
                binIdToModels[binId] = models

        return binIdToModels

    def createHmmModelFile(self, binId, markerFile):
        """Create HMM file for from a bin's marker set."""

        # determine type of marker set file
        markerFileType = self.markerFileType(markerFile)

        # create HMM file
        hmmModelFile = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        if markerFileType == BinMarkerSets.TAXONOMIC_MARKER_SET:
            binMarkerSets = self.parseTaxonomicMarkerSetFile(markerFile)
            self.__createMarkerHMMs(binMarkerSets, hmmModelFile, bReportProgress=False)
        elif markerFileType == BinMarkerSets.TREE_MARKER_SET:
            binIdToBinMarkerSets = self.parseLineageMarkerSetFile(markerFile)
            self.__createMarkerHMMs(binIdToBinMarkerSets[binId], hmmModelFile, bReportProgress=False)
        else:
            shutil.copyfile(markerFile, hmmModelFile)

        return hmmModelFile

    def __createLineageHmmModels(self, binIds, markerFile):
        """Create lineage-specific HMMs for each bin."""

        self.logger.info('  Extracting lineage-specific HMMs with %d threads:' % self.numThreads)

        workerQueue = mp.Queue()
        writerQueue = mp.Queue()

        for binId in binIds:
            workerQueue.put(binId)

        for _ in range(self.numThreads):
            workerQueue.put(None)

        binIdToModels = mp.Manager().dict()

        try:
            calcProc = [mp.Process(target=self.__fetchModelInfo, args=(binIdToModels, markerFile, workerQueue, writerQueue)) for _ in range(self.numThreads)]
            writeProc = mp.Process(target=self.__reportFetchProgress, args=(len(binIds), writerQueue))

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

        # create a standard dictionary from the managed dictionary
        d = {}
        for binId in binIdToModels.keys():
            d[binId] = binIdToModels[binId]

        return d

    def __fetchModelInfo(self, binIdToModels, markerFile, queueIn, queueOut):
        """Fetch HMM."""
        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            hmmModelFile = self.createHmmModelFile(binId, markerFile)

            modelParser = HmmModelParser(hmmModelFile)
            binIdToModels[binId] = modelParser.models()

            os.remove(hmmModelFile)

            queueOut.put(binId)

    def __reportFetchProgress(self, numBins, queueIn):
        """Report progress of extracted HMMs."""

        numProcessedBins = 0
        if self.logger.getEffectiveLevel() <= logging.INFO:
                statusStr = '    Finished extracting HMMs for %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        while True:
            binId = queueIn.get(block=True, timeout=None)
            if binId == None:
                break

            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                statusStr = '    Finished extracting HMMs for %d of %d (%.2f%%) bins.' % (numProcessedBins, numBins, float(numProcessedBins) * 100 / numBins)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

    def markerFileType(self, markerFile):
        """Determine type of marker file."""
        with open(markerFile, 'r') as f:
            header = f.readline()

        if DefaultValues.TAXON_MARKER_FILE_HEADER in header:
            return BinMarkerSets.TAXONOMIC_MARKER_SET
        elif DefaultValues.LINEAGE_MARKER_FILE_HEADER in header:
            return BinMarkerSets.TREE_MARKER_SET
        elif 'HMMER3' in header:
            return BinMarkerSets.HMM_MODELS_SET
        else:
            self.logger.error('Unrecognized file type: ' + markerFile)
            sys.exit()

    def __createMarkerHMMs(self, binMarkerSet, outputFile, bReportProgress=True):
        """Create HMM file for markers."""

        # get list of marker genes
        markerGenes = binMarkerSet.getMarkerGenes()

        # get all genes from the same clan as any marker gene
        pfam = PFAM(DefaultValues.PFAM_CLAN_FILE)
        genesInSameClan = pfam.genesInSameClan(markerGenes)

        # extract marker genes along with all genes from the same clan
        allMarkers = markerGenes | genesInSameClan

        if bReportProgress:
            self.logger.info("  There are %d genes in the marker set and %d genes from the same PFAM clan." % (len(markerGenes), len(genesInSameClan)))

        # create file with all model accession numbers
        keyFile = os.path.join(tempfile.gettempdir(), str(uuid.uuid4()))
        fout = open(keyFile, 'w')
        for modelAcc in allMarkers:
            fout.write(modelAcc + '\n')
        fout.close()

        # fetch specified models
        HF = HMMERRunner(mode='fetch')
        HF.fetch(DefaultValues.HMM_MODELS, keyFile, outputFile, bKeyFile=True)

        # index the HMM file
        if os.path.exists(outputFile + '.ssi'):
            os.remove(outputFile + '.ssi')
        HF.index(outputFile)

        # remove key file
        os.remove(keyFile)

    def parseTaxonomicMarkerSetFile(self, markerSetFile):
        """Parse marker set from a taxonomic-specific marker set file."""
        with open(markerSetFile) as f:
            f.readline()  # skip header

            binLine = f.readline()
            taxonId = binLine.split('\t')[0]
            binMarkerSets = BinMarkerSets(taxonId, BinMarkerSets.TAXONOMIC_MARKER_SET)
            binMarkerSets.read(binLine)

        return binMarkerSets

    def parseLineageMarkerSetFile(self, markerSetFile):
        """Parse marker sets from a lineage-specific marker set file."""

        # read all marker sets
        binIdToBinMarkerSets = {}
        with open(markerSetFile) as f:
            f.readline()  # skip header

            for line in f:
                lineSplit = line.split('\t')
                binId = lineSplit[0]

                binMarkerSets = BinMarkerSets(binId, BinMarkerSets.TREE_MARKER_SET)
                binMarkerSets.read(line)

                # determine selected marker set
                selectedMarkerSetMap = self.parseSelectedMarkerSetMap()
                binMarkerSets.setLineageSpecificSelectedMarkerSet(selectedMarkerSetMap)

                binIdToBinMarkerSets[binId] = binMarkerSets

        return binIdToBinMarkerSets

    def parseSelectedMarkerSetMap(self):
        selectedMarkerSetMap = {}
        for line in open(DefaultValues.SELECTED_MARKER_SETS):
            lineSplit = line.split('\t')
            internalID = lineSplit[0]
            selectedID = lineSplit[1].rstrip()

            selectedMarkerSetMap[internalID] = selectedID

        return selectedMarkerSetMap

    def writeBinModels(self, binIdToModels, filename):
        """Save HMM model info for each bin to file."""

        self.logger.info('  Saving HMM info to file.')

        with gzip.open(filename, 'wb') as output:
            pickle.dump(binIdToModels, output, pickle.HIGHEST_PROTOCOL)

    def loadBinModels(self, filename):
        """Read HMM model info for each bin from file."""

        self.logger.info('  Reading HMM info from file.')

        with gzip.open(filename, 'rb') as f:
            binIdToModels = pickle.load(f)

        return binIdToModels
