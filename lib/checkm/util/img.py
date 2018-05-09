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
from collections import defaultdict

from checkm.util.seqUtils import readFasta
from checkm.util.taxonomyUtils import ranksByLabel


class IMG(object):
    genomeDir = '/srv/whitlam/bio/db/img//07042014/genomes/'
    pfamExtension = '.pfam.tab.txt'
    tigrExtension = '.tigrfam.tab.txt'

    def __init__(self, imgMetadataFile, redundantTIGRFAMsFile):
        self.logger = logging.getLogger()

        self.metadataFile = imgMetadataFile
        self.redundantTIGRFAMs = redundantTIGRFAMsFile

        self.cachedGenomeSeqLens = None
        self.cachedGenomeFamilyPositions = None
        self.cachedGenomeFamilyScaffolds = None

    def filterGenomeIds(self, genomeIds, metadata, fieldToFilterOn, valueToRetain):
        filteredGenomeIds = set()
        for genomeId in genomeIds:
            if metadata[genomeId][fieldToFilterOn] == valueToRetain:
                filteredGenomeIds.add(genomeId)

        return filteredGenomeIds

    def geneIdToScaffoldId(self, genomeId):
        d = {}

        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.gff'):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')

            scaffoldId = lineSplit[0]
            geneId = lineSplit[8]
            geneId = geneId[geneId.find('=') + 1:geneId.find(';')]

            d[geneId] = scaffoldId

        return d

    def pfamIdToGeneId(self, genomeId):
        return self.clusterIdToGeneId(genomeId, self.pfamExtension)

    def tigrIdToGeneId(self, genomeId):
        return self.clusterIdToGeneId(genomeId, self.tigrExtension)

    def clusterIdToGeneId(self, genomeId, extension):
        d = {}

        bHeader = True
        geneAnnotationFile = os.path.join(self.genomeDir, genomeId, genomeId + extension)
        for line in open(geneAnnotationFile):
            if bHeader:
                bHeader = False
                continue

            lineSplit = line.split('\t')
            geneId = lineSplit[0]
            clusterId = lineSplit[8]

            d[clusterId] = d.get(clusterId, []) + [geneId]

        return d

    def genomeMetadata(self):
        return self.genomeMetadataFromFile(self.metadataFile)

    def genomeMetadataFromFile(self, metadataFile):
        metadata = {}

        bHeader = True
        for line in open(metadataFile):
            lineSplit = line.split('\t')
            lineSplit = [x.strip() for x in lineSplit]

            if bHeader:
                statusIndex = lineSplit.index('Status')
                scaffoldCountIndex = lineSplit.index('Scaffold Count')
                gcIndex = lineSplit.index('GC Count')
                genomeSizeIndex = lineSplit.index('Genome Size')
                geneCountIndex = lineSplit.index('Gene Count')
                codingBaseCountIndex = lineSplit.index('Coding Base Count')
                bioticRelationshipsIndex = lineSplit.index('Biotic Relationships')
                n50Index = lineSplit.index('N50')

                domainIndex = lineSplit.index('Domain')
                phylumIndex = lineSplit.index('Phylum')
                classIndex = lineSplit.index('Class')
                orderIndex = lineSplit.index('Order')
                familyIndex = lineSplit.index('Family')
                genusIndex = lineSplit.index('Genus')
                speciesIndex = lineSplit.index('Species')

                bHeader = False
                continue

            genomeId = lineSplit[0].strip()

            rDomain = lineSplit[domainIndex].strip()
            rPhylum = lineSplit[phylumIndex].strip()
            rClass = lineSplit[classIndex].strip()
            rOrder = lineSplit[orderIndex].strip()
            rFamily = lineSplit[familyIndex].strip()
            rGenus = lineSplit[genusIndex].strip()
            rSpecies = lineSplit[speciesIndex].strip()

            metadata[genomeId] = {}
            metadata[genomeId]['status'] = lineSplit[statusIndex]
            metadata[genomeId]['taxonomy'] = [rDomain, rPhylum, rClass, rOrder, rFamily, rGenus, rSpecies]
            metadata[genomeId]['scaffold count'] = int(lineSplit[scaffoldCountIndex])

            try:
                metadata[genomeId]['GC Count'] = int(lineSplit[gcIndex])
                metadata[genomeId]['GC %'] = float(lineSplit[gcIndex]) / int(lineSplit[genomeSizeIndex])
            except:
                metadata[genomeId]['GC Count'] = 'NA'
                metadata[genomeId]['GC %'] = 'NA'

            try:
                metadata[genomeId]['genome size'] = int(lineSplit[genomeSizeIndex])
            except:
                metadata[genomeId]['genome size'] = 'NA'

            try:
                metadata[genomeId]['gene count'] = int(lineSplit[geneCountIndex])
            except:
                metadata[genomeId]['gene count'] = 'NA'

            try:
                metadata[genomeId]['coding base count'] = int(lineSplit[codingBaseCountIndex])
            except:
                metadata[genomeId]['coding base count'] = 'NA'

            metadata[genomeId]['biotic relationships'] = lineSplit[bioticRelationshipsIndex]
            metadata[genomeId]['N50'] = int(lineSplit[n50Index])

        return metadata

    def genomesWithMissingData(self, genomeIds):
        missingPFAM = self.missingPfamData(genomeIds)
        missingTIGR = self.missingTigrData(genomeIds)

        return missingPFAM.union(missingTIGR)

    def missingPfamData(self, genomeIds):
        missing = set()
        for genomeId in genomeIds:
            if not os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + self.pfamExtension):
                missing.add(genomeId)

                # if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.genes.fna'):
                #  print '[Warning] ' + genomeId + ' contains ORF data, but not PFAM annotations.'

        return missing

    def missingTigrData(self, genomeIds):
        missing = set()
        for genomeId in genomeIds:
            if not os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + self.tigrExtension):
                missing.add(genomeId)

                # if os.path.exists(IMG.genomeDir + genomeId + '/' + genomeId + '.genes.fna'):
                #  print '[Warning] ' + genomeId + ' contains ORF data, but not TIGRFAM annotations.'

        return missing

    def genomeIdsByTaxonomy(self, taxonStr, metadata):
        searchTaxa = taxonStr.split(';')
        genomeIdsOfInterest = set()
        for genomeId in metadata:
            bKeep = True
            for r in xrange(0, len(searchTaxa)):
                if taxonStr == 'universal':
                    bKeep = True
                elif taxonStr == 'prokaryotes' and (metadata[genomeId]['taxonomy'][0] == 'Bacteria' or metadata[genomeId]['taxonomy'][0] == 'Archaea'):
                    bKeep = True
                elif searchTaxa[r].strip() == metadata[genomeId]['taxonomy'][r]:
                    bKeep = True
                else:
                    bKeep = False
                    break

            if bKeep:
                genomeIdsOfInterest.add(genomeId)

        return genomeIdsOfInterest

    def getGenomesByClade(self, rank, clade, metadata):
        rankIndex = ranksByLabel[rank]
        genomeIdsOfInterest = set()
        for genomeId in metadata:
            if metadata[genomeId]['taxonomy'][rankIndex] == clade:
                genomeIdsOfInterest.add(genomeId)

        return genomeIdsOfInterest

    def lineageStats(self, metadata, mostSpecificRank):
        stats = {}
        for r in xrange(0, mostSpecificRank + 1):
            for _, data in metadata.iteritems():
                taxaStr = ';'.join(data['taxonomy'][0:r + 1])
                stats[taxaStr] = stats.get(taxaStr, 0) + 1

        return stats

    def lineagesSorted(self, metadata, mostSpecificRank=6):
        lineages = []
        for r in xrange(0, mostSpecificRank + 1):
            taxa = set()
            for _, data in metadata.iteritems():
                if 'unclassified' not in data['taxonomy'][0:r + 1]:
                    taxa.add(';'.join(data['taxonomy'][0:r + 1]))

            lineages += sorted(list(taxa))

        return lineages

    def lineagesByCriteria(self, metadata, minGenomes, mostSpecificRank):
        l = []

        stats = self.lineageStats(metadata, mostSpecificRank)
        for lineage in self.lineagesSorted(metadata, mostSpecificRank):
            if stats[lineage] > minGenomes:
                l.append(lineage)

        return l

    def __readTable(self, table, genomeIds, extension, clusterIdIndex):
        for genomeId in genomeIds:
            count = {}
            bHeader = True

            geneIdToFamilyIds = defaultdict(set)
            for line in open(os.path.join(self.genomeDir, genomeId, genomeId + extension)):
                if bHeader:
                    bHeader = False
                    continue

                lineSplit = line.split('\t')

                geneId = lineSplit[0]
                clusterId = lineSplit[clusterIdIndex]

                # IMG may annotate multiple parts of a gene as coming
                # from the same cluster (PFAM, TIGRFAM), but this should
                # only count as 1 gene having this annotation
                if clusterId not in geneIdToFamilyIds[geneId]:
                    geneIdToFamilyIds[geneId].add(clusterId)
                    count[clusterId] = count.get(clusterId, 0) + 1

            for clusterId, c in count.iteritems():
                if clusterId not in table:
                    table[clusterId] = {}
                table[clusterId][genomeId] = c

    def geneCountTable(self, genomeIds):
        table = {}
        self.__readTable(table, genomeIds, self.pfamExtension, 8)
        self.__readTable(table, genomeIds, self.tigrExtension, 6)

        return table

    def filterGeneCountTable(self, genomeIds, table, ubiquityThreshold=0.9, singleCopyThreshold=0.9):
        idsToFilter = []
        for pfamId, genomeCounts in table.iteritems():
            ubiquity = 0
            singleCopy = 0
            for genomeId in genomeIds:
                count = genomeCounts.get(genomeId, 0)

                if count > 0:
                    ubiquity += 1

                if count == 1:
                    singleCopy += 1

            if (float(ubiquity) / len(genomeIds) < ubiquityThreshold) or (float(singleCopy) / len(genomeIds) < singleCopyThreshold):
                idsToFilter.append(pfamId)

        for clusterId in idsToFilter:
            table.pop(clusterId)

        return table

    def __genomeIdToClusterScaffold(self, genomeId):
        """Determine position of PFAM and TIGRFAM genes in genome."""

        # determine mapping from gene ids to PFAM/TIGRFAM ids
        familyFile = os.path.join(self.genomeDir, genomeId, genomeId)
        pfamIdToGeneIds = self.familyIdToGeneId(familyFile + self.pfamExtension, 8)
        tigrIdToGeneIds = self.familyIdToGeneId(familyFile + self.tigrExtension, 6)

        # determine scaffold of genes from GFF file
        gffFile = os.path.join(self.genomeDir, genomeId, genomeId + '.gff')

        genePosition = {}
        for line in open(gffFile):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            if len(lineSplit) != 9:
                continue  # line likely annotates a CRISPR

            seqId = lineSplit[0]

            geneId = lineSplit[8].split(';')[0]
            geneId = geneId[geneId.find('=') + 1:]

            genePosition[geneId] = seqId

        # create gene mapping table
        try:
            # In theory, every PFAM or TIGRFAM gene identified should have
            # an entry in the GFF file and thus a position. In practice, there
            # are a few cases where this isn't tree (?) so only PFAMs/TIGRFAMs
            # with GFF entries are considered.
            familyIdToScaffoldIds = {}
            for pfamId, geneIds in pfamIdToGeneIds.iteritems():
                scaffolds = []
                for geneId in geneIds:
                    scaffold = genePosition.get(geneId, None)
                    if scaffold != None:
                        scaffolds.append(scaffold)

                if scaffolds:
                    familyIdToScaffoldIds[pfamId] = scaffolds

            for tigrId, geneIds in tigrIdToGeneIds.iteritems():
                scaffolds = []
                for geneId in geneIds:
                    scaffold = genePosition.get(geneId, None)
                    if scaffold != None:
                        scaffolds.append(scaffold)

                if scaffold:
                    familyIdToScaffoldIds[tigrId] = scaffolds
        except:
            print '[BUG]: __genomeIdToClusterScaffold'
            print sys.exc_info()[0]
            print genomeId, geneId, tigrId, pfamId
            sys.exit()

        return familyIdToScaffoldIds

    def precomputeGenomeFamilyScaffolds(self, genomeIds):
        """Cache scaffold of PFAM and TIGRFAM genes in genomes."""

        # This function is intended to speed up functions, such as geneDistTable(),
        # that are called multiple times (typically during simulations)
        self.cachedGenomeFamilyScaffolds = {}
        for genomeId in genomeIds:
            self.cachedGenomeFamilyScaffolds[genomeId] = self.__genomeIdToClusterScaffold(genomeId)

        return self.cachedGenomeFamilyScaffolds

    def familyIdToGeneId(self, filename, clusterIdIndex):
        """Determine gene ids associated with PFAMs or TIGRFAMs."""
        familyIdToGeneId = defaultdict(set)
        with open(filename) as f:
            f.readline()
            for line in f:
                lineSplit = line.split('\t')

                geneId = lineSplit[0]
                familyId = lineSplit[clusterIdIndex]
                familyIdToGeneId[familyId].update([geneId])

        return familyIdToGeneId

    def __genomeSeqLens(self, genomeId):
        """Determine length of contigs/scaffolds comprising genome."""
        genomeFile = os.path.join(self.genomeDir, genomeId, genomeId + '.fna')
        seqs = readFasta(genomeFile)

        seqLens = {}
        for seqId, seq in seqs.iteritems():
            seqLens[seqId] = len(seq)

        return seqLens

    def precomputeGenomeSeqLens(self, genomeIds):
        """Cache the length of contigs/scaffolds for all genomes."""

        # This function is intended to speed up functions, such as geneDistTable(),
        # that are called multiple times (typically during simulations)

        self.cachedGenomeSeqLens = {}
        for genomeId in genomeIds:
            self.cachedGenomeSeqLens[genomeId] = self.__genomeSeqLens(genomeId)

        return self.cachedGenomeSeqLens

    def __genomeFamilyPositions(self, genomeId, seqLens, spacingBetweenContigs):
        """Determine position of PFAM and TIGRFAM genes in genome."""

        # determine mapping from gene ids to PFAM/TIGRFAM ids
        familyFile = os.path.join(self.genomeDir, genomeId, genomeId)
        pfamIdToGeneIds = self.familyIdToGeneId(familyFile + self.pfamExtension, 8)
        tigrIdToGeneIds = self.familyIdToGeneId(familyFile + self.tigrExtension, 6)

        # determine position of genes from GFF file
        gffFile = os.path.join(self.genomeDir, genomeId, genomeId + '.gff')

        genePosition = {}
        contigStart = 0
        curSeqId = None
        for line in open(gffFile):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            if len(lineSplit) != 9:
                continue  # line likely annotates a CRISPR

            # check if we've moved to the next contig
            if curSeqId == None:
                curSeqId = lineSplit[0]

            if curSeqId != lineSplit[0]:
                contigStart += spacingBetweenContigs + seqLens[curSeqId]
                curSeqId = lineSplit[0]

            geneId = lineSplit[8].split(';')[0]
            geneId = geneId[geneId.find('=') + 1:]

            start = int(lineSplit[3])
            end = int(lineSplit[4])

            genePosition[geneId] = [contigStart + start, contigStart + end]

        # create gene mapping table
        try:
            # In theory, every PFAM or TIGRFAM gene identified should have
            # an entry in the GFF file and thus a position. In practice, there
            # are a few cases where this isn't tree (?) so only PFAMs/TIGRFAMs
            # with GFF entries are considered.
            familyIdToGenomePositions = {}
            for pfamId, geneIds in pfamIdToGeneIds.iteritems():
                positions = []
                for geneId in geneIds:
                    position = genePosition.get(geneId, None)
                    if position != None:
                        positions.append(position)

                if positions:
                    familyIdToGenomePositions[pfamId] = positions

            for tigrId, geneIds in tigrIdToGeneIds.iteritems():
                positions = []
                for geneId in geneIds:
                    position = genePosition.get(geneId, None)
                    if position != None:
                        positions.append(position)

                if positions:
                    familyIdToGenomePositions[tigrId] = positions
        except:
            print '[BUG]: __genomeFamilyPositions'
            print sys.exc_info()[0]
            print genomeId, geneId, tigrId, pfamId
            sys.exit()

        return familyIdToGenomePositions

    def precomputeGenomeFamilyPositions(self, genomeIds, spacingBetweenContigs):
        """Cache position of PFAM and TIGRFAM genes in genomes."""

        # This function is intended to speed up functions, such as geneDistTable(),
        # that are called multiple times (typically during simulations)
        self.cachedGenomeFamilyPositions = {}
        for genomeId in genomeIds:
            self.cachedGenomeFamilyPositions[genomeId] = self.__genomeFamilyPositions(genomeId, self.cachedGenomeSeqLens[genomeId], spacingBetweenContigs)

    def geneDistTable(self, genomeIds, markerGenes, spacingBetweenContigs=0):
        """Create table indicating position of each marker gene in a genome."""

        # Note: genomes split into multiple contigs are treated as contiguous,
        # with a spacing between contigs as specified

        table = {}
        for genomeId in genomeIds:
            # read length of scaffolds/contigs in genome
            if self.cachedGenomeSeqLens:
                seqLens = self.cachedGenomeSeqLens[genomeId]
            else:
                seqLens = self.__genomeSeqLens(genomeId)

            # read position of protein families on genome
            if self.cachedGenomeFamilyPositions:
                genomeFamilyPositions = self.cachedGenomeFamilyPositions[genomeId]
            else:
                genomeFamilyPositions = self.__genomeFamilyPositions(genomeId, seqLens, spacingBetweenContigs)

            # create marker gene position table for genome
            clusterIdToGenomePositions = {}
            for markerGene in markerGenes:
                positions = genomeFamilyPositions.get(markerGene, None)
                if positions != None:
                    clusterIdToGenomePositions[markerGene] = genomeFamilyPositions[markerGene]

            table[genomeId] = clusterIdToGenomePositions

        return table

    def identifyMitochondrialChloroplastGenes(self, genomeId):
        # identify mitochondrial or chloroplast sequences
        mitoChloroSeqs = set()
        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.fna'):
            if line[0] == '>':
                if 'mitochondria' in line.lower() or 'mitochondrion' in line.lower() or 'chloroplast' in line.lower():
                    mitoChloroSeqs.add(line[1:].split()[0])

        # identify mitochondrial or chloroplast genes
        mitoChloroGenes = set()
        for line in open(self.genomeDir + genomeId + '/' + genomeId + '.gff'):
            if line[0] == '#':
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId in mitoChloroSeqs:
                desc = lineSplit[8]
                geneId = desc.split(';')[0].split('=')[1]
                mitoChloroGenes.add(geneId)

        return mitoChloroGenes

    def identifyRedundantPFAMs(self, markerGenes):
        pfamIdToTigrId = defaultdict(list)
        for line in open(self.redundantTIGRFAMs):
            lineSplit = line.split('\t')
            pfamId = lineSplit[0]
            tigrId = lineSplit[1].rstrip()

            pfamIdToTigrId[pfamId].append(tigrId)

        pfamToRemove = set()
        for markerGene in markerGenes:
            if markerGene in pfamIdToTigrId:
                for tigrId in pfamIdToTigrId[markerGene]:
                    if tigrId in markerGenes:
                        pfamToRemove.add(markerGene)

        return pfamToRemove

    def identifyRedundantTIGRFAMs(self, markerGenes):
        tigrIdToPfamId = {}
        for line in open(self.redundantTIGRFAMs):
            lineSplit = line.split('\t')
            pfamId = lineSplit[0]
            tigrId = lineSplit[1].rstrip()

            tigrIdToPfamId[tigrId] = tigrIdToPfamId.get(tigrId, []) + [pfamId]

        tigrToRemove = set()
        for markerGene in markerGenes:
            if markerGene in tigrIdToPfamId:
                for pfamId in tigrIdToPfamId[markerGene]:
                    if pfamId in markerGenes:
                        tigrToRemove.add(markerGene)

        return tigrToRemove
