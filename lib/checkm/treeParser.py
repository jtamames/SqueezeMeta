###############################################################################
#
# treeParser.py - parse genome tree and associated tree metadata
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
import json

import dendropy

import checkm.prettytable as prettytable

from checkm.defaultValues import DefaultValues
from checkm.markerSets import MarkerSet, BinMarkerSets
from checkm.common import checkDirExists, reassignStdOut, restoreStdOut, getBinIdsFromOutDir
from checkm.util.seqUtils import readFasta
from checkm.util.taxonomyUtils import taxonomicPrefixes


class TreeParser():
    """Parse genome tree and associated tree metadata."""
    def __init__(self):
        self.logger = logging.getLogger()

    def printSummary(self, outputFormat, outDir, resultsParser, bTabTable, outFile, binStats):
        if outputFormat == 1:
            self.reportBinTaxonomy(outDir, resultsParser, bTabTable, outFile, binStats, bLineageStatistics=False)
        elif outputFormat == 2:
            self.reportBinTaxonomy(outDir, resultsParser, bTabTable, outFile, binStats, bLineageStatistics=True)
        elif outputFormat == 3:
            self.reportNewickTree(outDir, outFile, None)
        elif outputFormat == 4:
            self.reportNewickTree(outDir, outFile, 'taxonomy')
        elif outputFormat == 5:
            self.reportFullMSA(outDir, outFile)
        else:
            self.logger.error("Unknown output format: %d", outputFormat)

    def reportFullMSA(self, outDir, outFile):
        """Create MSA with all reference and bin alignments."""

        # write bin alignments to file
        oldStdOut = reassignStdOut(outFile)
        for line in open(os.path.join(outDir, 'storage', 'tree', DefaultValues.PPLACER_CONCAT_SEQ_OUT)):
            print(line.rstrip())

        # read duplicate seqs
        duplicateNodes = self.__readDuplicateSeqs()

        # write reference alignments to file
        seqs = readFasta(os.path.join(DefaultValues.PPLACER_REF_PACKAGE_FULL, DefaultValues.GENOME_TREE_FASTA))
        for seqId, seq in seqs.iteritems():
            print('>' + seqId)
            print(seq)

            if seqId in duplicateNodes:
                for dupSeqId in duplicateNodes[seqId]:
                    print('>' + dupSeqId)
                    print(seq)

        restoreStdOut(outFile, oldStdOut)

    def readPlacementFile(self, placementFile):
        '''Read pplacer JSON placement file.'''
        jsonData = open(placementFile)

        data = json.load(jsonData)
        binIdToPP = {}
        for placementData in data['placements']:
            binId = placementData['nm'][0][0]

            topPP = 0
            for pp in placementData['p']:
                if pp[2] > topPP:
                    topPP = pp[2]

            binIdToPP[binId] = topPP

        jsonData.close()

        return binIdToPP

    def reportNewickTree(self, outDir, outFile, leafLabels=None):
        # read duplicate nodes
        duplicateSeqs = self.__readDuplicateSeqs()

        # read tree
        treeFile = os.path.join(outDir, 'storage', 'tree', DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)

        # clean up internal node labels
        for node in tree.internal_nodes():
            if node.label:
                labelSplit = node.label.split('|')

                label = labelSplit[0]
                if labelSplit[1] != '':
                    label += '|' + labelSplit[1]
                if labelSplit[2] != '':
                    label += '|' + labelSplit[2]

                node.label = label

        # insert duplicate nodes into tree
        for leaf in tree.leaf_nodes():
            duplicates = duplicateSeqs.get(leaf.taxon.label, None)
            if duplicates != None:
                newParent = leaf.parent_node.new_child(edge_length=leaf.edge_length)
                curLeaf = leaf.parent_node.remove_child(leaf)
                newParent.new_child(taxon=curLeaf.taxon, edge_length=0)
                for d in duplicates:
                    newParent.new_child(taxon=dendropy.Taxon(label=d), edge_length=0)

        # append taxonomy to leaf nodes
        if leafLabels == 'taxonomy':
            # read taxonomy string for each IMG genome
            taxonomy = {}
            for line in open(os.path.join(DefaultValues.GENOME_TREE_DIR, DefaultValues.GENOME_TREE_TAXONOMY)):
                lineSplit = line.split('\t')
                taxonomy[lineSplit[0]] = lineSplit[1].rstrip()

            # append taxonomy to leaf labels
            for leaf in tree.leaf_nodes():
                taxaStr = taxonomy.get(leaf.taxon.label, None)
                if taxaStr:
                    leaf.taxon.label += '|' + taxaStr

        # write out tree
        oldStdOut = reassignStdOut(outFile)
        print(tree.as_string(schema='newick', suppress_rooting=True))
        restoreStdOut(outFile, oldStdOut)

    def getInsertionBranchId(self, outDir, binIds):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        # read genome tree (if it exists)
        binIdToUID = {}
        treeFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)

        # find first parent of each bin with a taxonomic label
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                binIdToUID[binId] = 'NA'
                continue

            # find first node decorated with a UID string between leaf and root
            parentNode = node.parent_node
            while parentNode != None:
                if parentNode.label:
                    uid = parentNode.label.split('|')[0]
                    break

                parentNode = parentNode.parent_node

            binIdToUID[binId] = uid

        return binIdToUID

    def getBinTaxonomy(self, outDir, binIds):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        # read genome tree (if it exists)
        binIdToTaxonomy = {}
        treeFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)

        # find first parent of each bin with a taxonomic label
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                binIdToTaxonomy[binId] = 'NA'
                continue

            # find first node decorated with a taxon string between leaf and root
            taxaStr = None
            parentNode = node.parent_node
            while parentNode != None:
                if parentNode.label:
                    tokens = parentNode.label.split('|')

                    if tokens[1] != '':
                        if taxaStr:
                            taxaStr = tokens[1] + ';' + taxaStr
                        else:
                            taxaStr = tokens[1]

                parentNode = parentNode.parent_node

            if not taxaStr:
                domainNode = self.__findDomainNode(node)
                taxaStr = domainNode.label.split('|')[1] + ' (root)'

            binIdToTaxonomy[node.taxon.label] = taxaStr

        return binIdToTaxonomy

    def __findDomainNode(self, binNode):
        """Find node defining the domain. Assumes 'binNode' is the leaf node of a bin on either the archaeal or bacterial branch."""

        # bin is either on the bacterial or archaeal branch descendant from the root,
        # so descend tree to first internal node with a label. Note that there may
        # be multiple bins inserted in the tree creating unlabeled internal nodes.
        # Identical population bins can also cause internal bifurcations.

        # find first parent node which contains at least one IMG genome
        curNode = binNode.parent_node
        while True:
            found_ref_genome = False
            for leaf in curNode.leaf_nodes():
                if leaf.taxon.label.startswith('IMG_'):
                    found_ref_genome = True
                    break
                    
            if found_ref_genome:
                break

            curNode = curNode.parent_node

        # perform depth first search to find a labeled internal node
        queue = [curNode]
        while queue:
            curNode = queue.pop(0)

            if curNode.label:
                return curNode

            for child in curNode.child_nodes():
                if child.is_internal():
                    queue.append(child)

        self.logger.error('  [Error] Failed to associate bin with a domain. Please report this bug.')
        sys.exit()

    def getBinSisterTaxonomy(self, outDir, binIds):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        # read genome tree
        treeFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)

        # read taxonomy string for each IMG genome
        leafIdToTaxonomy = {}
        for line in open(os.path.join(DefaultValues.GENOME_TREE_DIR, DefaultValues.GENOME_TREE_TAXONOMY)):
            lineSplit = line.split('\t')
            leafIdToTaxonomy[lineSplit[0]] = lineSplit[1].rstrip()

        # find LCA of all labeled node in sister lineage
        binIdToSisterTaxonomy = {}
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)

            taxaStr = ''
            if node != None:
                # get taxonomic labels of all internal nodes in sister lineages
                sisterNodes = node.sister_nodes()
                internalTaxonomyLabels = set()
                leafTaxonomyLabels = set()
                for sn in sisterNodes:
                    for curNode in sn.postorder_iter():
                        if curNode.is_leaf():
                            if curNode.taxon.label:
                                taxonomy = leafIdToTaxonomy.get(curNode.taxon.label, None)
                                if taxonomy != None:  # inserted bins will not have an assigned taxonomy
                                    for taxa in taxonomy.split(';'):
                                        leafTaxonomyLabels.add(taxa.strip())
                        else:
                            if curNode.label:
                                tokens = curNode.label.split('|')
                                if tokens[1] != '':
                                    for taxa in tokens[1].split(';'):
                                        internalTaxonomyLabels.add(taxa)

                # find LCA of taxonomic labels in rank order;
                # only consider leaf node labels if there were no internal labels
                labels = internalTaxonomyLabels
                if len(labels) == 0:
                    labels = leafTaxonomyLabels

                for prefix in taxonomicPrefixes:
                    taxa = []
                    for taxon in labels:
                        if prefix in taxon:
                            taxa.append(taxon)

                    if len(taxa) == 1:
                        # unambiguous label at this rank
                        taxaStr += taxa[0] + ';'
                    elif len(taxa) > 1:
                        # unable to resolve taxonomy at this rank
                        break

            if not taxaStr:
                taxaStr = 'unresolved'
            binIdToSisterTaxonomy[binId] = taxaStr

        return binIdToSisterTaxonomy

    def __getNextNamedNode(self, node, uniqueIdToLineageStatistics):
        """Get first parent node with taxonomy information."""
        parentNode = node.parent_node
        while True:
            if parentNode == None:
                break  # reached the root node so terminate

            if parentNode.label:
                trustedUniqueId = parentNode.label.split('|')[0]
                trustedStats = uniqueIdToLineageStatistics[trustedUniqueId]
                if trustedStats['taxonomy'] != '':
                    return trustedStats['taxonomy']

            parentNode = parentNode.parent_node

        return 'root'

    def __getMarkerSet(self, parentNode, tree, uniqueIdToLineageStatistics,
                                    numGenomesMarkers, bootstrap,
                                    bForceDomain, bRequireTaxonomy):
        """Get marker set for next parent node meeting selection criteria."""

        # ascend tree to root finding first node meeting all selection criteria
        selectedParentNode = parentNode
        taxonomyStr = 'root'
        while True:
            if selectedParentNode.label:  # nodes inserted by PPLACER will not have a label
                trustedUniqueId = selectedParentNode.label.split('|')[0]
                nodeTaxonomy = selectedParentNode.label.split('|')[1]

                stats = uniqueIdToLineageStatistics[trustedUniqueId]
                if stats['# genomes'] >= numGenomesMarkers and stats['bootstrap'] >= bootstrap:
                    if not bForceDomain or nodeTaxonomy in ['k__Bacteria', 'k__Archaea']:
                        if not bRequireTaxonomy or stats['taxonomy'] != '':
                            # get closest taxonomic label
                            taxonomyStr = stats['taxonomy']
                            if not bRequireTaxonomy and stats['taxonomy'] == '':
                                taxonomyStr = self.__getNextNamedNode(selectedParentNode, uniqueIdToLineageStatistics)

                            # all criteria meet, so use marker set from this node
                            break

            if selectedParentNode.parent_node == None:
                break  # reached the root node so terminate

            selectedParentNode = selectedParentNode.parent_node

        # get marker set meeting all criteria required for a trusted marker set
        taxonomyStr = taxonomyStr.split(';')[-1]  # extract most specific taxonomy identifier
        markerSet = MarkerSet(trustedUniqueId, taxonomyStr, int(stats['# genomes']), eval(stats['marker set']))

        return selectedParentNode, markerSet

    def __refineMarkerSet(self, markerSet, binNode, tree, uniqueIdToLineageStatistics, numGenomesRefine):
        """Refine marker set to account for lineage-specific gene loss and duplication."""

        # lineage-specific refine is done with the sister lineage to where the bin is inserted

        # get lineage-specific marker set which will be used to refine the above marker set
        curNode = binNode.sister_nodes()[0]
        while True:
            if curNode.label:  # nodes inserted by PPLACER will not have a label
                uniqueId = curNode.label.split('|')[0]
                stats = uniqueIdToLineageStatistics[uniqueId]

                if stats['# genomes'] >= numGenomesRefine:
                    break

            curNode = curNode.parent_node
            if curNode == None:
                break  # reached the root node so terminate

        # get lineage-specific marker set
        lineageMarkerSet = eval(stats['marker set'])

        # refine marker set by finding the intersection between these two sets,
        # this removes markers that are not single-copy or ubiquitous in the
        # specific lineage of a bin
        # Note: co-localization information is taken from the trusted set

        # get all lineage-specific marker genes
        allLineageSpecificGenes = set()
        for m in lineageMarkerSet:
            for gene in m:
                allLineageSpecificGenes.add(gene)

        # remove genes not present in the lineage-specific gene set
        finalMarkerSet = []
        for ms in markerSet.markerSet:
            s = set()
            for gene in ms:
                if gene in allLineageSpecificGenes:
                    s.add(gene)

            if s:
                finalMarkerSet.append(s)

        refinedMarkerSet = MarkerSet(markerSet.UID, markerSet.lineageStr, markerSet.numGenomes, finalMarkerSet)

        return refinedMarkerSet

    def __readLineageSpecificGenesToRemove(self):
        """Get set of genes subject to lineage-specific gene loss and duplication."""

        self.lineageSpecificGenesToRemove = {}
        for line in open(os.path.join(DefaultValues.GENOME_TREE_DIR, DefaultValues.GENOME_TREE_MISSING_DUPLICATE)):
            lineSplit = line.split('\t')
            uid = lineSplit[0]
            missingGenes = eval(lineSplit[1])
            duplicateGenes = eval(lineSplit[2])
            self.lineageSpecificGenesToRemove[uid] = missingGenes.union(duplicateGenes)

    def __removeInvalidLineageMarkerGenes(self, markerSet, lineageSpecificMarkersToRemove):
        """Refine marker set to account for lineage-specific gene loss and duplication."""

        # refine marker set by removing marker genes subject to lineage-specific
        # gene loss and duplication
        #
        # Note: co-localization information is taken from the trusted set

        finalMarkerSet = []
        for ms in markerSet.markerSet:
            s = set()
            for gene in ms:
                geneIdToTest = gene
                if geneIdToTest.startswith('PF'):
                    geneIdToTest = gene.replace('PF', 'pfam')
                    geneIdToTest = geneIdToTest[0:geneIdToTest.rfind('.')]

                if geneIdToTest not in lineageSpecificMarkersToRemove:
                    s.add(gene)

            if s:
                finalMarkerSet.append(s)

        refinedMarkerSet = MarkerSet(markerSet.UID, markerSet.lineageStr, markerSet.numGenomes, finalMarkerSet)

        return refinedMarkerSet

    def getBinMarkerSets(self, outDir, markerFile,
                                    numGenomesMarkers,
                                    bootstrap, bNoLineageSpecificRefinement,
                                    bForceDomain, bRequireTaxonomy,
                                    resultsParser, minUnique, maxMulti):
        """Determine marker sets for each bin."""

        self.logger.info('  Determining marker sets for each genome bin.')

        # get all bin ids
        binIds = getBinIdsFromOutDir(outDir)

        # get statistics for internal nodes
        uniqueIdToLineageStatistics = self.readNodeMetadata()

        # determine marker set for each bin
        treeFile = os.path.join(outDir, 'storage', 'tree', DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)
        rootNode = tree.find_node(filter_fn=lambda n: n.parent_node == None)

        fout = open(markerFile, 'w')
        fout.write(DefaultValues.LINEAGE_MARKER_FILE_HEADER + '\n')

        numProcessedBins = 0
        statusStr = ''
        for binId in binIds:
            if self.logger.getEffectiveLevel() <= logging.INFO:
                numProcessedBins += 1
                sys.stderr.write(' ' * len(statusStr) + '\r')  # clear previous line
                statusStr = '    Finished processing %d of %d (%.2f%%) bins (current: %s).' % (numProcessedBins, len(binIds), float(numProcessedBins) * 100 / len(binIds), binId)
                sys.stderr.write('%s\r' % statusStr)
                sys.stderr.flush()

            node = tree.find_node_with_taxon_label(binId)
            binMarkerSets = BinMarkerSets(binId, BinMarkerSets.TREE_MARKER_SET)
            if node == None:
                # bin is not in tree
                node, markerSet = self.__getMarkerSet(rootNode, tree, uniqueIdToLineageStatistics,
                                                        numGenomesMarkers, bootstrap,
                                                        bForceDomain, bRequireTaxonomy)
                binMarkerSets.addMarkerSet(markerSet)
            else:
                # special case: if node is on the bacterial or archaeal branch descendant from the root,
                # then move down the tree to include the domain-specific marker set
                parentNode = node.parent_node
                while parentNode != None:
                    if parentNode.label:
                        bRoot = (parentNode.parent_node == None)
                        break

                    parentNode = parentNode.parent_node

                if bRoot:
                    # since the root is the first labeled node, we need to descend the
                    # tree to incorporate the domain-specific marker set
                    domainNode = self.__findDomainNode(node)
                    curNode = domainNode.child_nodes()[0]
                else:
                    curNode = node

                # get lineage specific refinement for first node with an id
                if not bNoLineageSpecificRefinement:
                    uniqueId = parentNode.label.split('|')[0]
                    self.__readLineageSpecificGenesToRemove()
                    lineageSpecificRefinement = self.lineageSpecificGenesToRemove[uniqueId]

                # ascend tree to root, recording all marker sets meeting selection criteria
                while curNode.parent_node != None:
                    uniqueHits, multiCopyHits = resultsParser.results[binId].countUniqueHits()
                    tempForceDomain = bForceDomain or (uniqueHits < minUnique) or (multiCopyHits > maxMulti)

                    curNode, markerSet = self.__getMarkerSet(curNode.parent_node, tree, uniqueIdToLineageStatistics,
                                                                numGenomesMarkers, bootstrap,
                                                                tempForceDomain, bRequireTaxonomy)

                    if not bNoLineageSpecificRefinement:
                        markerSet = self.__removeInvalidLineageMarkerGenes(markerSet, lineageSpecificRefinement)

                    binMarkerSets.addMarkerSet(markerSet)

            binMarkerSets.write(fout)

        if self.logger.getEffectiveLevel() <= logging.INFO:
            sys.stderr.write('\n')

        fout.close()

    def readNodeMetadata(self):
        """Read metadata for internal nodes."""

        uniqueIdToLineageStatistics = {}
        metadataFile = os.path.join(DefaultValues.GENOME_TREE_DIR, DefaultValues.GENOME_TREE_METADATA)
        with open(metadataFile) as f:
            f.readline()
            for line in f:
                lineSplit = line.rstrip().split('\t')

                uniqueId = lineSplit[0]

                d = {}
                d['# genomes'] = int(lineSplit[1])
                d['taxonomy'] = lineSplit[2]
                try:
                    d['bootstrap'] = float(lineSplit[3])
                except:
                    d['bootstrap'] = 'NA'
                d['gc mean'] = float(lineSplit[4])
                d['gc std'] = float(lineSplit[5])
                d['genome size mean'] = float(lineSplit[6]) / 1e6
                d['genome size std'] = float(lineSplit[7]) / 1e6
                d['gene count mean'] = float(lineSplit[8])
                d['gene count std'] = float(lineSplit[9])
                d['marker set'] = lineSplit[10].rstrip()

                uniqueIdToLineageStatistics[uniqueId] = d

        return uniqueIdToLineageStatistics

    def readLineageMetadata(self, outDir, binIds):
        """Get metadata for each bin."""

        uniqueIdToLineageStatistics = self.readNodeMetadata()

        # read genome tree
        treeFile = os.path.join(outDir, 'storage', 'tree', DefaultValues.PPLACER_TREE_OUT)
        tree = dendropy.Tree.get_from_path(treeFile, schema='newick', rooting="force-rooted", preserve_underscores=True)

        # find first parent of each bin with a label
        binIdToLineageStatistics = {}
        for binId in binIds:
            node = tree.find_node_with_taxon_label(binId)
            if node == None:
                d = {}
                d['# genomes'] = 'NA'
                d['taxonomy'] = 'unresolved'
                d['gc mean'] = 'NA'
                d['gc std'] = 'NA'
                d['genome size mean'] = 'NA'
                d['genome size std'] = 'NA'
                d['gene count mean'] = 'NA'
                d['gene count std'] = 'NA'
                d['marker set'] = 'NA'
                binIdToLineageStatistics[binId] = d
                continue

            # find first labeled parent node (nodes inserted by pplacer will be unlabeled)
            parentNode = node.parent_node
            uniqueId = None
            while parentNode != None:
                if parentNode.label:
                    uniqueId = parentNode.label.split('|')[0]
                    break

                parentNode = parentNode.parent_node

            if uniqueId:
                binIdToLineageStatistics[binId] = uniqueIdToLineageStatistics[uniqueId]
            else:
                self.logger.error('Failed to find lineage-specific statistics for inserted bin: ' + node.taxon.label)
                sys.exit(0)

        return binIdToLineageStatistics

    def reportBinTaxonomy(self, outDir, resultsParser, bTabTable, outFile, binStats, bLineageStatistics):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        # get all bin ids
        binIds = getBinIdsFromOutDir(outDir)

        # get taxonomy for each bin
        binIdToTaxonomy = self.getBinTaxonomy(outDir, binIds)

        # write table
        if not bLineageStatistics:
            self.__printSimpleSummaryTable(binIdToTaxonomy, resultsParser, bTabTable, outFile)
        else:
            # get taxonomy of sister lineage for each bin
            binIdToSisterTaxonomy = self.getBinSisterTaxonomy(outDir, binIds)
            binIdToUID = self.getInsertionBranchId(outDir, binIds)

            binIdToLineageStatistics = self.readLineageMetadata(outDir, binIds)
            self.__printFullTable(binIdToUID, binIdToTaxonomy, binIdToSisterTaxonomy, binIdToLineageStatistics, resultsParser, binStats, bTabTable, outFile)

    def __printSimpleSummaryTable(self, binIdToTaxonomy, resultsParser, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        arbitraryBinId = binIdToTaxonomy.keys()[0]
        markerCountLabel = '# unique markers (of %d)' % len(resultsParser.models[arbitraryBinId])
        header = ['Bin Id', markerCountLabel, '# multi-copy', 'Taxonomy']

        if bTabTable:
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.align['Taxonomy'] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(binIdToTaxonomy.keys()):
            uniqueHits, multiCopyHits = resultsParser.results[binId].countUniqueHits()

            row = [binId, uniqueHits, multiCopyHits, binIdToTaxonomy[binId]]

            if bTabTable:
                print('\t'.join(map(str, row)))
            else:
                pTable.add_row(row)

        if not bTabTable:
            print(pTable.get_string(sortby=markerCountLabel, reversesort=True))

        # restore stdout
        restoreStdOut(outFile, oldStdOut)

    def __printFullTable(self, binIdToUID, binIdToTaxonomy, binIdToSisterTaxonomy, binIdToLineageStatistics, resultsParser, binStats, bTabTable, outFile):
        # redirect output
        oldStdOut = reassignStdOut(outFile)

        arbitraryBinId = binIdToTaxonomy.keys()[0]
        markerCountLabel = '# unique markers (of %d)' % len(resultsParser.models[arbitraryBinId])
        header = ['Bin Id', markerCountLabel, "# multi-copy"]
        header += ['Insertion branch UID', 'Taxonomy (contained)', 'Taxonomy (sister lineage)']
        header += ['GC', 'Genome size (Mbp)', 'Gene count', 'Coding density', 'Translation table']
        header += ['# descendant genomes', 'Lineage: GC mean', 'Lineage: GC std']
        header += ['Lineage: genome size (Mbp) mean', 'Lineage: genome size (Mbp) std']
        header += ['Lineage: gene count mean', 'Lineage: gene count std']

        if bTabTable:
            pTable = None
            print('\t'.join(header))
        else:
            pTable = prettytable.PrettyTable(header)
            pTable.float_format = '.2'
            pTable.float_format['GC'] = '.1'
            pTable.float_format['Lineage: GC mean'] = '.1'
            pTable.float_format['Lineage: GC std'] = '.1'
            pTable.float_format['Lineage: gene count mean'] = '.0'
            pTable.float_format['Lineage: gene count std'] = '.0'
            pTable.align = 'c'
            pTable.align[header[0]] = 'l'
            pTable.align['Insertion branch UID'] = 'l'
            pTable.align['Taxonomy (contained)'] = 'l'
            pTable.align['Taxonomy (sister lineage)'] = 'l'
            pTable.hrules = prettytable.FRAME
            pTable.vrules = prettytable.NONE

        for binId in sorted(binIdToTaxonomy.keys()):
            uniqueHits, multiCopyHits = resultsParser.results[binId].countUniqueHits()

            truncSisterLineage = binIdToSisterTaxonomy[binId]
            for taxa in binIdToTaxonomy[binId].split(';'):
                truncSisterLineage = truncSisterLineage.replace(taxa + ';', '')

            if len(truncSisterLineage) == 0:
                truncSisterLineage = 'unresolved'
            elif truncSisterLineage[-1] == ';':
                truncSisterLineage = truncSisterLineage[0:-1]

            row = [binId, uniqueHits, multiCopyHits]
            row += [binIdToUID[binId], binIdToTaxonomy[binId], truncSisterLineage]
            row += [binStats[binId]['GC'] * 100]
            row += [float(binStats[binId]['Genome size']) / 1e6]
            row += [binStats[binId]['# predicted genes']]
            row += [binStats[binId]['Coding density']]
            row += [binStats[binId]['Translation table']]
            row += [binIdToLineageStatistics[binId]['# genomes']]
            row += [binIdToLineageStatistics[binId]['gc mean']]
            row += [binIdToLineageStatistics[binId]['gc std']]
            row += [binIdToLineageStatistics[binId]['genome size mean']]
            row += [binIdToLineageStatistics[binId]['genome size std']]
            row += [binIdToLineageStatistics[binId]['gene count mean']]
            row += [binIdToLineageStatistics[binId]['gene count std']]

            if bTabTable:
                print('\t'.join(map(str, row)))
            else:
                pTable.add_row(row)

        if not bTabTable:
            print(pTable.get_string(sortby=markerCountLabel, reversesort=True))

        # restore stdout
        restoreStdOut(outFile, oldStdOut)

    def __readDuplicateSeqs(self):
        """Parse file indicating duplicate sequence alignments."""
        duplicateSeqs = {}
        for line in open(os.path.join(DefaultValues.GENOME_TREE_DIR, DefaultValues.GENOME_TREE_DEREP)):
            lineSplit = line.rstrip().split()
            if len(lineSplit) > 1:
                duplicateSeqs[lineSplit[0]] = lineSplit[1:]

        return duplicateSeqs
