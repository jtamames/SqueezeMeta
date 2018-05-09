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

import numpy as np

from checkm.defaultValues import DefaultValues
from checkm.timeKeeper import TimeKeeper
from checkm.markerSets import MarkerSetParser, BinMarkerSets
from checkm.resultsParser import ResultsParser
from checkm.hmmerAligner import HmmerAligner
from checkm.markerGeneFinder import MarkerGeneFinder
from checkm.pplacer import PplacerRunner
from checkm.treeParser import TreeParser
from checkm.taxonParser import TaxonParser
from checkm.aminoAcidIdentity import AminoAcidIdentity
from checkm.binComparer import BinComparer
from checkm.binUnion import BinUnion
from checkm.binStatistics import BinStatistics
from checkm.coverage import Coverage
from checkm.coverageWindows import CoverageWindows
from checkm.genomicSignatures import GenomicSignatures
from checkm.unbinned import Unbinned
from checkm.merger import Merger
from checkm.profile import Profile
from checkm.binTools import BinTools
from checkm.ssuFinder import SSU_Finder
from checkm.PCA import PCA
from checkm.common import makeSurePathExists, checkFileExists, binIdFromFilename, reassignStdOut, restoreStdOut, getBinIdsFromOutDir, checkDirExists, checkEmptyDir

from checkm.plot.gcPlots import GcPlots
from checkm.plot.codingDensityPlots import CodingDensityPlots
from checkm.plot.tetraDistPlots import TetraDistPlots
from checkm.plot.distributionPlots import DistributionPlots
from checkm.plot.nxPlot import NxPlot
from checkm.plot.cumulativeLengthPlot import CumulativeLengthPlot
from checkm.plot.lengthHistogram import LengthHistogram
from checkm.plot.markerGenePosPlot import MarkerGenePosPlot
from checkm.plot.parallelCoordPlot import ParallelCoordPlot
from checkm.plot.binQAPlot import BinQAPlot
from checkm.plot.pcaPlot import PcaPlot
from checkm.plot.gcBiasPlots import GcBiasPlot

from checkm.util.seqUtils import checkNuclotideSeqs, checkProteinSeqs

from checkm.checkmData import DBManager

from checkm.test.test_ecoli import VerifyEcoli


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger()
        self.timeKeeper = TimeKeeper()

    def updateCheckM_DB(self, options):
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - data] Check for database updates. [%s]' % options.action[0])
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        DBM = DBManager()
        DBM.runAction(options.action)

    def binFiles(self, binFolder, binExtension):
        binFiles = []
        if binFolder is not None:
            all_files = os.listdir(binFolder)
            for f in all_files:
                if f.endswith(binExtension):
                    binFile = os.path.join(binFolder, f)
                    if os.stat(binFile).st_size == 0:
                        self.logger.warning("  [Warning] Skipping bin %s as it has a size of 0 bytes." % f)
                    else:
                        binFiles.append(binFile)

        if not binFiles:
            self.logger.error("  [Error] No bins found. Check the extension (-x) used to identify bins.")
            sys.exit()

        return sorted(binFiles)

    def tree(self, options):
        """Tree command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tree] Placing bins in reference genome tree.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        binFiles = self.binFiles(options.bin_folder, options.extension)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        # setup directory structure
        checkEmptyDir(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage'))

        # find phylogenetically informative genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                                 options.out_folder,
                                 DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                 DefaultValues.HMMER_PHYLO_OUT,
                                 DefaultValues.PHYLO_HMM_MODELS,
                                 options.bKeepAlignment,
                                 options.bNucORFs,
                                 options.bCalledGenes)

        # write model information to file
        markerSetParser = MarkerSetParser(options.threads)
        hmmModelInfoFile = os.path.join(options.out_folder, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        markerSetParser.writeBinModels(binIdToModels, hmmModelInfoFile)

        # calculate statistics for each genome bin
        self.logger.info('')
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.out_folder, DefaultValues.BIN_STATS_PHYLO_OUT)

        # align identified marker genes
        self.logger.info('')
        HA = HmmerAligner(options.threads)
        resultsParser = HA.makeAlignmentToPhyloMarkers(options.out_folder,
                                                            DefaultValues.PHYLO_HMM_MODELS,
                                                            DefaultValues.HMMER_TABLE_PHYLO_OUT,
                                                            binIdToModels,
                                                            False,
                                                            DefaultValues.E_VAL,
                                                            DefaultValues.LENGTH,
                                                            False,
                                                            os.path.join(options.out_folder, 'storage', 'tree')
                                                            )

        # place bins into genome tree
        self.logger.info('')
        pplacer = PplacerRunner(threads=options.pplacer_threads)  # fix at one thread to keep memory requirements reasonable
        pplacer.run(binFiles, resultsParser, options.out_folder, options.bReducedTree)

        self.timeKeeper.printTimeStamp()

    def treeQA(self, options):
        """QA command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tree_qa] Assessing phylogenetic markers found in each bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.tree_folder)

        # set HMM file for each bin
        markerSetParser = MarkerSetParser()
        hmmModelInfoFile = os.path.join(options.tree_folder, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        # calculate marker gene statistics
        RP = ResultsParser(binIdToModels)
        binStats = RP.analyseResults(options.tree_folder,
                                          DefaultValues.BIN_STATS_PHYLO_OUT,
                                          DefaultValues.HMMER_TABLE_PHYLO_OUT)

        # determine taxonomy of each bin
        self.logger.info('')
        treeParser = TreeParser()
        treeParser.printSummary(options.out_format, options.tree_folder, RP, options.bTabTable, options.file, binStats)

        if options.file != '':
            self.logger.info('  QA information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def lineageSet(self, options, db=None):
        """Lineage set command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - lineage_set] Inferring lineage-specific marker sets.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.tree_folder)

        # set HMM file for each bin
        markerSetParser = MarkerSetParser()
        hmmModelInfoFile = os.path.join(options.tree_folder, 'storage', DefaultValues.PHYLO_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        # calculate marker gene statistics
        resultsParser = ResultsParser(binIdToModels)
        resultsParser.analyseResults(options.tree_folder,
                                      DefaultValues.BIN_STATS_PHYLO_OUT,
                                      DefaultValues.HMMER_TABLE_PHYLO_OUT)

        # These options are incompatible with how the lineage-specific marker set is selected, so
        # the default values are currently hard-coded
        self.logger.info('')
        options.num_genomes_markers = 2
        options.bootstrap = 0
        options.bRequireTaxonomy = False

        treeParser = TreeParser()
        treeParser.getBinMarkerSets(options.tree_folder, options.marker_file,
                                    options.num_genomes_markers,
                                    options.bootstrap, options.bNoLineageSpecificRefinement,
                                    options.bForceDomain, options.bRequireTaxonomy,
                                    resultsParser, options.unique, options.multi)

        self.logger.info('')
        self.logger.info('  Marker set written to: ' + options.marker_file)

        self.timeKeeper.printTimeStamp()

    def taxonList(self, options, db=None):
        """Lineage set command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - taxon_list] Listing available taxonomic-specific marker sets.')
        self.logger.info('*******************************************************************************')

        taxonParser = TaxonParser()
        taxonParser.list(options.rank)

        self.timeKeeper.printTimeStamp()

    def taxonSet(self, options, db=None):
        """Taxon set command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - taxon_set] Generate taxonomic-specific marker set.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        path = os.path.split(options.marker_file)[0]
        if path:
            makeSurePathExists(path)

        taxonParser = TaxonParser()
        bValidSet = taxonParser.markerSet(options.rank, options.taxon, options.marker_file)

        if bValidSet:
            self.logger.info('')
            self.logger.info('  Marker set written to: ' + options.marker_file)

        self.timeKeeper.printTimeStamp()

    def analyze(self, options, db=None):
        """Analyze command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - analyze] Identifying marker genes in bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        binFiles = self.binFiles(options.bin_folder, options.extension)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        # setup directory structure
        makeSurePathExists(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'aai_qa'))

        # find marker genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                                 options.out_folder,
                                 DefaultValues.HMMER_TABLE_OUT,
                                 DefaultValues.HMMER_OUT,
                                 options.marker_file,
                                 options.bKeepAlignment,
                                 options.bNucORFs,
                                 options.bCalledGenes)

        markerSetParser = MarkerSetParser(options.threads)
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.out_folder,
                                                             getBinIdsFromOutDir(options.out_folder),
                                                             options.marker_file)

        hmmModelInfoFile = os.path.join(options.out_folder, 'storage', DefaultValues.CHECKM_HMM_MODEL_INFO)
        markerSetParser.writeBinModels(binIdToModels, hmmModelInfoFile)

        self.timeKeeper.printTimeStamp()

        # HMM model file
        if markerSetParser.markerFileType(options.marker_file) == BinMarkerSets.HMM_MODELS_SET:
            markerFile = options.marker_file
        else:
            markerFile = DefaultValues.HMM_MODELS

        # align marker genes with multiple hits within a bin
        self.logger.info('')
        HA = HmmerAligner(options.threads)
        HA.makeAlignmentsOfMultipleHits(options.out_folder,
                                          markerFile,
                                          DefaultValues.HMMER_TABLE_OUT,
                                          binIdToModels,
                                          binIdToBinMarkerSets,
                                          False,
                                          DefaultValues.E_VAL,
                                          DefaultValues.LENGTH,
                                          os.path.join(options.out_folder, 'storage', 'aai_qa')
                                          )

        self.timeKeeper.printTimeStamp()

        # calculate statistics for each genome bin
        self.logger.info('')
        binStats = BinStatistics(options.threads)
        binStats.calculate(binFiles, options.out_folder, DefaultValues.BIN_STATS_OUT)

        self.timeKeeper.printTimeStamp()

        # align top hit to each marker if requested
        if options.bAlignTopHit:
            alignmentOutputFolder = os.path.join(options.out_folder, 'storage', 'alignments')
            makeSurePathExists(alignmentOutputFolder)

            self.logger.info('')
            HA = HmmerAligner(options.threads)
            resultsParser = HA.makeAlignmentTopHit(options.out_folder,
                                                                options.marker_file,
                                                                DefaultValues.HMMER_TABLE_OUT,
                                                                binIdToModels,
                                                                False,
                                                                DefaultValues.E_VAL,
                                                                DefaultValues.LENGTH,
                                                                True,
                                                                alignmentOutputFolder
                                                                )

            # report marker gene data
            fout = open(os.path.join(alignmentOutputFolder, 'alignment_info.tsv'), 'w')
            fout.write('Marker Id\tLength (bp)\n')
            markerIds = resultsParser.models[resultsParser.models.keys()[0]].keys()
            for markerId in markerIds:
                fout.write('%s\t%d\n' % (markerId, resultsParser.models[resultsParser.models.keys()[0]][markerId].leng))
            fout.close()

            self.logger.info('')
            self.logger.info('  Alignments to top hits stored in: ' + alignmentOutputFolder)

            self.timeKeeper.printTimeStamp()

    def qa(self, options):
        """QA command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - qa] Tabulating genome statistics.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.analyze_folder)

        if options.exclude_markers:
            checkFileExists(options.exclude_markers)

        # calculate AAI between marks with multiple hits in a single bin
        aai = AminoAcidIdentity()
        aai.run(options.aai_strain, options.analyze_folder, options.alignment_file)

        # get HMM file for each bin
        self.logger.info('')
        markerSetParser = MarkerSetParser(options.threads)

        hmmModelInfoFile = os.path.join(options.analyze_folder, 'storage', DefaultValues.CHECKM_HMM_MODEL_INFO)
        binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.analyze_folder,
                                                             getBinIdsFromOutDir(options.analyze_folder),
                                                             options.marker_file,
                                                             options.exclude_markers)

        # get results for each bin
        RP = ResultsParser(binIdToModels)
        RP.analyseResults(options.analyze_folder,
                          DefaultValues.BIN_STATS_OUT,
                          DefaultValues.HMMER_TABLE_OUT,
                          bIgnoreThresholds=options.bIgnoreThresholds,
                          evalueThreshold=options.e_value,
                          lengthThreshold=options.length,
                          bSkipPseudoGeneCorrection=options.bSkipPseudoGeneCorrection,
                          bSkipAdjCorrection=options.bSkipAdjCorrection
                          )

        self.logger.info('')
        RP.printSummary(options.out_format, 
                        aai, binIdToBinMarkerSets, 
                        options.bIndividualMarkers, 
                        options.coverage_file, 
                        options.bTabTable, 
                        options.file, 
                        anaFolder=options.analyze_folder)
        RP.cacheResults(options.analyze_folder, 
                            binIdToBinMarkerSets, 
                            options.bIndividualMarkers)

        if options.file != '':
            self.logger.info('  QA information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def gcPlot(self, options):
        """GC plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - gc_plot] Creating GC histogram and delta-GC plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        plots = GcPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, options.distributions)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.gc_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def codingDensityPlot(self, options):
        """Coding density plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - coding_plot] Creating coding density histogram and delta-CD plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        plots = CodingDensityPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting coding density plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, options.distributions)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.coding_density_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def tetraDistPlot(self, options):
        """Tetranucleotide distance plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra_plot] Creating tetra-distance histogram and delta-TD plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)

        plots = TetraDistPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting tetranuclotide distance plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)

            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def distributionPlots(self, options):
        """Reference distribution plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - dist_plot] Creating GC, CD, and TD distribution plots.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        genomicSignatures = GenomicSignatures(K=4, threads=1)
        tetraSigs = genomicSignatures.read(options.tetra_profile)

        plots = DistributionPlots(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting reference distribution plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            binId = binIdFromFilename(f)
            plots.plot(f, tetraSigs, options.distributions)

            outputFile = os.path.join(options.plot_folder, binId) + '.ref_dist_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def tetraPcaPlot(self, options):
        """PCA plot of tetranucleotide signatures"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra_pca] Creating PCA plot of tetranucleotide signatures.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        self.logger.info('  Computing PCA of tetranuclotide signatures.\n')
        pca = PCA()
        seqIds, pc, variance = pca.pcaFile(options.tetra_profile, fraction=1.0, bCenter=True, bScale=False)

        plots = PcaPlot(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting PCA of tetranuclotide signatures for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, seqIds, pc, variance)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.tetra_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def coveragePcaPlot(self, options):
        """PCA plot of coverage profiles"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - cov_pca] Creating PCA plot of coverage profiles.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        checkFileExists(options.coverage_file)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        coverage = Coverage(threads=1)
        coverageStats = coverage.parseCoverage(options.coverage_file)

        seqIds = []
        coverageProfiles = []
        for binId, seqDict in coverageStats.iteritems():
            for seqId, bamDict in seqDict.iteritems():
                seqIds.append(seqId)

                coverages = []
                for _, coverage in bamDict.iteritems():
                    coverages.append(coverage)

                coverageProfiles.append(coverages)

        coverageProfiles = np.array(coverageProfiles)
        if coverageProfiles.shape[1] < 2:
            self.logger.error('  [Error] Coverage profile is 1 dimensional. PCA requires at least 2 dimensions.')
            sys.exit()

        self.logger.info('  Computing PCA of coverage profiles.\n')
        pca = PCA()
        pc, variance = pca.pcaMatrix(coverageProfiles, fraction=1.0, bCenter=True, bScale=False)

        plots = PcaPlot(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting PCA of coverage profiles for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, seqIds, pc, variance)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.cov_pca_plots.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def gcBiasPlot(self, options):
        """GC bias plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - gc_bias_plot] Plotting bin coverage as a function of GC.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        coverageWindows = CoverageWindows(options.threads)
        coverageProfile = coverageWindows.run(binFiles, options.bam_file, options.all_reads, options.min_align, options.max_edit_dist, options.window_size)

        plots = GcBiasPlot(options)
        filesProcessed = 1
        for f in binFiles:
            self.logger.info('  Plotting GC plots for %s (%d of %d)' % (f, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plots.plot(f, coverageProfile)

            binId = binIdFromFilename(f)
            outputFile = os.path.join(options.plot_folder, binId) + '.gc_bias_plot.' + options.image_type
            plots.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def nxPlot(self, options):
        """Nx-plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - nx_plot] Creating Nx-plots.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        nx = NxPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting Nx-plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            nx.plot(f)

            outputFile = os.path.join(options.plot_folder, binId) + '.nx_plot.' + options.image_type
            nx.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def cumulativeLengthPlot(self, options):
        """Cumulative sequence length plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - len_plot] Creating cumulative sequence length plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        plot = CumulativeLengthPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting cumulative sequence length plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            plot.plot(f)

            outputFile = os.path.join(options.plot_folder, binId) + '.seq_len_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def lengthHistogram(self, options):
        """Sequence length histogram command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - len_hist] Creating sequence length histogram.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        plot = LengthHistogram(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting sequence length histogram for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1
            plot.plot(f)

            outputFile = os.path.join(options.plot_folder, binId) + '.len_hist.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def markerPlot(self, options):
        """Marker gene position plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - marker_plot] Creating marker gene position plot.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        # generate plot for each bin
        binFiles = self.binFiles(options.bin_folder, options.extension)

        resultsParser = ResultsParser(None)
        markerGeneStats = resultsParser.parseMarkerGeneStats(options.out_folder)
        binStats = resultsParser.parseBinStatsExt(options.out_folder)

        plot = MarkerGenePosPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting marker gene position plot for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1

            if binId not in markerGeneStats or binId not in binStats:
                continue  # bin has no marker genes

            bPlotted = plot.plot(f, markerGeneStats[binId], binStats[binId])

            if bPlotted:
                outputFile = os.path.join(options.plot_folder, binId) + '.marker_pos_plot.' + options.image_type
                plot.savePlot(outputFile, dpi=options.dpi)
                self.logger.info('    Plot written to: ' + outputFile)
            else:
                self.logger.info('    No marker genes found in bin.')

        self.timeKeeper.printTimeStamp()

    def parallelCoordPlot(self, options):
        """Parallel coordinate plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - par_plot] Creating parallel coordinate plot of GC and coverage.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)
        checkFileExists(options.coverage_file)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        # read coverage stats file
        coverage = Coverage(threads=1)
        coverageStats = coverage.parseCoverage(options.coverage_file)

        # calculate sequence stats for all bins
        self.logger.info('  Calculating sequence statistics for each bin.')
        binStats = BinStatistics()
        seqStats = {}
        for f in binFiles:
            binId = binIdFromFilename(f)
            seqStats[binId] = binStats.sequenceStats(options.out_folder, f)

        # create plot for each bin
        self.logger.info('')
        plot = ParallelCoordPlot(options)
        filesProcessed = 1
        for f in binFiles:
            binId = binIdFromFilename(f)
            self.logger.info('  Plotting parallel coordinates for %s (%d of %d)' % (binId, filesProcessed, len(binFiles)))
            filesProcessed += 1

            plot.plot(binId, seqStats, coverageStats)

            outputFile = os.path.join(options.plot_folder, binId) + '.paralel_coord_plot.' + options.image_type
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('    Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def binQAPlot(self, options):
        """Bin QA plot command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - bin_qa_plot] Creating bar plot of bin quality.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(options.plot_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        # read model info
        # hmmModelInfoFile = os.path.join(options.analyze_folder, 'storage', DefaultValues.CHECKM_HMM_MODEL_INFO)
        # binIdToModels = markerSetParser.loadBinModels(hmmModelInfoFile)

        # read sequence stats file
        resultsParser = ResultsParser(None)
        binStatsExt = resultsParser.parseBinStatsExt(options.out_folder)

        # create plot for each bin
        plot = BinQAPlot(options)
        bMakePlot = True
        if not options.bIgnoreHetero:
            aai = AminoAcidIdentity()
            aai.run(options.aai_strain, options.out_folder, None)
            bMakePlot = plot.plot(binFiles, binStatsExt, options.bIgnoreHetero, aai.aaiHetero)
        else:
            bMakePlot = plot.plot(binFiles, binStatsExt, options.bIgnoreHetero, None)

        if bMakePlot:
            outputFile = os.path.join(options.plot_folder, 'bin_qa_plot.' + options.image_type)
            plot.savePlot(outputFile, dpi=options.dpi)
            self.logger.info('')
            self.logger.info('  Plot written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def unbinned(self, options):
        """Unbinned Command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - unbinned] Identify unbinned sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        unbinned = Unbinned()
        unbinned.run(binFiles, options.seq_file, options.output_seq_file, options.output_stats_file, options.min_seq_len)

        self.logger.info('')
        self.logger.info('  Unbinned sequences written to: ' + options.output_seq_file)
        self.logger.info('  Unbinned sequences statistics written to: ' + options.output_stats_file)

        self.timeKeeper.printTimeStamp()

    def coverage(self, options):
        """Coverage command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - coverage] Calculating coverage of sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        makeSurePathExists(os.path.dirname(options.output_file))

        binFiles = self.binFiles(options.bin_folder, options.extension)

        coverage = Coverage(options.threads)
        coverage.run(binFiles, options.bam_files, options.output_file, options.all_reads,
                        options.min_align, options.max_edit_dist, options.min_qc)

        self.logger.info('  Coverage information written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def tetraSignatures(self, options):
        """Tetranucleotide signature command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - tetra] Calculating tetranucleotide signature of sequences.')
        self.logger.info('*******************************************************************************')

        checkFileExists(options.seq_file)
        makeSurePathExists(os.path.dirname(options.output_file))

        self.logger.info('')
        tetraSig = GenomicSignatures(4, options.threads)
        tetraSig.calculate(options.seq_file, options.output_file)

        self.logger.info('  Tetranucletoide signatures written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def profile(self, options):
        """Profile command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - profile] Calculating percentage of reads mapped to each bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkFileExists(options.coverage_file)

        profile = Profile()
        profile.run(options.coverage_file, options.file, options.bTabTable)

        if options.file != '':
            self.logger.info('  Profile information written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def merge(self, options):
        """Merge command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - merge] Identifying bins with complementary sets of marker genes.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)

        binFiles = self.binFiles(options.bin_folder, options.extension)

        if not options.bCalledGenes:
            if not checkNuclotideSeqs(binFiles):
                return
        else:
            if not checkProteinSeqs(binFiles):
                return

        markerSetParser = MarkerSetParser()
        if markerSetParser.markerFileType(options.marker_file) == BinMarkerSets.TREE_MARKER_SET:
            self.logger.error('  [Error] Merge command requires a taxonomic-specific marker set or a user-defined HMM file.\n')
            return

        # setup directory structure
        makeSurePathExists(options.out_folder)
        makeSurePathExists(os.path.join(options.out_folder, 'bins'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage'))
        makeSurePathExists(os.path.join(options.out_folder, 'storage', 'hmms'))

        binIds = []
        for binFile in binFiles:
            binIds.append(binIdFromFilename(binFile))

        # find marker genes in genome bins
        mgf = MarkerGeneFinder(options.threads)
        binIdToModels = mgf.find(binFiles,
                         options.out_folder,
                         "merger.table.txt",
                         "merger.hmmer3",
                         options.marker_file,
                         False,
                         False,
                         options.bCalledGenes)

        # get HMM file for each bin
        markerSetParser = MarkerSetParser()
        binIdToBinMarkerSets = markerSetParser.getMarkerSets(options.out_folder, binIds, options.marker_file)

        # compare markers found in each bin
        self.logger.info('')
        merger = Merger()
        outputFile = merger.run(binFiles, options.out_folder, "merger.table.txt", binIdToModels, binIdToBinMarkerSets,
                                options.delta_comp, options.delta_cont, options.merged_comp, options.merged_cont)

        self.logger.info('\n  Merger information written to: ' + outputFile)

        self.timeKeeper.printTimeStamp()

    def outliers(self, options):
        """Outlier command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - outlier] Identifying outliers in bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder)
        checkFileExists(options.tetra_profile)
        makeSurePathExists(os.path.dirname(options.output_file))

        binFiles = self.binFiles(options.bin_folder, options.extension)

        binTools = BinTools()
        binTools.identifyOutliers(options.out_folder, binFiles, options.tetra_profile, options.distributions, options.report_type, options.output_file)

        self.logger.info('\n  Outlier information written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def joinTables(self, options):
        """Join tables command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - join_tables] Joining tables containing bin information.')
        self.logger.info('*******************************************************************************')

        # read all tables
        headers = {}
        rows = defaultdict(dict)
        binIds = set()
        for f in options.tables:
            with open(f) as fin:
                headers[f] = [x.strip() for x in fin.readline().split('\t')][1:]

                for line in fin:
                    lineSplit = [x.strip() for x in line.split('\t')]

                    binId = lineSplit[0]
                    binIds.add(binId)

                    for i, header in enumerate(headers[f]):
                        rows[binId][header] = lineSplit[i + 1]

        # write merge table
        oldStdOut = reassignStdOut(options.file)

        row = 'Bin Id'
        for f in options.tables:
            row += '\t' + '\t'.join(headers[f])
        print(row)

        for binId in binIds:
            row = binId
            for f in options.tables:
                for header in headers[f]:
                    row += '\t' + rows[binId].get(header, '')
            print(row)

        restoreStdOut(options.file, oldStdOut)

        if options.file:
            self.logger.info('\n  Joined table written to: ' + options.file)

        self.timeKeeper.printTimeStamp()

    def modify(self, options):
        """Modify command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [CheckM - modify] Modifying sequences in bin.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        makeSurePathExists(os.path.dirname(options.output_file))

        if not (options.add or options.remove or options.outlier_file):
            self.logger.warning('  [Warning] No modification to bin requested.\n')
            sys.exit()

        if (options.add or options.remove) and options.outlier_file:
            self.logger.warning("  [Warning] The 'outlier_file' option cannot be specified with 'add' or 'remove'.\n")
            sys.exit()

        binTools = BinTools()

        if options.add or options.remove:
            binTools.modify(options.bin_file, options.seq_file, options.add, options.remove, options.output_file)
        elif options.outlier_file:
            binTools.removeOutliers(options.bin_file, options.outlier_file, options.output_file)

        self.logger.info('  Modified bin written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - unique] Ensuring no sequences are assigned to multiple bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        binFiles = self.binFiles(options.bin_folder, options.extension)

        binTools = BinTools()
        binTools.unique(binFiles)

        self.timeKeeper.printTimeStamp()

    def ssuFinder(self, options):
        """SSU finder command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - ssu_finder] Identifying SSU (16S/18S) rRNAs in sequences.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        binFiles = self.binFiles(options.bin_folder, options.extension)

        checkFileExists(options.seq_file)
        makeSurePathExists(options.out_folder)

        ssuFinder = SSU_Finder(options.threads)
        ssuFinder.run(options.seq_file, binFiles, options.out_folder, options.evalue, options.concatenate)

        self.timeKeeper.printTimeStamp()

    def binCompare(self, options):
        """Bin compare command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - bin_compare] Comparing two sets of bins.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        checkDirExists(options.bin_folder1)
        checkDirExists(options.bin_folder2)

        binFiles1 = self.binFiles(options.bin_folder1, options.extension1)
        binFiles2 = self.binFiles(options.bin_folder2, options.extension2)

        binComparer = BinComparer()
        binComparer.report(binFiles1, binFiles2, options.seq_file, options.output_file)

        self.logger.info('')
        self.logger.info('  Detailed bin comparison written to: ' + options.output_file)

        self.timeKeeper.printTimeStamp()

    def binUnion(self, options):
        """Bin union command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - bin_union] Redundancy reduce multiple sets of bins into a single set.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        output_dir = options.output_dir
        makeSurePathExists(output_dir)

        bin_folders = []
        checkmQaTsvs = []
        for i, arg in enumerate(options.bin_or_checkm_qa_table):
            if i % 2 == 0:
                checkDirExists(arg)
                bin_folders.append(arg)
            else:
                checkFileExists(arg)
                checkmQaTsvs.append(arg)

        if len(bin_folders) < 2:
            self.logger.error("   [Error] Need to specify at least two bin folders, found %i: " % len(bin_folders))
            sys.exit()
        if len(bin_folders) != len(checkmQaTsvs):
            self.logger.error("   [Error] Need to specify the same number of bin folders as checkm_qa_tsv files, found %i and %i, respectively: " % (len(bin_folders), len(checkmQaTsvs)))
            sys.exit()

        binFileSets = []
        for bin_folder in bin_folders:
            self.logger.info("   Reading fasta files with extension %s from bin folder %s" % (options.extension, bin_folder))
            binFileSets.append(self.binFiles(bin_folder, options.extension))

        binUnion = BinUnion()

        contigConflictsOutputFile = os.path.join(output_dir, 'contigConflicts.csv')
        unionBinOutputFile = os.path.join(output_dir, 'union.txt')
        binUnion.report(bin_folders, binFileSets, checkmQaTsvs, unionBinOutputFile, contigConflictsOutputFile, options.min_completeness, options.max_contamination)

        self.logger.info('')

    def test(self, options):
        """Quick test of CheckM"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[CheckM - Test] Processing E.coli K12-W3310 to verify operation of CheckM.')
        self.logger.info('*******************************************************************************')
        self.logger.info('')

        verifyEcoli = VerifyEcoli()
        verifyEcoli.run(self, options.output_dir)

        self.timeKeeper.printTimeStamp()

    def parseOptions(self, options):
        """Parse user options and call the correct pipeline(s)"""
        try:
            if options.bVerbose:
                logging.basicConfig(format='', level=logging.DEBUG)
            elif options.bQuiet:
                logging.basicConfig(format='', level=logging.ERROR)
            else:
                logging.basicConfig(format='', level=logging.INFO)
        except:
            logging.basicConfig(format='', level=logging.INFO)

        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == "data"):
            self.updateCheckM_DB(options)
        elif(options.subparser_name == 'tree'):
            self.tree(options)
        elif(options.subparser_name == 'tree_qa'):
            self.treeQA(options)
        elif(options.subparser_name == 'lineage_set'):
            self.lineageSet(options)
        elif(options.subparser_name == 'taxon_list'):
            self.taxonList(options)
        elif(options.subparser_name == 'taxon_set'):
            self.taxonSet(options)
        elif(options.subparser_name == 'analyze'):
            self.analyze(options)
        elif(options.subparser_name == 'qa'):
            self.qa(options)
        elif(options.subparser_name == 'lineage_wf'):
            options.marker_file = os.path.join(options.out_folder, 'lineage.ms')
            options.tree_folder = options.out_folder
            options.analyze_folder = options.out_folder
            options.out_format = 1
            options.bAlignTopHit = False
            options.exclude_markers = None
            options.coverage_file = None

            self.tree(options)
            self.lineageSet(options)
            self.analyze(options)
            self.qa(options)
        elif(options.subparser_name == 'taxonomy_wf'):
            options.marker_file = os.path.join(options.out_folder, options.taxon + '.ms')
            options.analyze_folder = options.out_folder
            options.out_format = 1
            options.bAlignTopHit = False
            options.exclude_markers = None

            self.taxonSet(options)
            self.analyze(options)
            self.qa(options)
        elif(options.subparser_name == 'gc_plot'):
            self.gcPlot(options)
        elif(options.subparser_name == 'coding_plot'):
            self.codingDensityPlot(options)
        elif(options.subparser_name == 'tetra_plot'):
            self.tetraDistPlot(options)
        elif(options.subparser_name == 'dist_plot'):
            self.distributionPlots(options)
        elif(options.subparser_name == 'nx_plot'):
            self.nxPlot(options)
        elif(options.subparser_name == 'len_plot'):
            self.cumulativeLengthPlot(options)
        elif(options.subparser_name == 'len_hist'):
            self.lengthHistogram(options)
        elif(options.subparser_name == 'marker_plot'):
            self.markerPlot(options)
        elif(options.subparser_name == 'par_plot'):
            self.parallelCoordPlot(options)
        elif(options.subparser_name == 'tetra_pca'):
            self.tetraPcaPlot(options)
        elif(options.subparser_name == 'cov_pca'):
            self.coveragePcaPlot(options)
        elif(options.subparser_name == 'gc_bias_plot'):
            self.gcBiasPlot(options)
        elif(options.subparser_name == 'bin_qa_plot'):
            self.binQAPlot(options)
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif(options.subparser_name == 'coverage'):
            self.coverage(options)
        elif(options.subparser_name == 'tetra'):
            self.tetraSignatures(options)
        elif(options.subparser_name == 'profile'):
            self.profile(options)
        elif(options.subparser_name == 'merge'):
            self.merge(options)
        elif(options.subparser_name == 'outliers'):
            self.outliers(options)
        elif(options.subparser_name == 'join_tables'):
            self.joinTables(options)
        elif(options.subparser_name == 'modify'):
            self.modify(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'ssu_finder'):
            self.ssuFinder(options)
        elif(options.subparser_name == 'bin_compare'):
            self.binCompare(options)
        elif(options.subparser_name == 'bin_union'):
            self.binUnion(options)
        elif(options.subparser_name == 'test'):
            options.bCalledGenes = False
            options.pplacer_threads = 1
            self.test(options)
        else:
            self.logger.error('  [Error] Unknown CheckM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
