###############################################################################
#
# pplacer.py - runs pplacer and provides functions for parsing output
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
import stat
import shutil
import subprocess
import logging
from collections import defaultdict

from checkm.defaultValues import DefaultValues

from checkm.common import checkDirExists
from checkm.util.seqUtils import readFasta, writeFasta


class PplacerRunner():
    """Wrapper for running pplacer."""
    def __init__(self, threads):
        self.logger = logging.getLogger()
        self.numThreads = threads

        # make sure pplace and guppy are on the system path
        self.__checkForPplacer()
        self.__checkForGuppy()

    def run(self, binFiles, resultsParser, outDir, bReducedTree):
        # make sure output and tree directories exist
        checkDirExists(outDir)
        alignOutputDir = os.path.join(outDir, 'storage', 'tree')
        checkDirExists(alignOutputDir)

        treeFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_TREE_OUT)
        pplacerJsonOut = os.path.join(alignOutputDir, DefaultValues.PPLACER_JSON_OUT)
        pplacerOut = os.path.join(alignOutputDir, DefaultValues.PPLACER_OUT)

        # create concatenated alignment file for each bin
        concatenatedAlignFile = self.__createConcatenatedAlignment(binFiles, resultsParser, alignOutputDir)

        pplacerRefPkg = DefaultValues.PPLACER_REF_PACKAGE_FULL
        if bReducedTree:
            pplacerRefPkg = DefaultValues.PPLACER_REF_PACKAGE_REDUCED

        # check if concatenated alignment file is empty
        # (this can occur when all genomes have no phylogenetically informative marker genes)
        if os.stat(concatenatedAlignFile)[stat.ST_SIZE] == 0:
            self.logger.info('  No genomes were identified that could be placed in the reference genome tree.')
            shutil.copyfile(os.path.join(pplacerRefPkg, DefaultValues.GENOME_TREE), treeFile)
            return

        # run pplacer to place bins in reference genome tree
        self.logger.info('  Placing %d bins into the genome tree with pplacer (be patient).' % len(binFiles))
        cmd = 'pplacer -j %d -c %s -o %s %s > %s' % (self.numThreads,
                                                     pplacerRefPkg,
                                                     pplacerJsonOut,
                                                     concatenatedAlignFile,
                                                     pplacerOut)
        os.system(cmd)

        # extract tree
        cmd = 'guppy tog -o %s %s' % (treeFile, pplacerJsonOut)
        os.system(cmd)

    def __createConcatenatedAlignment(self, binFiles, resultsParser, alignOutputDir):
        """Create a concatenated alignment of marker genes for each bin."""

        # read alignment files
        self.logger.info('  Reading marker alignment files.')
        alignments = defaultdict(dict)
        files = os.listdir(alignOutputDir)
        binIds = set()
        for f in files:
            if f.endswith('.masked.faa'):
                markerId = f[0:f.find('.masked.faa')]
                seqs = readFasta(os.path.join(alignOutputDir, f))

                for seqId, seq in seqs.iteritems():
                    binId = seqId[0:seqId.find(DefaultValues.SEQ_CONCAT_CHAR)]

                    alignments[markerId][binId] = seq
                    binIds.add(binId)

        # get all markers and their lengths
        markerIds = resultsParser.models[resultsParser.models.keys()[0]].keys()
        markerIdLens = {}
        for markerId in markerIds:
            markerIdLens[markerId] = resultsParser.models[resultsParser.models.keys()[0]][markerId].leng

        # create concatenated alignment
        self.logger.info('  Concatenating alignments.')
        concatenatedSeqs = {}
        for markerId in sorted(markerIds):
            seqs = alignments[markerId]

            for binId in binIds:
                if binId in seqs:
                    # append alignment
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + seqs[binId]
                else:
                    # missing gene
                    concatenatedSeqs[binId] = concatenatedSeqs.get(binId, '') + '-' * markerIdLens[markerId]

        # save concatenated alignment
        concatenatedAlignFile = os.path.join(alignOutputDir, DefaultValues.PPLACER_CONCAT_SEQ_OUT)
        writeFasta(concatenatedSeqs, concatenatedAlignFile)

        return concatenatedAlignFile

    def __checkForPplacer(self):
        """Check to see if pplacer is on the system before we try to run it."""

        # Assume that a successful pplacer -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['pplacer', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure pplacer is on your system path.")
            sys.exit()

    def __checkForGuppy(self):
        """Check to see if guppy is on the system before we try to run it."""

        # Assume that a successful pplacer -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['guppy', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            self.logger.error("  [Error] Make sure guppy, which is part of the pplacer package, is on your system path.")
            sys.exit()
