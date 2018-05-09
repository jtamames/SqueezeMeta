###############################################################################
#
# common.py - utility functions used in many places in CheckM
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
import errno
import sys
import ast
import logging

import numpy as np

from checkm.defaultValues import DefaultValues


def getBinIdsFromOutDir(outDir):
    """Get bin ids."""
    binIds = []

    binDir = os.path.join(outDir, 'bins')

    for f in os.listdir(binDir):
        if os.path.isdir(os.path.join(binDir, f)) and f != 'storage':
            binIds.append(f)

    return binIds


def readDistribution(prefix):
    """Read distribution file."""
    distFile = os.path.join(DefaultValues.DISTRIBUTION_DIR, prefix + '.txt')
    checkFileExists(distFile)

    with open(distFile, 'r') as f:
        s = f.read()
        d = ast.literal_eval(s)

    return d


def findNearest(array, value):
    '''Find nearest array element to a given value.'''
    idx = (np.abs(np.array(array) - value)).argmin()
    return array[idx]


def checkEmptyDir(inputDir):
    """Check the the specified directory is empty and create it if necessary."""
    if not os.path.exists(inputDir):
        makeSurePathExists(inputDir)
    else:
        # check if directory is empty
        files = os.listdir(inputDir)
        if len(files) != 0:
            logger = logging.getLogger()
            logger.warning('  [Warning] Output directory must be empty: ' + inputDir + '\n')
            sys.exit()


def checkFileExists(inputFile):
    """Check if file exists."""
    if not os.path.exists(inputFile):
        logger = logging.getLogger()
        logger.error('  [Error] Input file does not exists: ' + inputFile + '\n')
        sys.exit()


def checkDirExists(inputDir):
    """Check if directory exists."""
    if not os.path.exists(inputDir):
        logger = logging.getLogger()
        logger.error('  [Error] Input directory does not exists: ' + inputDir + '\n')
        sys.exit()


def makeSurePathExists(path):
    """Create directory if it does not exist."""
    if not path:
        return

    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            logger = logging.getLogger()
            logger.error('  [Error] Specified path does not exist: ' + path + '\n')
            sys.exit()


def binIdFromFilename(filename):
    """Extract bin id from bin filename."""
    binId = os.path.basename(filename)
    binId = os.path.splitext(binId)[0]

    return binId


def reassignStdOut(outFile):
    """Redirect standard out to a file."""
    oldStdOut = sys.stdout
    if(outFile != ''):
        try:
            # redirect stdout to a file
            sys.stdout = open(outFile, 'w')
        except:
            logger = logging.getLogger()
            logger.error("   [Error] Error diverting stdout to file: " + outFile)
            sys.exit()

    return oldStdOut


def restoreStdOut(outFile, oldStdOut):
    """Redirect standard out back to system standard out."""
    if(outFile != ''):
        try:
            # redirect stdout to a file
            sys.stdout.close()
            sys.stdout = oldStdOut
        except:
            logger = logging.getLogger()
            logger.error("   [Error] Error restoring stdout ", outFile)
            sys.exit()
