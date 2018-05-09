###############################################################################
#
# seqUtils.py - Common functions for interacting with sequences
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
import gzip
import logging

nucleotide_bases = {'a', 'c', 'g', 't'}
insertion_bases = {'-', '.'}


def isNucleotide(seq_file, req_perc=0.9, max_seqs_to_read=10):
    """Check if a file contains sequences in nucleotide space.

    The check is performed by looking for the characters in
    {a,c,g,t,n,.,-} and confirming that these comprise the
    majority of a sequences. A set number of sequences are
    read and the file assumed to be not be in nucleotide space
    if none of these sequences are comprised primarily of the
    defined nucleotide set.

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.
    req_perc : float
        Percentage of bases in {a,c,g,t,n,.,-} before
        declaring the sequences as being in nucleotide
        space.
    max_seqs_to_read : int
        Maximum sequences to read before declaring
        sequence file to not be in nucleotide space.

    Returns
    -------
    boolean
        True is sequences are in nucleotide space, or file
        contains no sequences.
    """

    seqs = readFasta(seq_file)
    if len(seqs) == 0:
        return True

    seq_count = 0
    for _seq_id, seq in seqs.iteritems():
        seq = seq.lower()

        nt_bases = 0
        for c in (nucleotide_bases | {'n'} | insertion_bases):
            nt_bases += seq.count(c)

        if float(nt_bases) / len(seq) >= req_perc:
            return True

        seq_count += 1
        if seq_count == max_seqs_to_read:
            break

    return False


def queryYesNo(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    http://stackoverflow.com/questions/3041986/python-command-line-yes-no-input

    Parameters
    ----------
    question : str
        Prompt presented to the user.
    default : str
        Presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    Returns
    -------
    boolean
        True for "yes", False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def checkNuclotideSeqs(seq_files):
    """Check if files contain sequences in nucleotide space.

    Parameters
    ----------
    seq_files : iterable
        Sequence files to check.

    Returns
    -------
    boolean
        True if files can be treated as containing nucleotide sequences.
    """

    for seq_file in seq_files:
        if os.stat(seq_file).st_size == 0:
            continue
            
        if not isNucleotide(seq_file):
            print('Expected all files to contain sequences in nucleotide space.')
            print('File %s appears to contain amino acids sequences.' % seq_file)

            yes_response = queryYesNo('Do all files contain only nucleotide sequences?', default='no')
            if not yes_response:
                return False

    return True


def checkProteinSeqs(seq_files):
    """Check if files contain sequences in amino acid space.

    Parameters
    ----------
    seq_files : iterable
        Sequence files to check.

    Returns
    -------
    boolean
        True if files can be treated as containing amino acid sequences.
    """

    for seq_file in seq_files:
        if os.stat(seq_file).st_size == 0:
            continue
        
        if isNucleotide(seq_file):
            print('Expected all files to contain sequences in amino acid space.')
            print('File %s appears to contain nucleotide sequences.' % seq_file)

            yes_response = queryYesNo('Do all files contain only amino acid sequences?', default='no')
            if not yes_response:
                return False

    return True


def readFasta(fastaFile, trimHeader=True):
    '''Read sequences from FASTA file.'''
    try:
        if fastaFile.endswith('.gz'):
            openFile = gzip.open
        else:
            openFile = open

        seqs = {}
        for line in openFile(fastaFile):
            # skip blank lines
            if not line.strip():
                continue

            if line[0] == '>':
                if trimHeader:
                    seqId = line[1:].split(None, 1)[0]
                else:
                    seqId = line[1:].rstrip()
                seqs[seqId] = []
            else:
                seqs[seqId].append(line[0:-1])

        for seqId, seq in seqs.iteritems():
            seqs[seqId] = ''.join(seq)
    except:
        logger = logging.getLogger()
        logger.error("  [Error] Failed to process sequence file: " + fastaFile)
        sys.exit()

    return seqs


def readFastaSeqIds(fastaFile):
    '''Read sequence ids from FASTA file.'''
    if fastaFile.endswith('.gz'):
        openFile = gzip.open
    else:
        openFile = open

    seqIds = []
    for line in openFile(fastaFile):
        if line[0] == '>':
            seqId = line[1:].split(None, 1)[0]
            seqIds.append(seqId)

    return seqIds


def readFastaBases(fastaFile):
    '''Determine number of bases in FASTA file.'''
    if fastaFile.endswith('.gz'):
        openFile = gzip.open
    else:
        openFile = open

    bases = 0
    for line in openFile(fastaFile):
        if line[0] != '>':
            bases += len(line.rstrip())

    return bases


def readGenomicSeqsFromFasta(fastaFile, seqToIgnore=None):
    '''Read genomic sequences from FASTA file. Explicitly ignores sequences marked as plasmids.'''
    seqs = {}
    bRead = False
    for line in open(fastaFile):
        if line[0] == '>':
            if 'plasmid' in line.lower():
                bRead = False
            else:
                seqId = line[1:].split(None, 1)[0]
                seqs[seqId] = []
                bRead = True
        elif bRead:
            seqs[seqId].append(line[0:-1])

    for seqId, seq in seqs.iteritems():
        seqs[seqId] = ''.join(seq)

    return seqs


def writeFasta(seqs, outputFile):
    '''write sequences to FASTA file'''
    if outputFile.endswith('.gz'):
        fout = gzip.open(outputFile, 'wb')
    else:
        fout = open(outputFile, 'w')

    for seqId, seq in seqs.iteritems():
        fout.write('>' + seqId + '\n')
        fout.write(seq + '\n')
    fout.close()


def baseCount(seq):
    testSeq = seq.upper()
    a = testSeq.count('A')
    c = testSeq.count('C')
    g = testSeq.count('G')
    t = testSeq.count('T') + testSeq.count('U')

    return a, c, g, t


def calculateN50(seqLens):
    thresholdN50 = sum(seqLens) / 2.0

    seqLens.sort(reverse=True)

    testSum = 0
    for seqLen in seqLens:
        testSum += seqLen
        if testSum >= thresholdN50:
            N50 = seqLen
            break

    return N50
