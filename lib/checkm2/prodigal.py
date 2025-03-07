import os
import sys
import stat
import subprocess
import logging
import numpy as np
import shutil
import gzip
import tempfile

from checkm2 import sequenceClasses
from checkm2 import fileManager

'''Prodigal module taken from CheckM1.'''

#TODO dont use meta mode for prodigal if a translation table was provided
#TODO: take provided translation table


class ProdigalRunner():
    """Wrapper for running prodigal."""

    def __init__(self, out_dir, bin_file):

        self.file_basename = os.path.splitext(os.path.basename(bin_file))[0]

        # make sure prodigal is installed
        self.checkForProdigal()
        self.faa_directory = out_dir

        self.aaGeneFile = os.path.join(out_dir, "{}{}".format(self.file_basename, '.faa'))
        self.ntGeneFile = os.path.join(out_dir,  "{}{}".format(self.file_basename, '.fna'))
        self.gffFile = os.path.join(out_dir,  "{}{}".format(self.file_basename, '.gff'))

    def __calculate_N50(self, list_of_lengths):

        if np.array(list_of_lengths).mean() == 0:
            return 0

        tmp = []
        for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
        tmp.sort()

        if (len(tmp) % 2) == 0:
            median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
        else:
            median = tmp[int(len(tmp) / 2)]

        return median


    def run(self, query, supplied_coding_table=None):
        bNucORFs = True
        prodigal_input = query
                  
        # decompress archive input files                
        if prodigal_input.endswith('.gz'):
            tmp_dir = tempfile.mkdtemp()
            prodigal_input = os.path.join(tmp_dir, os.path.basename(prodigal_input[0:-3]) + '.fna')  
            with gzip.open(query, 'rb') as f_in:
                with open(prodigal_input, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                
        seqs = sequenceClasses.SeqReader().read_nucleotide_sequences(prodigal_input)

        totalBases = 0

        contig_lengths = []
        GC = 0
        AT = 0
        for seqId, seq in seqs.items():
            totalBases += len(seq)
            contig_lengths.append(len(seq))
            GC += sum(seq.upper().count(x) for x in ("G", "C"))
            AT += sum(seq.upper().count(x) for x in ("A", "T"))

        GC = float(GC/(AT + GC + 1))

        # call ORFs with different translation tables and select the one with the highest coding density
        tableCodingDensity = {}

        if supplied_coding_table is None:
            ttables_to_check = [4, 11]
        else:
            ttables_to_check = [supplied_coding_table]


        for translationTable in ttables_to_check:
            aaGeneFile = self.aaGeneFile + '.' + str(translationTable)
            ntGeneFile = self.ntGeneFile + '.' + str(translationTable)
            gffFile = self.gffFile + '.' + str(translationTable)

            # check if there is sufficient bases to calculate prodigal parameters
            if totalBases < 100000:
                procedureStr = 'meta'  # use best precalculated parameters
            else:
                procedureStr = 'single'  # estimate parameters from data

            if bNucORFs:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -d %s -i %s > %s' % (procedureStr,
                                                                                                  translationTable,
                                                                                                  aaGeneFile,
                                                                                                  ntGeneFile,
                                                                                                  prodigal_input,
                                                                                                  gffFile))
            else:
                cmd = ('prodigal -p %s -q -m -f gff -g %d -a %s -i %s > %s' % (procedureStr,
                                                                                            translationTable,
                                                                                            aaGeneFile,
                                                                                            prodigal_input,
                                                                                            gffFile))
            os.system(cmd)

            if not self.__areORFsCalled(aaGeneFile) and procedureStr == 'single':
                # prodigal will fail to learn a model if the input genome has a large number of N's
                # so try gene prediction with 'meta'
                cmd = cmd.replace('-p single', '-p meta')
                try:
                    os.system(cmd)
                except Exception as e:
                    logging.error('An error occured while running prodigal: {}'.format(e))
                    sys.exit(1)

            # determine coding density
            prodigalParser = ProdigalGeneFeatureParser(gffFile)

            codingBases = 0



            for seqId, seq in seqs.items():
                codingBases += prodigalParser.codingBases(seqId)

            if totalBases != 0:
                codingDensity = float(codingBases) / totalBases
            else:
                codingDensity = 0
            tableCodingDensity[translationTable] = codingDensity

        # determine best translation table

        if supplied_coding_table is not None:
            bestTranslationTable = supplied_coding_table
        else:
            bestTranslationTable = 11
            if (tableCodingDensity[4] - tableCodingDensity[11] > 0.05) and tableCodingDensity[4] > 0.7:
                bestTranslationTable = 4


        shutil.copyfile(self.aaGeneFile + '.' + str(bestTranslationTable), self.aaGeneFile)
        shutil.copyfile(self.gffFile + '.' + str(bestTranslationTable), self.gffFile)
        if bNucORFs:
            shutil.copyfile(self.ntGeneFile + '.' + str(bestTranslationTable), self.ntGeneFile)

        # clean up redundant prodigal results
        for translationTable in ttables_to_check:
            os.remove(self.aaGeneFile + '.' + str(translationTable))
            os.remove(self.gffFile + '.' + str(translationTable))
            if bNucORFs:
                os.remove(self.ntGeneFile + '.' + str(translationTable))

        os.remove(self.ntGeneFile)
        os.remove(self.gffFile)


        gene_lengths = []
        cds_count = 0
        aa_seqs = sequenceClasses.SeqReader().read_nucleotide_sequences(self.aaGeneFile)
        for seqId, seq in aa_seqs.items():
            gene_lengths.append(len(seq))
            cds_count += 1


#        if prodigal_input.endswith('.gz'):
#            shutil.rmtree(tmp_dir)

        maxContigLen = np.array(contig_lengths).max()
        totalContigs = len(contig_lengths)

        return self.file_basename, bestTranslationTable, tableCodingDensity[bestTranslationTable], \
               self.__calculate_N50(contig_lengths), np.array(gene_lengths).mean(), totalBases,\
               cds_count, GC, totalContigs, maxContigLen

    def __areORFsCalled(self, aaGeneFile):
        return os.path.exists(aaGeneFile) and os.stat(aaGeneFile)[stat.ST_SIZE] != 0

    def areORFsCalled(self, bNucORFs):
        # if requested, check if nucleotide gene sequences have been generated
        if bNucORFs:
            return os.path.exists(self.ntGeneFile) and os.stat(self.ntGeneFile)[stat.ST_SIZE] != 0

        # otherwise, only the amino acid gene sequences are required
        return os.path.exists(self.aaGeneFile) and os.stat(self.aaGeneFile)[stat.ST_SIZE] != 0

    def checkForProdigal(self):
        """Check to see if Prodigal is on the system before we try to run it."""

        # Assume that a successful prodigal -h returns 0 and anything
        # else returns something non-zero
        try:
            subprocess.call(['prodigal', '-h'], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        except:
            logging.error("Make sure prodigal is on your system path.")
            sys.exit(1)

class ProdigalFastaParser():
    """Parses prodigal FASTA output."""

    def __init__(self):
        pass

    def genePositions(self, filename):
        fileManager.check_if_file_exists(filename)

        gp = {}
        for line in open(filename):
            if line[0] == '>':
                lineSplit = line[1:].split()

                geneId = lineSplit[0]
                startPos = int(lineSplit[2])
                endPos = int(lineSplit[4])

                gp[geneId] = [startPos, endPos]

        return gp

class ProdigalGeneFeatureParser():
    """Parses prodigal FASTA output."""

    def __init__(self, filename):
        fileManager.check_if_file_exists(filename)

        self.genes = {}
        self.lastCodingBase = {}

        self.__parseGFF(filename)

        self.codingBaseMasks = {}
        for seqId in self.genes:
            self.codingBaseMasks[seqId] = self.__buildCodingBaseMask(seqId)

    def __parseGFF(self, filename):
        """Parse genes from GFF file."""
        self.translationTable = None
        for line in open(filename):
            if line.startswith('# Model Data') and not self.translationTable:
                lineSplit = line.split(';')
                for token in lineSplit:
                    if 'transl_table' in token:
                        self.translationTable = int(token[token.find('=') + 1:])

            if line[0] == '#' or line.strip() == '"':
                # work around for Prodigal having lines with just a
                # quotation on it when FASTA files have Windows style
                # line endings
                continue

            lineSplit = line.split('\t')
            seqId = lineSplit[0]
            if seqId not in self.genes:
                geneCounter = 0
                self.genes[seqId] = {}
                self.lastCodingBase[seqId] = 0

            geneId = seqId + '_' + str(geneCounter)
            geneCounter += 1

            start = int(lineSplit[3])
            end = int(lineSplit[4])

            self.genes[seqId][geneId] = [start, end]
            self.lastCodingBase[seqId] = max(self.lastCodingBase[seqId], end)

    def __buildCodingBaseMask(self, seqId):
        """Build mask indicating which bases in a sequences are coding."""

        # safe way to calculate coding bases as it accounts
        # for the potential of overlapping genes
        codingBaseMask = np.zeros(self.lastCodingBase[seqId])
        for pos in self.genes[seqId].values():
            codingBaseMask[pos[0]:pos[1] + 1] = 1

        return codingBaseMask

    def codingBases(self, seqId, start=0, end=None):
        """Calculate number of coding bases in sequence between [start, end)."""

        # check if sequence has any genes
        if seqId not in self.genes:
            return 0

        # set end to last coding base if not specified
        if end == None:
            end = self.lastCodingBase[seqId]

        return np.sum(self.codingBaseMasks[seqId][start:end])
