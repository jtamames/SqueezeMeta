###############################################################################
#
# hmmer.py - runs HMMER and provides functions for parsing output
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
import logging
from re import split as re_split

from biolib.external.execute import check_on_path


class FormatError(BaseException):
    pass


class HMMERError(BaseException):
    pass


class HMMMERModeError(BaseException):
    pass


class HmmModelError(Exception):
    pass


class HMMER():
    """Wrapper for running HMMER3."""
    def __init__(self, mode="dom"):
        self.logger = logging.getLogger('timestamp')

        # make sure HMMER is installed
        check_on_path('hmmsearch')

        # set the mode
        if mode == "dom":
            self.mode = 'domtblout'
        elif mode == "tbl":
            self.mode = 'tblout'
        elif mode == 'align':
            self.mode = 'align'
        elif mode == 'fetch':
            self.mode = 'fetch'
        else:
            raise HMMMERModeError("Mode %s not understood" % mode)

    def search(self, db, query, tableOut, hmmerOut, cmdlineOptions='', bKeepOutput=True):
        """Run hmmsearch"""
        # make the output dir and files
        if self.mode != 'domtblout' and self.mode != 'tblout':
            raise HMMMERModeError("Mode %s not compatible with search" % self.mode)

        if not bKeepOutput:
            hmmerOut = '/dev/null'

        cmd = ('hmmsearch --%s %s %s %s %s > %s' % (self.mode, tableOut, cmdlineOptions, db, query, hmmerOut))
        os.system(cmd)

    def align(self, db, query, outputFile, writeMode='>', outputFormat='PSIBLAST', trim=True):
        """Run hmmalign"""
        if self.mode != 'align':
            raise HMMMERModeError("Mode %s not compatible with align" % self.mode)

        # build up the option string
        opts = ''
        if trim:
            opts += ' --trim '

        # run hmmer
        cmd = 'hmmalign --allcol %s --outformat %s %s %s %s %s' % (opts, outputFormat, db, query, writeMode, outputFile)
        rtn = os.system(cmd)
        if rtn == 256:
            # assume user has a newer version of HMMER (>= 3.1b1) and the allcol parameter is no longer valid
            cmd = cmd.replace('--allcol', '')
            os.system(cmd)

    def fetch(self, db, key, fetchFileName, bKeyFile=False):
        """Run hmmfetch"""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        keyFileOpt = ''
        if bKeyFile:
            keyFileOpt = '-f'

        os.system('hmmfetch ' + keyFileOpt + ' %s %s > %s' % (db, key, fetchFileName))

    def index(self, hmmModelFile):
        """Index a HMM file."""
        if self.mode != 'fetch':
            raise HMMMERModeError("Mode %s not compatible with fetch" % self.mode)

        os.system('hmmfetch --index %s > /dev/null' % hmmModelFile)


class HMMERParser():
    """Parses tabular output."""
    def __init__(self, fileHandle, mode='dom'):
        self.handle = fileHandle
        if mode == 'dom':
            self.mode = 'domtblout'
        elif mode == 'tbl':
            self.mode = 'tblout'
        else:
            raise HMMERError("Mode %s not understood, please use 'dom' or 'tbl'" % mode)

    def __next__(self):
        """Process each hit in the file."""
        while 1:
            if self.mode == 'domtblout':
                hit = self.read_hits_domain_table()
            elif self.mode == 'tblout':
                hit = self.read_hits_table()
            else:
                raise HMMERError("Mode %s not understood" % self.mode)

            if hit == {}:
                return None
            else:
                return hit

    def read_hits_table(self):
        """Process single hit in tblout format."""
        """
We expect line to look like:
NODE_110054_length_1926_cov_24.692627_41_3 -          Ribosomal_S9         PF00380.14   5.9e-48  158.7   0.0   6.7e-48  158.5   0.0   1.0   1   0   0   1   1   1   1 # 1370 # 1756 # 1 # ID=41_3;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None
        """
        while (1):
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split(r'\s+', line.rstrip())
                    if len(dMatch) < 19:
                        raise FormatError("Error processing line:\n%s" % (line))
                    refined_match = dMatch[0:18] + [" ".join([str(i) for i in dMatch[18:]])]
                    return HmmerHitTBL(refined_match)
            except IndexError:
                return {}

    def read_hits_domain_table(self):
        """Process single hit in domtblout format."""
        """
We expect the line to look like:
NODE_925902_length_6780_cov_18.428171_754_2 -            399 PGK                  PF00162.14   384  2.2e-164  543.7   0.1   1   1  1.3e-167  2.5e-164  543.5   0.1     1   384     9   386     9   386 1.00 # 1767 # 2963 # -1 # ID=754_2;partial=00;start_type=ATG;rbs_motif=AGGA;rbs_spacer=5-10bp
        """
        while (1):
            line = self.handle.readline().rstrip()
            try:
                if line[0] != '#' and len(line) != 0:
                    dMatch = re_split(r'\s+', line.rstrip())
                    if len(dMatch) < 23:
                        raise FormatError("Error processing line:\n%s" % (line))
                    refined_match = dMatch[0:22] + [" ".join([str(i) for i in dMatch[22:]])]
                    return HmmerHitDOM(refined_match)
            except IndexError:
                return {}


class HmmerHitTBL():
    """Encapsulate a HMMER hit given in tblout format."""
    def __init__(self, values):
        if len(values) == 19:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.query_name = values[2]

            self.query_accession = values[3]
            if self.query_accession == '-':
                self.query_accession = self.query_name

            self.full_e_value = float(values[4])
            self.full_score = float(values[5])
            self.full_bias = float(values[6])
            self.best_e_value = float(values[7])
            self.best_score = float(values[8])
            self.best_bias = float(values[9])
            self.exp = float(values[10])
            self.reg = int(values[11])
            self.clu = int(values[12])
            self.ov = int(values[13])
            self.env = int(values[14])
            self.dom = int(values[15])
            self.rep = int(values[16])
            self.inc = int(values[17])
            self.target_description = values[18]

    def __str__(self):
        return "\t".join([self.target_name,
                          self.target_accession,
                          self.query_name,
                          self.query_accession,
                          str(self.full_e_value),
                          str(self.full_score),
                          str(self.full_bias),
                          str(self.best_e_value),
                          str(self.best_score),
                          str(self.best_bias),
                          str(self.exp),
                          str(self.reg),
                          str(self.clu),
                          str(self.ov),
                          str(self.env),
                          str(self.dom),
                          str(self.rep),
                          str(self.inc),
                          self.target_description
                          ]
                         )


class HmmerHitDOM():
    """Encapsulate a HMMER hit given in domtblout format."""
    def __init__(self, values):
        if len(values) == 23:
            self.target_name = values[0]
            self.target_accession = values[1]
            self.target_length = int(values[2])
            self.query_name = values[3]

            self.query_accession = values[4]
            if self.query_accession == '-':
                self.query_accession = self.query_name

            self.query_length = int(values[5])
            self.full_e_value = float(values[6])
            self.full_score = float(values[7])
            self.full_bias = float(values[8])
            self.dom = int(values[9])
            self.ndom = int(values[10])
            self.c_evalue = float(values[11])
            self.i_evalue = float(values[12])
            self.dom_score = float(values[13])
            self.dom_bias = float(values[14])
            self.hmm_from = int(values[15])
            self.hmm_to = int(values[16])
            self.ali_from = int(values[17])
            self.ali_to = int(values[18])
            self.env_from = int(values[19])
            self.env_to = int(values[20])
            self.acc = float(values[21])
            self.target_description = values[22]

    def __str__(self):
        return "\t".join([self.target_name,
                          self.target_accession,
                          str(self.target_length),
                          self.query_name,
                          self.query_accession,
                          str(self.query_length),
                          str(self.full_e_value),
                          str(self.full_score),
                          str(self.full_bias),
                          str(self.dom),
                          str(self.ndom),
                          str(self.c_evalue),
                          str(self.i_evalue),
                          str(self.dom_score),
                          str(self.dom_bias),
                          str(self.hmm_from),
                          str(self.hmm_to),
                          str(self.ali_from),
                          str(self.ali_to),
                          str(self.env_from),
                          str(self.env_to),
                          str(self.acc),
                          self.target_description]
                         )


class HmmModel(object):
    """Store HMM parameters."""
    def __init__(self, keys):
        setattr(self, 'ga', None)
        setattr(self, 'tc', None)
        setattr(self, 'nc', None)

        if 'acc' not in keys:
            setattr(self, 'acc', keys['name'])

        for key, value in list(keys.items()):
            setattr(self, key, value)


class HmmModelParser(object):
    """Parse HMM file."""
    def __init__(self, hmmFile):
        self.hmmFile = open(hmmFile)

    def models(self):
        """Parse all models from HMM file."""
        models = {}
        for model in self.simple_parse():
            key = model.acc
            models[key] = model

        return models

    def simple_parse(self):
        """Parse simplified description of single model from HMM file."""
        headerKeys = dict()
        for line in self.hmmFile:
            if line.startswith("HMMER"):
                headerKeys['format'] = line.rstrip()
            elif line.startswith("HMM"):
                # beginning of the hmm; iterate through till end of model
                for line in self.hmmFile:
                    if line.startswith("//"):
                        yield HmmModel(headerKeys)
                        break
            else:
                # header sections
                fields = line.rstrip().split(None, 1)
                if len(fields) != 2:
                    raise HmmModelError
                else:
                    # transform some data based on some of the header tags
                    if fields[0] == 'ACC' or fields[0] == 'NAME' or fields[0] == 'DESC':
                        headerKeys[fields[0].lower()] = fields[1]
                    elif fields[0] == "LENG":
                        headerKeys[fields[0].lower()] = int(fields[1])
                    elif fields[0] == "GA" or fields[0] == "TC" or fields[0] == "NC":
                        params = fields[1].split()
                        if len(params) != 2:
                            raise HmmModelError
                        headerKeys[fields[0].lower()] = (float(params[0].replace(';', '')), float(params[1].replace(';', '')))
                    else:
                        pass

    def parse(self):
        """Parse full description of single model from HMM file."""
        fields = []
        headerKeys = dict()
        for line in self.hmmFile:
            if line.startswith("HMMER"):
                headerKeys['format'] = line.rstrip()
            elif line.startswith("HMM"):
                # beginning of the hmm; iterate through till end of model
                for line in self.hmmFile:
                    if line.startswith("//"):
                        yield HmmModel(headerKeys)
                        break
            else:
                # header sections
                fields = line.rstrip().split(None, 1)
                if len(fields) != 2:
                    raise HmmModelError
                else:
                    # transform some data based on some of the header tags
                    if fields[0] == "LENG" or fields[0] == "NSEQ" or fields[0] == "CKSUM":
                        headerKeys[fields[0].lower()] = int(fields[1])
                    elif fields[0] == "RF" or fields[0] == "CS" or fields[0] == "MAP":
                        if fields[1].lower() == "no":
                            headerKeys[fields[0].lower()] = False
                        else:
                            headerKeys[fields[0].lower()] = True
                    elif fields[0] == "EFFN":
                        headerKeys[fields[0].lower()] = float(fields[1])
                    elif fields[0] == "GA" or fields[0] == "TC" or fields[0] == "NC":
                        params = fields[1].split()
                        if len(params) != 2:
                            raise HmmModelError
                        headerKeys[fields[0].lower()] = (float(params[0].replace(';', '')), float(params[1].replace(';', '')))
                    elif fields[0] == "STATS":
                        params = fields[1].split()
                        if params[0] != "LOCAL":
                            raise HmmModelError
                        if params[1] == "MSV" or params[1] == "VITERBI" or params[1] == "FORWARD":
                            headerKeys[(fields[0] + "_" + params[0] + "_" + params[1]).lower()] = (float(params[2]), float(params[3]))
                        else:
                            print(("'" + params[1] + "'"))
                            raise HmmModelError
                    else:
                        headerKeys[fields[0].lower()] = fields[1]
