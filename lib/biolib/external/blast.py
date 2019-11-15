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

__author__ = "Donovan Parks"
__copyright__ = "Copyright 2015"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import logging
from collections import namedtuple

from biolib.external.execute import check_dependencies

"""
To do:
 - this class and the diamond class should mirror each other
   to the extent possible
"""


class Blast():
    """Wrapper for running blast."""

    def __init__(self, cpus, silent=False):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """

        self.logger = logging.getLogger('timestamp')

        check_dependencies(['blastn', 'blastp', 'makeblastdb'])

        self.cpus = cpus
        self.silent = silent

        self.output_fmt = {'standard': '6',
                            'custom': '6 qseqid qlen sseqid stitle slen length pident evalue bitscore'}
        self.blastp_tasks = {'blastp', 'blastp-fast', 'blastp-short'}
        self.blastn_tasks = {'blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn'}

        self.BlastHit = namedtuple('BlastHit', """query_id
                                                subject_id
                                                perc_identity
                                                aln_length
                                                mismatch_count
                                                gap_open_count
                                                query_start
                                                query_end
                                                subject_start
                                                subject_end
                                                evalue
                                                bitscore""")

        self.BlastHitCustom = namedtuple('BlastHitCustom', """query_id
                                                                query_len
                                                                subject_id
                                                                subject_annotation
                                                                subject_len
                                                                alignment_len
                                                                perc_identity
                                                                evalue
                                                                bitscore""")

        self.BlastHitHomologs = namedtuple('BlastHitHomologs', """query_id
                                                                subject_id
                                                                subject_annotation
                                                                perc_identity
                                                                query_perc_aln_len
                                                                subject_perc_aln_len
                                                                evalue
                                                                bitscore""")

    def blastp(self, query_seqs, prot_db, output_file, evalue=1e-3, max_matches=500, output_fmt='standard', task='blastp'):
        """Apply blastp to query file.

        Finds homologs to query sequences using blastp homology search
        against a protein database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        prot_db : str
            File containing blastp formatted database.
        output_file : str
            Output file containing blastp results.
        evalue : float
            E-value threshold used to identify homologs.
        max_matches : int
            Maximum hits per query sequence.
        output_fmt : str
            Specified output format of blast table: standard or custom.
        """

        assert(output_fmt in list(self.output_fmt.keys()))
        assert(task in self.blastp_tasks)

        cmd = "blastp -num_threads %d" % self.cpus
        cmd += " -query %s -db %s -out %s -evalue %g" % (query_seqs, prot_db, output_file, evalue)
        cmd += " -max_target_seqs %d" % max_matches
        cmd += " -task %s" % task
        cmd += " -outfmt '%s'" % self.output_fmt[output_fmt]
        os.system(cmd)

    def blastn(self, query_seqs, nucl_db, output_file, evalue=1e-3, max_matches=500, output_fmt='standard', task='megablast'):
        """Apply blastn to query file.

        Finds homologs to query sequences using blastn homology search
        against a nucleotide database. Hit can be reported using  either
        the 'standard' table 6 format or the following 'custom' format:
            qseqid qlen sseqid slen length pident evalue bitscore


        Parameters
        ----------
        query_seqs : str
            File containing query sequences.
        nucl_db : str
            File containing blastn formatted database.
        output_file : str
            Output file containing blastn results.
        evalue : float
            E-value threshold used to identify homologs.
        max_matches : int
            Maximum hits per query sequence.
        output_fmt : str
            Specified output format of blast table: standard or custom.
        """

        assert(output_fmt in list(self.output_fmt.keys()))
        assert(task in self.blastn_tasks)

        cmd = "blastn -num_threads %d" % self.cpus
        cmd += " -query %s -db %s -out %s -evalue %g" % (query_seqs, nucl_db, output_file, evalue)
        cmd += " -max_target_seqs %d" % max_matches
        cmd += " -task %s" % task
        cmd += " -outfmt '%s'" % self.output_fmt[output_fmt]
        os.system(cmd)

    def create_blastn_db(self, fasta_file):
        """Create nucleotide database."""

        if self.silent:
            os.system('makeblastdb -dbtype nucl -in %s > /dev/null' % fasta_file)
        else:
            os.system('makeblastdb -dbtype nucl -in %s' % fasta_file)
        
    def create_blastp_db(self, fasta_file):
        """Create protein database."""

        if self.silent:
            os.system('makeblastdb -dbtype prot -in %s > /dev/null' % fasta_file)
        else:
            os.system('makeblastdb -dbtype prot -in %s' % fasta_file)

    def read_hit(self, table, table_fmt):
        """Generator function to read hits from a blast output table.

        Parameters
        ----------
        table : str
            Name of table to read.
        table_fmt : str
            Specified output format of blast table: standard or custom.

        Yields
        ------
        namedtuple
            Information about blast hit.
        """

        assert(table_fmt in self.output_fmt)

        if table.endswith('.gz'):
            open_file = gzip.open
        else:
            open_file = open

        if table_fmt == 'standard':
            for line in open_file(table):
                line_split = line.split('\t')
                hit = self.BlastHit(query_id=line_split[0],
                                subject_id=line_split[1],
                                perc_identity=float(line_split[2]),
                                aln_length=int(line_split[3]),
                                mismatch_count=int(line_split[4]),
                                gap_open_count=int(line_split[5]),
                                query_start=int(line_split[6]),
                                query_end=int(line_split[7]),
                                subject_start=int(line_split[8]),
                                subject_end=int(line_split[9]),
                                evalue=float(line_split[10]),
                                bitscore=float(line_split[11]))

                yield hit
        else:
            for line in open_file(table):
                line_split = line.split('\t')
                hit = self.BlastHitCustom(query_id=line_split[0],
                                            query_len=int(line_split[1]),
                                            subject_id=line_split[2],
                                            subject_annotation=line_split[3],
                                            subject_len=int(line_split[4]),
                                            alignment_len=int(line_split[5]),
                                            perc_identity=float(line_split[6]),
                                            evalue=float(line_split[7]),
                                            bitscore=float(line_split[8]))

                yield hit

    def identify_homologs(self,
                          custom_blast_table,
                          evalue_threshold,
                          perc_identity_threshold,
                          perc_aln_len_threshold):
        """Identify homologs among blast hits based on specified criteria.

        Parameters
        ----------
        custom_blast_table : str
            File containing blast hits in the custom tabular format.
        evalue_threshold : float
            E-value threshold used to define homologs.
        perc_identity_threshold : float
            Percent identity threshold used to define a homologs.
        perc_aln_len_threshold : float
            Alignment length threshold used to define a homologs.

        Returns
        -------
        dict : d[subject_id] -> BlastHitCustom named tuple
            Dictionary with information about blast hits to homologs.
        """

        homologs = {}
        for line in open(custom_blast_table):
            line_split = line.split('\t')

            query_id = line_split[0]
            query_len = int(line_split[1])
            subject_id = line_split[2]
            subject_title = line_split[3]
            subject_len = int(line_split[4])
            align_len = int(line_split[5])
            perc_identity = float(line_split[6])
            evalue = float(line_split[7])
            bitscore = float(line_split[8])
            
            if evalue <= evalue_threshold and perc_identity >= perc_identity_threshold:
                query_perc_aln_len = align_len * 100.0 / query_len
                subject_perc_aln_len = align_len * 100.0 / subject_len

                if query_perc_aln_len >= perc_aln_len_threshold and subject_perc_aln_len >= perc_aln_len_threshold:
                    prev_hit = homologs.get(subject_id, None)
                    if not prev_hit or bitscore > prev_hit.bitscore:
                        homologs[subject_id] = self.BlastHitHomologs(query_id=query_id,
                                                                    subject_id=subject_id,
                                                                    subject_annotation=subject_title,
                                                                    perc_identity=perc_identity,
                                                                    query_perc_aln_len=query_perc_aln_len,
                                                                    subject_perc_aln_len=subject_perc_aln_len,
                                                                    evalue=evalue,
                                                                    bitscore=bitscore)

        return homologs
