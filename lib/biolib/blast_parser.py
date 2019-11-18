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

import gzip
from collections import namedtuple


class BlastParser():
    """Parses output files produced with Blast."""

    def __init__(self):
        """Initialization."""
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

    def read_hit(self, table):
        """Generator function to read hits from a blast output table.

        The table should be in blast format 6. This is
        also the format used by Diamond.

        Parameters
        ----------
        table : str
            Name of table to read.

        Yields
        ------
        namedtuple
            Information about blast hit.
        """

        if table.endswith('.gz'):
            open_file = gzip.open
        else:
            open_file = open

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

    def identify_homologs(self,
                          blast_table,
                          evalue_threshold,
                          per_identity_threshold,
                          per_aln_len_threshold,
                          seq_lens):
        """Identify homologs among blast hits.

        Identifies hits satisfying the criteria required for a
        gene to be considered a homolog. The table should be in
        blast format 6.

        Parameters
        ----------
        blast_table : str
            File containing blast hits in the custom tabular format produced by BlastRunner.
        evalue_threshold : float
            E-value threshold used to define homologous gene.
        per_identity_threshold : float
            Percent identity threshold used to define a homologous gene.
        per_aln_len_threshold : float
            Alignment length threshold used to define a homologous gene.
        seq_lens : dict
            Length of sequences indexed by their unique id.

        Returns
        -------
        set
            Identifiers for homologous genes.
        """

        homologs = set()
        for hit in self.read_hit(blast_table):

            if hit.evalue <= evalue_threshold and hit.perc_identity >= per_identity_threshold:
                query_len = seq_lens[hit.query_id]
                per_aln_len = hit.aln_length * 100.0 / query_len

                if per_aln_len >= per_aln_len_threshold:
                    homologs.add(hit.subject_id)

        return homologs
