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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import tempfile
import logging

from biolib.external.execute import check_on_path


class Diamond(object):
    """Wrapper for running diamond."""

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        check_on_path('diamond')

        self.cpus = cpus
        
        self.output_fmt = { 'diamond_daa': '100',
                            'standard': '6',
                            'custom': '6 qseqid qlen sseqid stitle slen length pident evalue bitscore'}

    def make_database(self, prot_file, db_file):
        """Make diamond database.

        Parameters
        ----------
        prot_file : str
            Fasta file with protein sequences.
        db_file : str
            Desired name of Diamond database.
        """

        cmd = 'diamond makedb --quiet -p %d --in %s -d %s' % (self.cpus, prot_file, db_file)
        os.system(cmd)

    def blastp(self, prot_file, 
                    db_file, 
                    evalue, 
                    per_identity, 
                    per_aln_len, 
                    max_target_seqs,
                    sensitive,
                    output_file, 
                    output_fmt='standard', 
                    tmp_dir=None, 
                    chunk_size=None, 
                    block_size=None):
        """Apply diamond blastp to a set of protein sequences.

        Parameters
        ----------
        prot_file : str
            Fasta file with protein sequences.
        db_file : str
            Diamond database of protein sequences.
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast [0, 100].
        per_aln_len : float
            Percent query coverage threshold for reporting hits [0, 100].
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        sensitive : boolean
            Run DIAMOND in sensitive mode.
        output_file : str
            Desired name of output file.
        output_fmt : str
            Desired output format (diamond_daa, standard, custom)
        tmp_dir : str
            Directory to store temporary files.
        chunk_size : int
            Number of chunks for index processing.
        block_size : int
            Sequence block size in billions of letters.
        """
        
        assert(output_fmt in self.output_fmt)

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        args = ''
        if tmp_dir:
             args += '-t %s' % tmp_dir

        if chunk_size:
             args += ' -c %d' % chunk_size
             
        if block_size:
            args += ' -b %d' % block_size
            
        if sensitive:
            args += ' --sensitive'

        cmd = "diamond blastp --quiet -p %d -q %s -d %s -e %g --id %f --query-cover %f -k %d -o %s -f %s %s" % (self.cpus,
                                                                                                                    prot_file,
                                                                                                                    db_file,
                                                                                                                    evalue,
                                                                                                                    per_identity,
                                                                                                                    per_aln_len,
                                                                                                                    max_target_seqs,
                                                                                                                    output_file,
                                                                                                                    self.output_fmt[output_fmt],
                                                                                                                    args)

        os.system(cmd)

    def blastx(self, nt_file, db_file, evalue, per_identity, max_target_seqs, diamond_daa_file):
        """Apply diamond blastx to a set of nucleotide sequences.

        Parameters
        ----------
        nt_file : str
            Fasta file with nucleotide sequences.
        db_file : str
            Diamond database of protein sequences.
        evalue : float
            E-value threshold used by blast.
        per_identity : float
            Percent identity threshold used by blast [0, 100].
        max_target_seqs : int
            Maximum number of hits to report per sequence.
        diamond_daa_file : str
            Desired name of Diamond data file.
        """

        if db_file.endswith('.dmnd'):
            db_file = db_file[0:db_file.rfind('.dmnd')]

        cmd = 'diamond blastx --quiet -p %d -t %s -q %s -d %s -e %f --id %f -k %d -a %s' % (self.cpus,
                                                                                                tempfile.gettempdir(),
                                                                                                nt_file,
                                                                                                db_file,
                                                                                                evalue,
                                                                                                per_identity,
                                                                                                max_target_seqs,
                                                                                                diamond_daa_file)

        os.system(cmd)

    def view(self, diamond_daa_file, output_table, compress=False):
        """Generate flat file from diamond DAA file.

        Parameters
        ----------
        diamond_daa_file : str
            Diamond DAA file.
        output_table : str
            Diamond database of protein sequeces.
        compress : boolean
            Flag indicating if output table should be compressed.
        """

        cmd = 'diamond view --quiet -p %d -a %s -o %s --compress %d' % (self.cpus,
                                                                            diamond_daa_file,
                                                                            output_table,
                                                                            compress)
        os.system(cmd)
