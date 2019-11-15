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
__copyright__ = 'Copyright 2016'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import subprocess
import logging
import tempfile
import shutil

import biolib.seq_io as seq_io
from biolib.common import concatenate_files, remove_extension, make_sure_path_exists
from biolib.parallel import Parallel
from biolib.external.diamond import Diamond
from biolib.external.blast import Blast


class SimilaritySearch(object):
    """Performs similarity search of gene sequences."""

    def __init__(self, cpus):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus
        
    def _prefix_gene_identifiers(self, gene_files, keep_headers, output_file):
        """Prefix all gene IDs with genome IDs: <genome_id>~<gene_id>.
        
        Parameters
        ----------
        gene_files : list of str
            Genes in fasta files to modify.
        keep_headers : boolean
            If True, indicates FASTA headers already have the format <genome_id>~<gene_id>.
        output_file : str
            Name of FASTA file to contain modified genes.
        """
        
        fout = open(output_file, 'w')
        for gf in gene_files:           
            genome_id = remove_extension(gf)
            if genome_id.endswith('_genes'):
                genome_id = genome_id[0:genome_id.rfind('_genes')]
                
            for seq_id, seq, annotation in seq_io.read_fasta_seq(gf, keep_annotation=True):
                if keep_headers:
                    fout.write('>' + seq_id  + ' ' + annotation + '\n')
                else:
                    fout.write('>' + genome_id + '~' + seq_id  + ' ' + annotation + '\n')
                fout.write(seq + '\n')
        fout.close()
        
    def _sort_hit_table(self, input_hit_table, output_hit_table):
        """Sort hit table.
        
        Parameters
        ----------
        input_hit_table : str
            Table to sort.
        output_hit_table : str
            Name of sorted output table. 
        """
        
        self.logger.info('Sorting table with hits (be patient!).')
        os.system("LC_ALL=C sed -i 's/~/\t/g' %s" % input_hit_table)
        os.system("LC_ALL=C sort --parallel=8 -o %s -k1,1 -k3,3 %s" % (input_hit_table, input_hit_table))
        os.system('mv %s %s' % (input_hit_table, output_hit_table))
     
    def _run_self_blastp(self, query_gene_file, 
                                evalue, 
                                per_identity, 
                                per_aln_len,
                                max_hits,
                                tmp_dir,
                                output_dir):
        """Perform similarity search of query genes against themselves.

        Parameters
        ----------
        query_gene_file : str
            File with all query sequences.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        max_hits : int
            Maximum number of hits to report per query sequences.
        tmp_dir : str
            Directory to store temporary files.
        output_dir : str
            Directory to store blast results.
        """
        
        # concatenate all gene files and create a single diamond database
        self.logger.info('Creating BLASTP database (be patient!).')
        
        blast = Blast(self.cpus, silent=True)
        blast.create_blastp_db(query_gene_file)
        
        # create temporary hits table
        if tmp_dir:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_hits_table.close()

        # blast all genes against the database
        self.logger.info('Performing sequence similarity search between all query genomes (be patient!).')
        hits_daa_file = os.path.join(output_dir, 'query_hits')
        blast.blastp(query_gene_file, query_gene_file, tmp_hits_table.name, evalue, max_hits, task='blastp-fast')
        
        # sort hit table
        hits_table_file = os.path.join(output_dir, 'hits_sorted.tsv')
        self._sort_hit_table(tmp_hits_table.name, hits_table_file)
                
    def _run_self_diamond(self, query_gene_file, 
                                evalue, 
                                per_identity, 
                                per_aln_len,
                                max_hits,
                                sensitive,
                                high_mem,
                                tmp_dir,
                                output_dir):
        """Perform similarity search of query genes against themselves.

        Parameters
        ----------
        query_gene_file : str
            File with all query sequences.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        max_hits : int
            Maximum number of hits to report per query sequences.
        tmp_dir : str
            Directory to store temporary files.
        output_dir : str
            Directory to store blast results.
        """
        
        self.logger.info('Creating DIAMOND database (be patient!).')
        
        diamond_db = os.path.join(output_dir, 'query_genes')
        diamond = Diamond(self.cpus)
        diamond.make_database(query_gene_file, diamond_db)
            
        # create flat hits table
        if tmp_dir:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            tmp_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_hits_table.close()

        # blast all genes against the database
        self.logger.info('Performing self similarity sequence between genomes (be patient!).')

        if high_mem:
            diamond.blastp(query_gene_file, 
                            diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_hits_table.name, 
                            'standard', 
                            tmp_dir, 
                            chunk_size=1, 
                            block_size=8)
        else:
            diamond.blastp(query_gene_file, 
                            diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_hits_table.name, 
                            'standard', 
                            tmp_dir)

        # sort hit table
        hits_table_file = os.path.join(output_dir, 'hits_sorted.tsv')
        self._sort_hit_table(tmp_hits_table.name, hits_table_file)
        
        
    def _run_reciprocal_diamond(self, query_gene_file,
                                        target_gene_file,
                                        evalue, 
                                        per_identity, 
                                        per_aln_len,
                                        max_hits,
                                        sensitive,
                                        high_mem,
                                        tmp_dir,
                                        output_dir):
        """Perform similarity search of query genes against target genes, and reciprocal hits.

        Parameters
        ----------
        query_gene_file : str
            File with all query proteins.
        target_gene_file : str
            File with all target proteins.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        max_hits : int
            Maximum number of hits to report per query sequences.
        tmp_dir : str
            Directory to store temporary files.
        output_dir : str
            Directory to store blast results.
        """
        
        self.logger.info('Creating DIAMOND database of query proteins (be patient!).')
        diamond = Diamond(self.cpus)
        query_diamond_db = os.path.join(output_dir, 'query_genes')
        diamond.make_database(query_gene_file, query_diamond_db)
        
        self.logger.info('Creating DIAMOND database of target proteins (be patient!).')
        target_diamond_db = os.path.join(output_dir, 'target_genes')
        diamond.make_database(target_gene_file, target_diamond_db)

        # blast query genes against target proteins
        self.logger.info('Performing similarity sequence between query and target proteins (be patient!).')
        
        if tmp_dir:
            tmp_query_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            tmp_query_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_query_hits_table.close()
        
        query_hits_daa_file = os.path.join(output_dir, 'query_hits')
        
        if high_mem:
            diamond.blastp(query_gene_file, 
                            target_diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_query_hits_table.name, 
                            'standard', 
                            tmp_dir, 
                            chunk_size=1, 
                            block_size=8)
        else:
            diamond.blastp(query_gene_file, 
                            target_diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_query_hits_table.name, 
                            'standard', 
                            tmp_dir)
                
        # get target genes hit by one or more query proteins
        self.logger.info('Creating file with target proteins with similarity to query proteins.')
        target_hit = set()
        for line in open(tmp_query_hits_table.name):
            line_split = line.split('\t')
            target_hit.add(line_split[1])

        target_genes_hits = os.path.join(output_dir, 'target_genes_hit.faa')
        fout = open(target_genes_hits, 'w')
        for seq_id, seq in seq_io.read_seq(target_gene_file):
            if seq_id in target_hit:
                fout.write('>' + seq_id + '\n')
                fout.write(seq + '\n')
        fout.close()
        
        self.logger.info('Identified %d target proteins to be used in reciprocal search.' % len(target_hit))
        
        # perform reciprocal blast
        self.logger.info('Performing reciprocal similarity sequence between target and query proteins (be patient!).')

        if tmp_dir:
            tmp_target_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', dir=tmp_dir, delete=False)
        else:
            tmp_target_hits_table = tempfile.NamedTemporaryFile(prefix='comparem_hits_', delete=False)
        tmp_target_hits_table.close()
        
        if high_mem:
            diamond.blastp(target_genes_hits, 
                            query_diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_target_hits_table.name, 
                            'standard', 
                            tmp_dir, 
                            chunk_size=1, 
                            block_size=8)
        else:
            diamond.blastp(target_genes_hits, 
                            query_diamond_db, 
                            evalue, 
                            per_identity, 
                            per_aln_len, 
                            max_hits,
                            sensitive,
                            tmp_target_hits_table.name, 
                            'standard', 
                            tmp_dir)
                
        # combine hit tables and sort
        os.system('cat %s >> %s' % (tmp_target_hits_table.name, tmp_query_hits_table.name))
        os.remove(tmp_target_hits_table.name)
        hits_table_file = os.path.join(output_dir, 'hits_sorted.tsv')
        self._sort_hit_table(tmp_query_hits_table.name, hits_table_file)

    def run(self, query_gene_files, 
                    target_gene_files,
                    evalue, 
                    per_identity, 
                    per_aln_len,
                    high_mem,
                    tmp_dir,
                    blastp,
                    sensitive,
                    keep_headers,
                    output_dir):
        """Perform similarity search of query genes against target genes.

        Parameters
        ----------
        query_gene_files : list of str
            Query genes in fasta files to process.
        target_gene_files : list of str
            Query genes in fasta files to process.
        evalue : float
            E-value threshold for reporting hits.
        per_identity : float
            Percent identity threshold for reporting hits.
        per_aln_len : float
            Percent query coverage threshold for reporting hits.
        tmp_dir : str
            Directory to store temporary files.
        blastp : boolean
            If True, blasp-fast is used instead of DIAMOND.
        sensitive : boolean
            If True, the sensitive mode of DIAMOND is used.
        keep_headers : boolean
            If True, indicates FASTA headers already have the format <genome_id>~<gene_id>.
        output_dir : str
            Directory to store blast results.
        """
           
        # modify gene ids to include genome ids in order to ensure
        # all gene identifiers are unique across the set of genomes
        self.logger.info('Appending genome identifiers to query genes.')
        query_gene_file = os.path.join(output_dir, 'query_genes.faa')
        modified_query_gene_files = self._prefix_gene_identifiers(query_gene_files, 
                                                                    keep_headers,
                                                                    query_gene_file)
        
 
            
        if query_gene_files == target_gene_files:
            target_gene_file = query_gene_file
        else:
            self.logger.info('Appending genome identifiers to target genes.')
            target_gene_file = os.path.join(output_dir, 'target_genes.faa')
            modified_target_gene_files = self._prefix_gene_identifiers(target_gene_files, 
                                                                        keep_headers, 
                                                                        target_gene_file)

        if blastp:
            if target_gene_file == query_gene_file:
                self._run_self_blastp(query_gene_file, 
                                        evalue, 
                                        per_identity, 
                                        per_aln_len,
                                        len(target_gene_files) * 10,
                                        tmp_dir,
                                        output_dir)
            else:
                self.logger.info('NOT YET IMPLEMENTED!')
                sys.exit()
                self._run_reciprocal_blastp(query_gene_file, 
                                            target_gene_file, 
                                            evalue, 
                                            per_identity, 
                                            per_aln_len,
                                            len(target_gene_files) * 10,
                                            high_mem,
                                            tmp_dir,
                                            output_dir)
        else:
            if target_gene_file == query_gene_file:
                self._run_self_diamond(query_gene_file, 
                                        evalue, 
                                        per_identity, 
                                        per_aln_len,
                                        len(target_gene_files) * 10,
                                        sensitive,
                                        high_mem,
                                        tmp_dir,
                                        output_dir)

            
            else:
                self._run_reciprocal_diamond(query_gene_file, 
                                                target_gene_file, 
                                                evalue, 
                                                per_identity, 
                                                per_aln_len,
                                                len(target_gene_files) * 10,
                                                sensitive,
                                                high_mem,
                                                tmp_dir,
                                                output_dir)
