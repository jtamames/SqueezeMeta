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

from comparem.similarity_search import SimilaritySearch
from comparem.aai_calculator import AAICalculator
from comparem.classify import Classify
from comparem.codon_usage import CodonUsage
from comparem.amino_acid_usage import AminoAcidUsage
from comparem.kmer_usage import KmerUsage
from comparem.lgt_dinucleotide import LgtDinucleotide
from comparem.lgt_codon import LgtCodon
from comparem.hierarchical_clustering import HierarchicalCluster
#from comparem.PCoA import PCoA
from comparem.plots.heatmap import Heatmap

import biolib.seq_io as seq_io
from biolib.external.prodigal import Prodigal
from biolib.common import (remove_extension,
                             make_sure_path_exists,
                             check_dir_exists,
                             check_file_exists,
                             concatenate_files)
                             
from scipy.spatial.distance import (pdist as scipy_pdist,
                                    squareform as scipy_squareform)


class OptionsParser():
    def __init__(self):
        self.logger = logging.getLogger('timestamp')

    def _input_files(self, input, file_ext):
        """Identify genomes files.

        Parameters
        input : str
            File or directory specifying input files to process.
        file_ext : str
            Extension of files to process.

        Returns
        -------
        list
            Name of files to process.
        """
        
        input_files = []
        if os.path.isfile(input):
            for line in  open(input):
                input_file = line.strip().split('\t')[0]
                if not os.path.exists(input_file):
                    self.logger.error('Specified input file does not exist: %s' % input_file)
                    sys.exit()
                input_files.append(input_file)
                
            if not input_files:
                self.logger.warning('No genomes found in file: %s. Check that the file has the correct format.' % input)
                sys.exit()
        elif os.path.isdir(input):
            for f in os.listdir(input):
                if f.endswith(file_ext):
                    input_files.append(os.path.join(input, f))
                    
            if not input_files:
                self.logger.warning('No genomes found in directory: %s. Check the --file_ext flag used to identify genomes.' % input)
                sys.exit()
        else:
            self.logger.error('Specified input file or directory does not exists: %s' % input)
            sys.exit()

        return input_files

    def _write_usage_profile(self, genome_usage, feature_set, counts, output_file):
        """Write out occurrence of specified features for each genome.

        Parameters
        ----------
        genome_usage : d[genome_id][feature] -> count
            Occurrence of genomic feature in genome
        feature_set : iterable
            All genomic features.
        counts : boolean
            Write raw counts if True, else write frequencies.
        output_file : str
            File to produce.
        """

        sorted_feature_set = sorted(feature_set)

        fout = open(output_file, 'w')
        fout.write('Genome ID')
        for feature in sorted_feature_set:
            fout.write('\t' + feature)
        fout.write('\n')

        totals = defaultdict(int)
        for genome_id, features in genome_usage.items():
            for feature in sorted_feature_set:
                totals[genome_id] += features.get(feature, 0)

        for genome_id, features in genome_usage.items():
            fout.write(genome_id)

            for feature in sorted_feature_set:
                if counts:
                    fout.write('\t%d' % features.get(feature, 0))
                else:
                    fout.write('\t%f' % (features.get(feature, 0) * 100.0 / totals[genome_id]))
            fout.write('\n')

    def ani(self, options):
        """ANI command"""

        make_sure_path_exists(options.output_dir)

        genome_files = self._genome_files(options.genome_dir, options.file_ext)

        self.logger.info('Average nucleotide identity information written to: %s' % options.output_dir)

    def call_genes(self, options):
        """Call genes command"""

        make_sure_path_exists(options.output_dir)
        
        genome_files = self._input_files(options.input_genomes, options.file_ext)

        prodigal = Prodigal(options.cpus, not options.silent)
        summary_stats = prodigal.run(genome_files, 
                                        options.output_dir, 
                                        called_genes=False, 
                                        translation_table=options.force_table, 
                                        meta=False,
                                        closed_ends=True)

        # write gene calling summary
        fout = open(os.path.join(options.output_dir, 'call_genes.summary.tsv'), 'w')
        fout.write('Genome Id\tSelected translation table\tTable 4 coding density\tTable 11 coding density\n')
        for genome_id, stats in summary_stats.items():
            fout.write('%s\t%d\t%.2f%%\t%.2f%%\n' % (genome_id,
                                                     stats.best_translation_table,
                                                     stats.coding_density_4,
                                                     stats.coding_density_11))
        fout.close()

        self.logger.info('Identified genes written to: %s' % options.output_dir)

    def similarity(self, options):
        """Perform sequence similarity search between genes"""

        make_sure_path_exists(options.output_dir)
        
        query_gene_files = self._input_files(options.query_proteins, options.file_ext)
        target_gene_files = self._input_files(options.target_proteins, options.file_ext)
        
        ss = SimilaritySearch(options.cpus)
        ss.run(query_gene_files, 
                target_gene_files,
                options.evalue, 
                options.per_identity, 
                options.per_aln_len,
                True,
                options.tmp_dir,
                options.blastp,
                options.sensitive,
                options.keep_headers,
                options.output_dir)

        self.logger.info('Sequence similarity results written to: %s' % options.output_dir)
        
    def aai(self, options):
        """AAI command"""
        check_file_exists(options.sorted_hit_table)
        make_sure_path_exists(options.output_dir)

        aai_calculator = AAICalculator(options.cpus)
        aai_output_file, rbh_output_file = aai_calculator.run(options.query_gene_file,
                                                                None,
                                                                options.sorted_hit_table,
                                                                options.evalue,
                                                                options.per_identity,
                                                                options.per_aln_len,
                                                                options.keep_rbhs,
                                                                options.output_dir)

        if rbh_output_file:
            self.logger.info('Identified reciprocal best hits written to: %s' % rbh_output_file)
            
        self.logger.info('AAI between genomes written to: %s' % aai_output_file)
        
    def classify(self, options):
        """Classify genomes based on AAI values."""
        check_file_exists(options.sorted_hit_table)
        make_sure_path_exists(options.output_dir)
        
        classify = Classify(options.cpus)
        results_file = classify.run(options.query_gene_file,
                                        options.target_gene_file,
                                        options.sorted_hit_table,
                                        options.evalue,
                                        options.per_identity,
                                        options.per_aln_len,
                                        options.num_top_targets,
                                        options.taxonomy_file,
                                        options.keep_rbhs,
                                        options.output_dir)
        
        self.logger.info('Classification results written to: %s' % results_file)

    def aa_usage(self, options):
        """Amino acid usage command"""

        gene_files = self._input_files(options.protein_gene_files, options.file_ext)

        # calculate amino acid usage
        amino_acid_usage = AminoAcidUsage(options.cpus)
        genome_aa_usage, aa_set = amino_acid_usage.run(gene_files)

        # write out results
        self._write_usage_profile(genome_aa_usage, aa_set, options.counts, options.output_file)

        self.logger.info('Amino acid usage written to: %s' % options.output_file)

    def codon_usage(self, options):
        """Codon usage command"""

        gene_files = self._input_files(options.nucleotide_gene_files, options.file_ext)

        # calculate amino acid usage
        codon_usage = CodonUsage(options.cpus, options.keep_ambiguous)
        genome_codon_usage, codon_set, _mean_length = codon_usage.run(gene_files)

        # write out results
        self._write_usage_profile(genome_codon_usage, codon_set, options.counts, options.output_file)

        self.logger.info('Codon usage written to: %s' % options.output_file)

    def stop_usage(self, options):
        """Stop codon usage command"""

        gene_files = self._input_files(options.nucleotide_gene_files, options.file_ext)

        # calculate amino acid usage
        codon_usage = CodonUsage(options.cpus, keep_ambiguous=False, stop_codon_only=True)
        genome_codon_usage, codon_set, mean_gene_length = codon_usage.run(gene_files)

        # write out results
        if not options.mean_gene_length:
            self._write_usage_profile(genome_codon_usage, codon_set, options.counts, options.output_file)
        else:
            fout = open(options.output_file, 'w')
            for codon in codon_set:
                fout.write('\t' + codon)
                if mean_gene_length:
                    fout.write('\t' + codon + ': avg. seq. length')
            fout.write('\n')

            for genome_id, codons in genome_codon_usage.items():
                fout.write(genome_id)

                for codon in codon_set:
                    fout.write('\t%d' % codons.get(codon, 0))

                    if mean_gene_length:
                        mean_len = mean_gene_length[genome_id].get(codon, None)
                        if mean_len:
                            fout.write('\t%.1f' % mean_len)
                        else:
                            fout.write('\tna')
                fout.write('\n')

        self.logger.info('Stop codon usage written to: %s' % options.output_file)

    def kmer_usage(self, options):
        """Kmer usage command"""

        if options.k > 10 or options.k <= 0:
            self.logger.warning('CompareM only support kmers with k <= 10.')
            sys.exit(0)

        genome_files = self._input_files(options.genome_files, options.file_ext)

        # calculate amino acid usage
        kmer_usage = KmerUsage(options.k, options.cpus)
        genome_kmer_usage, kmer_set = kmer_usage.run(genome_files)

        # write out results
        self.logger.info('Writing kmer profiles to file (be patient!).')
        self._write_usage_profile(genome_kmer_usage, kmer_set, options.counts, options.output_file)

        self.logger.info('Kmer usage written to: %s' % options.output_file)

    def lgt_di(self, options):
        """LGT dinucleotide usage command"""

        gene_files = self._input_files(options.nucleotide_gene_files, options.file_ext)

        lgt_dinucleotide = LgtDinucleotide(options.cpus)
        lgt_dinucleotide.run(gene_files, options.crit_value, options.output_dir)

        self.logger.info('Dinucleotide usage written to directory: %s' % options.output_dir)

    def lgt_codon(self, options):
        """LGT dinucleotide usage command"""

        gene_files = self._input_files(options.nucleotide_gene_files, options.file_ext)

        lgt_codon = LgtCodon(options.cpus)
        lgt_codon.run(gene_files, options.output_dir)

        self.logger.info('Codon usage written to directory: %s' % options.output_dir)
        
    def diss(self, options):
        """Calculate dissimilarity between usage profiles."""
        
        check_file_exists(options.profile_file)
        
        genome_ids = []
        profiles = []
        with open(options.profile_file) as f:
            f.readline() # burn header
            
            for line in f:
                line_split = line.rstrip().split('\t')
                genome_id = line_split[0]
                profile = [float(v) for v in line_split[1:]]
                
                genome_ids.append(genome_id)
                profiles.append(profile)
                
        # calculate dissimilarity between genomes
        d = scipy_pdist(profiles, metric=options.metric)

        fout = open(options.output_file, 'w')
        if not options.full_matrix:
            # write out lower triangle from condense dissimilarity matrix,
            # in pairwise fashion
            fout.write('Genome A\tGenome B\tDissimilarity\n')
            condensed_idx = lambda i,j,n: n*j - j*(j+1)/2 + i - 1 - j
            for i in range(1, len(genome_ids)):
                for j in range(i):
                    fout.write('%s\t%s\t%f\n' % (genome_ids[i], genome_ids[j], d[condensed_idx(i, j, len(genome_ids))]))
        else:
            # write out full dissimilarity matrix
            ds = scipy_squareform(d)
            for genome_id in genome_ids:
                fout.write('\t' + genome_id)
            fout.write('\n')
            
            for i, genome_id in enumerate(genome_ids):
                fout.write(genome_id)
                for j in range(len(genome_ids)):
                    fout.write('\t%f' % ds[i,j])
                fout.write('\n')
        
        fout.close()
        
        self.logger.info('Dissimilarity values written to: %s' % options.output_file)

    def pcoa_plot(self, options):
        """PCoA command"""

        self.logger.info('Performing PCoA.')
        #pcoa = PCoA()
        #pcoa.plot(options.aai_summary_file)

    def hclust(self, options):
        """Hierarchical clustering command"""
        
        if options.similarity and not options.max_sim_value:
            self.logger.error("The 'max_sim_value' must be specified for similarity values.")
            sys.exit(-1)

        self.logger.info('Performing hierarchical clustering.')
        hclust = HierarchicalCluster()
        hclust.run(options.pairwise_value_file,
                    options.method, 
                    options.similarity,
                    options.max_sim_value,
                    options.name_col1,
                    options.name_col2,
                    options.value_col,
                    options.output_tree)
        
        self.logger.info('Hierarchical cluster tree written to: %s' % options.output_tree)

    def heatmap(self, options):
        """Heatmap command"""

        self.logger.info('Making heatmap.')
        heatmapper = Heatmap(options.aai_summary_file, options.output_file)
        heatmapper.plot(options.cluster, options.method, options.metric)

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        try:
            if options.file == "stdout":
                options.file = ''
        except:
            pass

        if(options.subparser_name == 'call_genes'):
            self.call_genes(options)
        elif(options.subparser_name == 'similarity'):
            self.similarity(options)
        elif(options.subparser_name == 'aai'):
            self.aai(options)
        elif(options.subparser_name == 'classify'):
            self.classify(options)
        elif(options.subparser_name == 'aai_wf'):
            root_dir = options.output_dir
            make_sure_path_exists(root_dir)
            
            if options.proteins:
                if options.file_ext == 'fna':
                    self.logger.warning("Changing file extension from 'fna' to 'faa' since 'proteins' flag was given.")
                    options.file_ext = 'faa'
                options.query_proteins = options.input_files
                options.target_proteins = options.input_files
            else:
                options.input_genomes = options.input_files
                options.output_dir = os.path.join(root_dir, 'genes')
                self.call_genes(options)
                options.query_proteins = os.path.join(root_dir, 'genes')
                options.target_proteins = os.path.join(root_dir, 'genes')
                options.file_ext = 'faa'

            options.output_dir = os.path.join(root_dir, 'similarity')
            self.similarity(options)
            
            options.query_gene_file = os.path.join(options.output_dir, 'query_genes.faa')
            options.sorted_hit_table = os.path.join(options.output_dir, 'hits_sorted.tsv')
            options.output_dir = os.path.join(root_dir, 'aai')
            self.aai(options)
        elif(options.subparser_name == 'classify_wf'):
            root_dir = options.output_dir
            make_sure_path_exists(root_dir)
            
            if options.query_files == options.target_files:
                self.logger.error("The 'query_files' and 'target_files' arguments must be different.")
                sys.exit()
            
            if options.proteins:
                if options.file_ext == 'fna':
                    self.logger.warning("Changing file extension from 'fna' to 'faa' since 'proteins' flag was given.")
                    options.file_ext = 'faa'
                options.query_proteins = options.query_files
                options.target_proteins = options.target_files
            else:
                options.input_genomes = options.query_files
                options.output_dir = os.path.join(root_dir, 'query_genes')
                self.call_genes(options)
                
                options.input_genomes = options.target_files
                options.output_dir = os.path.join(root_dir, 'target_genes')
                self.call_genes(options)
                
                options.query_proteins = os.path.join(root_dir, 'query_genes')
                options.target_proteins = os.path.join(root_dir, 'target_genes')
                options.file_ext = 'faa'

            options.output_dir = os.path.join(root_dir, 'similarity')
            self.similarity(options)
            
            options.query_gene_file = os.path.join(options.output_dir, 'query_genes.faa')
            options.target_gene_file = os.path.join(options.output_dir, 'target_genes.faa')
            options.sorted_hit_table = os.path.join(options.output_dir, 'hits_sorted.tsv')
            options.output_dir = os.path.join(root_dir, 'classify')
            self.classify(options)
        elif(options.subparser_name == 'aa_usage'):
            self.aa_usage(options)
        elif(options.subparser_name == 'codon_usage'):
            self.codon_usage(options)
        elif(options.subparser_name == 'kmer_usage'):
            self.kmer_usage(options)
        elif(options.subparser_name == 'stop_usage'):
            self.stop_usage(options)
        elif(options.subparser_name == 'lgt_di'):
            self.lgt_di(options)
        elif(options.subparser_name == 'lgt_codon'):
            self.lgt_codon(options)
        elif(options.subparser_name == 'diss'):
            self.diss(options)
        elif(options.subparser_name == 'hclust'):
            self.hclust(options)
        elif(options.subparser_name == 'pcoa_plot'):
            self.pcoa_plot(options)
        elif(options.subparser_name == 'heatmap'):
            self.heatmap(options)
        else:
            self.logger.error('  [Error] Unknown CompareM command: "' + options.subparser_name + '"\n')
            sys.exit()

        return 0
