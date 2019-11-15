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
import logging
import ntpath
import operator
from collections import defaultdict

import biolib.seq_io as seq_io
import biolib.seq_tk as seq_tk
from biolib.parallel import Parallel

from numpy import mean, std, array, dot, sqrt
from numpy.linalg import inv
from scipy.stats import f


class LgtDinucleotide(object):
    """Calculate dinucleotide usage of genes within genomes.

    The dinucleotides are formed from the 3rd and succeeding
    1st position nucleotides. This has been suggested for
    identifying lateral gene transfer as this dinucleotide
    patten is minimally restricted by amino acid preference
    and codon usage:

    Hooper SD, Berg OG. 2002. Detection of genes with atypical
        nucleotide sequence in microbial genomes. J. Mol. Evol.
        54:365-75
    """

    def __init__(self, cpus=1):
        """Initialization.

        Parameters
        ----------
        cpus : int
            Number of cpus to use.
        """
        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

    def _hotelling_statistic(self, gene_di_usage,
                                    gene_n3,
                                    gene_n1,
                                    genome_di_usage,
                                    genome_n3,
                                    genome_n1):
        """Calculate Hotelling's T^2-statistic.

        Parameters
        ----------
        gene_di_usage : dict[gene_id][dinucleotide] -> count
            Occurrence of dinucleotides in each gene.
        gene_n3 : dict[gene_id][dinucleotide] -> count
            Occurrence of nucleotides in the 3rd codon position in each gene.
        gene_n1 : dict[gene_id][dinucleotide] -> count
            Occurrence of nucleotides in the 1st codon position in each gene.
        genome_di_usage : dict[gene_id][dinucleotide] -> count
            Occurrence of dinucleotides across the genome.
        genome_n3 : dict[gene_id][dinucleotide] -> count
            Occurrence of nucleotides in the 3rd codon position across the genome.
        genome_n1 : dict[gene_id][dinucleotide] -> count
            Occurrence of nucleotides in the 1st codon position across the genome.

        Returns
        -------
        dict : d[gene_id] ->  Hotelling's T^2-statistic
             Hotelling's T^2-statistic for each gene.
        """

        # calculate mean and std. dev. for each gene
        nucleotides = ['A', 'C', 'G', 'T']
        gene_di_mean = defaultdict(lambda: defaultdict(float))
        gene_di_std = defaultdict(lambda: defaultdict(float))
        for gene_id, dinucleotides in gene_di_usage.items():
            n1 = gene_n1[gene_id]
            n3 = gene_n3[gene_id]
            N = sum(dinucleotides.values())
            for a in nucleotides:
                for b in nucleotides:
                    m = float(n3[a] * n1[b]) / N
                    gene_di_mean[gene_id][a + b] = m
                    var = m * (1.0 - min(n3[a], n1[b]) / N) * ((N - max(n3[a], n1[b])) / (N - 1.0))
                    gene_di_std[gene_id][a + b] = sqrt(var)

        # calculate mean and std. dev. for the entire genome
        genome_di_mean = defaultdict(float)
        genome_di_std = defaultdict(float)
        N = sum(genome_di_usage.values())
        for a in nucleotides:
            for b in nucleotides:
                m = float(genome_n3[a] * genome_n1[b]) / N
                genome_di_mean[a + b] = m
                var = m * (1.0 - (min(genome_n3[a], genome_n1[b]) / N)) * ((N - max(genome_n3[a], genome_n1[b])) / (N - 1.0))
                genome_di_std[a + b] = sqrt(var)

        # calculate dinculotide bias for each gene and the entire genome
        gene_di_bias = defaultdict(lambda: defaultdict(float))
        for gene_id in gene_di_usage:
            for di in genome_di_usage:
                if gene_di_std[gene_id][di] != 0:
                    gene_di_bias[gene_id][di] = float(gene_di_usage[gene_id].get(di, 0) - gene_di_mean[gene_id][di]) / gene_di_std[gene_id][di]
                else:
                    gene_di_bias[gene_id][di] = 0

        genome_di_bias = defaultdict(float)
        for di in genome_di_usage:
            genome_di_bias[di] = float(genome_di_usage.get(di, 0) - genome_di_mean[di]) / genome_di_std[di]

        # calculate covariance matrix
        cov_matrix = []
        for di_j in list(genome_di_bias.keys())[0:-1]:
            row = []
            for di_k in list(genome_di_bias.keys())[0:-1]:
                s = 0
                for dinucleotides in list(gene_di_bias.values()):
                    s += (dinucleotides[di_j] - genome_di_bias[di_j]) * (dinucleotides[di_k] - genome_di_bias[di_k])
                s /= (len(gene_di_bias) - 1)

                row.append(s)
            cov_matrix.append(row)
        inv_cov_matrix = inv(array(cov_matrix))

        # calculate Hotelling's T^2-statistic for each gene
        T2 = {}
        for gene_id, dinucleotides in gene_di_bias.items():
            x = []
            for di in list(genome_di_bias.keys())[0:-1]:
                x.append(dinucleotides[di] - genome_di_bias[di])
            x = array(x)

            T2[gene_id] = dot(dot(x, inv_cov_matrix), x)

        return T2

    def _manhattan(self, gene_di_usage, genome_di_usage):
        """Calculate Manhattan distance.

        Parameters
        ----------
        gene_di_usage : dict[gene_id][dinucleotide] -> count
            Occurrence of dinucleotides in each gene.
        genome_di_usage : dict[gene_id][dinucleotide] -> count
            Occurrence of dinucleotides across the genome.

        Returns
        -------
        dict : d[gene_id] ->  [distance, std. dev. from mean]
             Manhattan distance and standard deviation from the mean for each gene.
        """

        dist = {}
        genome_sum_di = sum(genome_di_usage.values())
        for gene_id, dinucleotides in gene_di_usage.items():
            d = 0
            gene_sum_di = sum(dinucleotides.values())
            for di in genome_di_usage:
                d += abs(dinucleotides.get(di, 0) * 100.0 / gene_sum_di - genome_di_usage.get(di, 0) * 100.0 / genome_sum_di)
            dist[gene_id] = [d]

        # model all distances as a normal distribution
        m = mean(list(dist.values()))
        s = std(list(dist.values()))

        # calculate standard deviations from the mean
        for gene_id, d in dist.items():
            dist[gene_id].append((d - m) / s)

        return dist

    def dinucleotide_usage(self, seqs, genome_id):
        """Calculate dinucleotide (n3:n1) usage statistics for genome.

        Parameters
        ----------
        seqs : dict[seq_id] -> seq
            Sequences indexed by sequence id.
        genome_id : str
            Unique id of genome used to create output file.
        """

        # calculate dinucleotide usage for each gene and the genome as a while
        gene_di_usage = defaultdict(lambda: defaultdict(int))
        gene_n1 = defaultdict(lambda: defaultdict(int))
        gene_n3 = defaultdict(lambda: defaultdict(int))

        genome_di_usage = defaultdict(int)
        genome_n1 = defaultdict(int)
        genome_n3 = defaultdict(int)
        gc = {}
        for gene_id, seq in seqs.items():
            gc[gene_id] = seq_tk.gc(seq)

            for i in range(2, len(seq) - 2, 3):
                dinucleotide = seq[i:i + 2].upper()
                if 'N' not in dinucleotide:
                    gene_di_usage[gene_id][dinucleotide] += 1
                    gene_n3[gene_id][dinucleotide[0]] += 1
                    gene_n1[gene_id][dinucleotide[1]] += 1

                    genome_di_usage[dinucleotide] += 1
                    genome_n3[dinucleotide[0]] += 1
                    genome_n1[dinucleotide[1]] += 1

        # calculate Manhattan distance for each gene
        manhattan_dist = self._manhattan(gene_di_usage, genome_di_usage)

        # identify deviant genes under Hotelling T-squared statistic
        t2_stats = self._hotelling_statistic(gene_di_usage, gene_n3, gene_n1,
                                                genome_di_usage, genome_n3, genome_n1)

        # correction from T-squared distribution to F-distribution approximation
        # http://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
        num_genes = len(gene_di_usage)
        f_dist_correction = float(num_genes - 15 + 1) / (15 * num_genes)
        deviant_threshold = f.isf(self.critical_value, 15, num_genes - 15 + 1) / f_dist_correction

        # report dinucleotide usage of each gene
        di_set_sorted = sorted(genome_di_usage.keys())
        gene_ids_sorted = sorted(list(t2_stats.items()), key=operator.itemgetter(1), reverse=True)

        output_file = os.path.join(self.output_dir, genome_id + '.di_usage.tsv')
        fout = open(output_file, 'w')

        fout.write('Gene Id\tGC\tLength (bp)\t# dinucleotides\tHotelling T-squared statistic\tDeviant\tManhattan distance\tDeviations from mean')
        for di in di_set_sorted:
            fout.write('\t' + di)
        fout.write('\n')

        genome_gc = seq_tk.gc(''.join(list(seqs.values())))
        genome_sum_di = sum(genome_di_usage.values())
        fout.write('%s\t%.2f\t%d\t%d' % ('<complete genome>', genome_gc * 100.0, sum([len(x) for x in list(seqs.values())]), genome_sum_di))
        fout.write('\t%s\t%s\t%.1f\t%.1f' % ('na', 'na', 0, 0))
        for di in di_set_sorted:
            fout.write('\t%.2f' % (genome_di_usage.get(di, 0) * 100.0 / genome_sum_di))
        fout.write('\n')

        for gene_id, t2_stat in gene_ids_sorted:
            dinucleotides = gene_di_usage[gene_id]
            sum_di = sum(dinucleotides.values())
            fout.write('%s\t%.2f\t%d\t%d' % (gene_id, gc[gene_id] * 100, len(seqs[gene_id]), sum_di))
            fout.write('\t%.2f\t%s' % (t2_stat, 'yes' if t2_stat > deviant_threshold else 'no'))
            fout.write('\t%.2f\t%.2f' % (manhattan_dist[gene_id][0], manhattan_dist[gene_id][1]))

            for di in di_set_sorted:
                fout.write('\t%.2f' % (dinucleotides.get(di, 0) * 100.0 / sum_di))
            fout.write('\n')
        fout.close()

    def _producer(self, gene_file):
        """Calculates dinucleotide usage statistics of a genome.

        Parameters
        ----------
        gene_file : str
            Fasta file containing amino acid sequences.
        """

        genome_id = ntpath.basename(gene_file)
        genome_id = genome_id.replace('.genes.fna', '')
        genome_id = os.path.splitext(genome_id)[0]

        seqs = seq_io.read_fasta(gene_file)
        self.dinucleotide_usage(seqs, genome_id)

        return True

    def _progress(self, processed_items, total_items):
        """Report progress of consumer processes.

        Parameters
        ----------
        processed_items : int
            Number of genomes processed.
        total_items : int
            Total number of genomes to process.

        Returns
        -------
        str
            String indicating progress of data processing.
        """

        return '    Finished processing %d of %d (%.2f%%) genomes.' % (processed_items, total_items, float(processed_items) * 100 / total_items)

    def run(self, gene_files, critical_value, output_dir):
        """Calculate dinucleotide usage over a set of genomes.

        Parameters
        ----------
        gene_files : list
            Fasta files containing called genes in nucleotide space.
        critical_value : float
            Critical value used to define a deviant gene (i.e., potential LGT event).
        output_dir : str
            Directory to store results.
        """

        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.critical_value = critical_value

        self.logger.info('Calculating dinucleotide usage for each genome.')
        
        progress_func = self._progress
        if self.logger.is_silent:
            progress_func = None

        parallel = Parallel(self.cpus)
        parallel.run(self._producer, None, gene_files, progress_func)
