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
__copyright__ = "Copyright 2017"
__credits__ = ["Donovan Parks"]
__license__ = "GPL3"
__maintainer__ = "Donovan Parks"
__email__ = "donovan.parks@gmail.com"
__status__ = "Development"

import os
import logging
import time
import shutil
import tempfile
import ntpath

import biolib.seq_io as seq_io
from biolib.parallel import Parallel
from biolib.bootstrap import bootstrap_alignment, bootstrap_support
from biolib.external.execute import check_on_path


class RAxML():
    """Wrapper for running RAxML."""

    def __init__(self, cpus):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

        self.cpus = cpus

        if self.cpus > 1:
            check_on_path('raxmlHPC-PTHREADS-SSE3')
        else:
            check_on_path('raxmlHPC-SSE3')

    def bootstrap(self, input_tree, msa_file, model_str, num_replicates, output_dir, cpus):
        """Perform non-parametric bootstrapping.

        Parameters
        ----------
        input_tree : str
            File containing newick tree to decorate with bootstraps.
        msa_file : str
            Fasta file containing multiple sequence alignment.
        model_str : str
            Specified either the 'WAG' or 'LG' model.
        num_replicates : int
            Number of replicates to perform.
        output_dir: str
            Output directory to contain bootstrap trees.
        cpus : int
            Number of cpus to use.
        """
        
        check_on_path('seqmagick')

        assert(model_str.upper() in ['WAG', 'LG'])

        self.output_dir = output_dir
        self.model = model_str
        self.msa = seq_io.read(msa_file)

        # calculate replicates
        parallel = Parallel(cpus)
        parallel.run(self._bootstrap, None, range(num_replicates), None)

        # calculate support values
        rep_tree_files = []
        for rep_index in range(num_replicates):
            rep_tree_files.append(os.path.join(output_dir, 'rep_%d' % rep_index, 'RAxML_bestTree.support'))

        tree_name = os.path.splitext(os.path.basename(input_tree))[0]
        output_tree = os.path.join(output_dir, tree_name + '.bootstrap.tree')
        bootstrap_support(input_tree, rep_tree_files, output_tree)
        
        return output_tree

    def _bootstrap(self, rep_num):
        """Infer tree from bootstrapped multiple sequence alignment.

        Parameters
        ----------
        replicated_num : int
          Unique replicate number.
        cpus : int
            Number of cpus to use.
        """
        
        rep_dir = os.path.join(self.output_dir, 'rep_%d' % rep_num)
        if not os.path.exists(rep_dir):
            os.makedirs(rep_dir)

        output_msa = os.path.join(rep_dir, 'bootstrap.afa')
        bootstrap_alignment(self.msa, output_msa)
        
        phylip_msa_file = output_msa.replace('.afa', '.phyx')
        cmd = 'seqmagick convert %s %s' % (output_msa, phylip_msa_file)
        os.system(cmd)

        log_file = os.path.join(rep_dir, 'raxml.log')
        self.run(phylip_msa_file, self.model, rep_dir, log_file, output_suffix='support')

        return True

    def run(self, msa_file, model_str, output_dir, log_file=None, output_suffix=None):
        """Infer tree using RAxML.

        All trees are inferred using the GAMMA distribution to model
        rate heterogeneity. 

        Parameters
        ----------
        msa_file : str
            Phylip file containing multiple sequence alignment.
        model_str : str
            Specified either the 'WAG', 'LG', or 'AUTO' model.
        output_dir: str
            Output directory for inferred tree and auxillary information.
        log_file: str
            Output file containing information about inferred tree.
            
        Returns
        -------
        str
            Path to inferred tree.
        """

        assert(model_str.upper() in ['WAG', 'LG', 'AUTO'])

 
        output_dir = os.path.abspath(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        if not log_file:
            log_file = os.path.join(output_dir, 'raxml.log')
            
        if not output_suffix:
            output_suffix = time.strftime("%Y%m%d")
            
        model_str = 'PROTGAMMA%s' % model_str.upper()
        cmd = '-s %s -m %s -p 12345 -w %s -n %s > %s' % (msa_file,
                                                            model_str,
                                                            output_dir,
                                                            output_suffix,
                                                            log_file)
        if self.cpus > 1:
            cmd = 'raxmlHPC-PTHREADS-SSE3 -T %d ' % self.cpus + cmd
        else:
            cmd = 'raxmlHPC-SSE3 ' + cmd
        os.system(cmd)

        return os.path.join(output_dir, 'RAxML_bestTree.%s' % output_suffix)