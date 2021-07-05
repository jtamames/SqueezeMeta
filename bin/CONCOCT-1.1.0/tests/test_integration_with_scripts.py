#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest, assert_false
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as p
import glob

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path, "..", "test_data", "integration_test_data"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'

CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        os.chdir(test_dir_path)

    def tearDown(self):
        """remove temporary output files"""
        for d in os.listdir(tmp_dir_path):
            d_path = os.path.join(tmp_dir_path,d)
            try:
                os.remove(d_path)
            except:
                for f in os.listdir(d_path):
                    f_path = os.path.join(d_path,f)
                    os.remove(f_path)
                os.rmdir(d_path)
        assert os.listdir(tmp_dir_path) == []


    def run_command(self, call_list):
        call_string = " ".join(call_list) + " 2> /dev/null"
        self.op = subprocess.check_output(
                        call_string,
                        shell=True)

    def run_all(self, cov_file='coverage', comp_file='composition.fa',
                    tags=[], basename='nose_tmp_output/1'):

        original_contigs = os.path.join(data_path, "velvet_71.fa")
        contigs_bed = os.path.join(tmp_dir_path, "contigs.bed")
        cutup_contigs = os.path.join(tmp_dir_path, "contigs_c10K.fa")
        bam_files = glob.glob(os.path.join(data_path, "map", "*.bam"))
        coverage_table = os.path.join(tmp_dir_path, "coverage_table.tsv")
        concoct_output = os.path.join(tmp_dir_path, "concoct_output") + "/"
        clustering_file = os.path.join(concoct_output, "clustering_gt1000.csv")
        merged_clustering = os.path.join(concoct_output, "clustering_merged.csv")
        fasta_bins_dir = os.path.join(tmp_dir_path, "fasta_bins")

        cutup_call = ["cut_up_fasta.py", original_contigs, "-c", "10000", "-o", "0", "--merge_last", "-b", contigs_bed, ">", cutup_contigs]
        self.run_command(cutup_call)

        coverage_table_call = ["concoct_coverage_table.py", contigs_bed] + bam_files + [">", coverage_table]
        self.run_command(coverage_table_call)

        concoct_call = ["concoct", "--composition_file", cutup_contigs, "--coverage_file", coverage_table, "-b", concoct_output]
        self.run_command(concoct_call)

        merge_cutup_call = ["merge_cutup_clustering.py", clustering_file, ">", merged_clustering]
        self.run_command(merge_cutup_call)

        _ = subprocess.check_output("mkdir {}".format(fasta_bins_dir), shell=True)

        extract_call = ["extract_fasta_bins.py", original_contigs, merged_clustering, "--output_path", fasta_bins_dir]
        self.run_command(extract_call)

    def test_directory_creation(self):
        self.run_all()
        fasta_bins_dir = os.path.join(tmp_dir_path, "fasta_bins")
        print(os.listdir(fasta_bins_dir))
        assert_true(len(glob.glob(os.path.join(fasta_bins_dir, "*.fa"))) > 2, "Too few bins were created")
