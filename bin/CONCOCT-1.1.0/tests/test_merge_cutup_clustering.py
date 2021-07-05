#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
from Bio import SeqIO
from collections import defaultdict
import re

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
data_path = os.path.abspath(os.path.join(test_dir_path,"test_data", "integration_test_data"))
tmp_dir_path = os.path.join(test_dir_path, 'nose_tmp_output')
script_path = os.path.join(test_dir_path, '..', 'scripts', 'merge_cutup_clustering.py')

CWD = os.getcwd()

CONTIG_PART_EXPR = re.compile("(.*)\.concoct_part_([0-9]*)")

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


    def run_command(self, clustering_file, output_file):

        call = ["python", script_path, clustering_file]

        call += ['>', output_file]

        self.c = 0
        print(" ".join(call))
        self.op = subprocess.check_output(
                " ".join(call) + " 2> /dev/null",
                shell=True)

    def file_len(self,fh):
        i = 0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_basic(self):
        clustering_file = os.path.join(data_path, 'clustering_gt1000.csv')
        output_file = os.path.join(tmp_dir_path, 'clustering_merged.csv')

        self.run_command(
                clustering_file=clustering_file,
                output_file=output_file
                )

        clusters_per_contig = defaultdict(set)
        with open(clustering_file) as fhandle:
            first = True
            for line in fhandle:
                # First line is header
                if first:
                    first = False
                    continue
                line = line.strip()
                subcontig_id, cluster_id = line.split(',')
                contig_id, part_nr = CONTIG_PART_EXPR.match(subcontig_id).group(1,2)
                clusters_per_contig[contig_id].add(cluster_id)
                if len(clusters_per_contig[contig_id]) > 1:
                    print(contig_id)

        assert_equal(self.c, 0,
                     msg="Command exited with nonzero status")

        assert_equal(self.file_len(output_file), 1+len(clusters_per_contig))
        with open(output_file) as fhandle:
            for line in fhandle:
                line = line.strip()
                contig_id, cluster_id = line.split(',')
                if len(clusters_per_contig[contig_id]) == 1:
                    assert_true(cluster_id in clusters_per_contig[contig_id])
                else:
                    print(contig_id)
