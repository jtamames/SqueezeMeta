#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as p
import hashlib

file_path = os.path.realpath(__file__)
data_path = os.path.abspath(os.path.join(file_path,"..","..","data/"))
test_dir_path = os.path.dirname(file_path)
tmp_dir_path = test_dir_path + '/nose_tmp_output'
tmp_basename_dir = tmp_dir_path + '/1'
script_path = os.path.join(test_dir_path, '..', 'scripts', 'COG_table.py')

CWD = os.getcwd()

class TestCMD(object):
    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()
        os.mkdir(tmp_basename_dir)
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


    def run_command(self, 
            cluster_file = 'prodigal_example/clustering_gt1000.csv', 
            blastoutfile = 'prodigal_example/blast_output.out',
            marker_file = '../../scgs/scg_cogs_min0.97_max1.03_unique_genera.txt',
            cdd_cog_file = '../../scgs/cdd_to_cog.tsv',
            gfffile = None,
            scovs_threshold = None, 
            pident_threshold = None,
            output_file = os.path.join(tmp_dir_path, 'cog_table_output.csv')):
        
        call = ["python", script_path, 
                "--cluster_file", "test_data/{0}".format(cluster_file),
                "--blastoutfile", "test_data/{0}".format(blastoutfile),
                "--marker_file", "test_data/{0}".format(marker_file),
                "--cdd_cog_file", "test_data/{0}".format(cdd_cog_file)]

        if gfffile:
            call += ['--gfffile', "test_data/{0}".format(gfffile)]
        
        if scovs_threshold:
            call += ['--scovs-threshold', scovs_threshold]
        
        if pident_threshold:
            call += ['--pident-threshold', pident_threshold]
        
        call += [">", output_file]

        self.c = 0
        try:
            self.op = subprocess.check_output(
                " ".join(call) + " 2> /dev/null",
                shell=True)
        except subprocess.CalledProcessError as exc:
            self.c = exc.returncode

    def file_len(self,fh):
        i=0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def md5sum(self, filename):
        with open(filename, 'rb') as fh:
             content = fh.read()
        m = hashlib.md5() 
        m.update(content)
        return m.hexdigest()
 
    def test_without_gff(self):
        self.run_command()
        assert_equal(self.c, 0,
                     msg = "Command exited with nonzero status")

        new_output = self.md5sum(os.path.join(tmp_dir_path, 'cog_table_output.csv'))
        old_output = self.md5sum(os.path.join(test_dir_path, 'test_data', 'prodigal_example', 'prodigal_reference_out.tsv'))
        
        assert_equal(new_output, old_output,
                msg = "Output not the same as reference")

    def test_with_gff(self):
        self.run_command(cluster_file = 'prokka_example/clustering_gt1000.csv',
                blastoutfile = 'prokka_example/prokka_rpsblast.tsv',
                gfffile = 'prokka_example/prokka.gff')
        assert_equal(self.c, 0,
                msg = "Command exited with nonzero status")
        new_output = self.md5sum(os.path.join(tmp_dir_path, 'cog_table_output.csv'))
        old_output = self.md5sum(os.path.join(test_dir_path, 'test_data', 'prokka_example', 'prokka_reference_out.tsv'))
        assert_equal(new_output, old_output,
                msg = "Output not the same as reference")

