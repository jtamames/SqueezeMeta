#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
from os.path import isdir,isfile
from os import listdir
import os
import sys
import subprocess
import pandas as pd

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
data_path = os.path.abspath(os.path.join(test_dir_path,"test_data","map"))
tmp_dir_path = os.path.join(test_dir_path, 'nose_tmp_output')
tmp_basename_dir = os.path.join(tmp_dir_path, '1')
script_path = os.path.join(test_dir_path, '..', 'scripts', 'gen_input_table.py')

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


    def run_command(self, sample_names = None, contigs_file=None, bam_files=None, bed_files=None, output_file=None):
        
        assert not (bam_files and bed_files)

        
        call = ["python", script_path, 
                "--samplenames", sample_names,
                contigs_file]

        if bam_files:
            call.append(bam_files)
        if bed_files:
            call+= ["--isbedfiles", bed_files]

        call += ['>', output_file]

        self.c = 0
        #try:
        self.op = subprocess.check_output(
                " ".join(call) + " 2> /dev/null",
                shell=True)
        #except subprocess.CalledProcessError as exc:
        #    self.c = exc.returncode

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
 
    def test_with_bamfiles(self):
        bam_files = '*/*-s.bam',
            
        self.run_command(bam_files = data_path + '/*/*-s.bam',
                contigs_file = os.path.join(data_path, 'two_contigs.fa'),
                output_file = os.path.join(tmp_dir_path, 'inputtable.tsv'),
                sample_names = os.path.join(data_path, 'sample_names'))
        assert_equal(self.c, 0,
                     msg = "Command exited with nonzero status")
        
        new_output =  os.path.join(tmp_dir_path, 'inputtable.tsv')
        df = pd.read_csv(new_output, sep='\t', index_col=0)
        assert_almost_equal(df['cov_mean_sample_ten_reads'].loc['contig-75000034'], 10*100.0/1615, 5)
        assert_almost_equal(df['cov_mean_sample_ten_reads'].loc['contig-21000001'], 10*100.0/9998, 5)
        assert_almost_equal(df['cov_mean_sample_twenty_reads'].loc['contig-75000034'], 20*100.0/1615, 5)
        assert_almost_equal(df['cov_mean_sample_twenty_reads'].loc['contig-21000001'], 20*100.0/9998, 5)

            
        #assert_equal(new_output, old_output,
        #        msg = "Output not the same as reference")

    def test_with_bedfiles(self):
        bed_files = '*/*.coverage'
            
        self.run_command(bed_files = data_path + '/*/*.coverage',
                contigs_file = os.path.join(data_path, 'two_contigs.fa'),
                output_file = os.path.join(tmp_dir_path, 'inputtable.tsv'),
                sample_names = os.path.join(data_path, 'sample_names'))
        assert_equal(self.c, 0,
                     msg = "Command exited with nonzero status")
        
        new_output =  os.path.join(tmp_dir_path, 'inputtable.tsv')
        df = pd.read_csv(new_output, sep='\t', index_col=0)
        assert_almost_equal(df['cov_mean_sample_ten_reads'].loc['contig-75000034'], 10*100.0/1615, 5)
        assert_almost_equal(df['cov_mean_sample_ten_reads'].loc['contig-21000001'], 10*100.0/9998, 5)
        assert_almost_equal(df['cov_mean_sample_twenty_reads'].loc['contig-75000034'], 20*100.0/1615, 5)
        assert_almost_equal(df['cov_mean_sample_twenty_reads'].loc['contig-21000001'], 20*100.0/9998, 5)

