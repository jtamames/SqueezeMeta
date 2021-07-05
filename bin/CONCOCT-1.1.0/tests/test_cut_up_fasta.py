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
data_path = os.path.abspath(os.path.join(test_dir_path,"test_data"))
tmp_dir_path = os.path.join(test_dir_path, 'nose_tmp_output')
tmp_basename_dir = os.path.join(tmp_dir_path, '1')
script_path = os.path.join(test_dir_path, '..', 'scripts', 'cut_up_fasta.py')

CWD = os.getcwd()

CONTIG_PART_EXPR = re.compile("(.*)\.concoct_part_([0-9]*)")

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


    def run_command(self, contigs_file=None, chunk_size=None, overlap_size=None, merge_last=None, bedfile=None, output_file=None):

        call = ["python", script_path, contigs_file,
                "--chunk_size", str(chunk_size),
                "--overlap_size", str(overlap_size)]

        if merge_last:
            call.append("--merge_last")

        if bedfile:
            call += ["--bedfile", bedfile]

        call += ['>', output_file]

        self.c = 0
        print(" ".join(call))
        self.op = subprocess.check_output(
                " ".join(call) + " 2> /dev/null",
                shell=True)

    def file_len(self,fh):
        i=0
        with open(fh) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_basic(self):
        contigs_file = os.path.join(data_path, 'composition.fa')

        output_file = os.path.join(tmp_dir_path, 'contigs_cutup.fa')

        self.run_command(
                contigs_file=contigs_file,
                chunk_size=200,
                overlap_size=0,
                merge_last=True,
                bedfile=None,
                output_file=output_file)

        assert_equal(self.c, 0,
                     msg = "Command exited with nonzero status")

        result_dict = {}
        nr_parts_per_original = defaultdict(int)
        part_lengths_per_original = defaultdict(int)
        reconstructed_original = defaultdict(str)

        with open(output_file) as fhandle:
            for rec in SeqIO.parse(fhandle, "fasta"):
                assert_true(".concoct_part_" in rec.id)
                original_id, part_nr = CONTIG_PART_EXPR.match(rec.id).group(1,2)
                assert_true(".concoct_part_" not in original_id)
                nr_parts_per_original[original_id] += 1
                part_lengths_per_original[original_id] += len(rec.seq)
                assert_true((len(rec.seq) == 200) or (len(rec.seq) > 200 and len(rec.seq) < 400))

                reconstructed_original[original_id] += str(rec.seq)

        with open(contigs_file) as fhandle:
            for rec in SeqIO.parse(fhandle, "fasta"):
                assert_true(rec.id in reconstructed_original.keys())
                assert_equal(str(rec.seq), reconstructed_original[rec.id])
                assert_equal(int(len(rec.seq) / 200), \
                             nr_parts_per_original[rec.id])


    def test_with_bedfile(self):
        bedfile = os.path.join(tmp_dir_path, 'contigs_cutup.bed'),
