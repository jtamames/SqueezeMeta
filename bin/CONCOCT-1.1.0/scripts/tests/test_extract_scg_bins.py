import sys
import os
from os.path import join as ospj
from nose.tools import assert_equal, ok_
import pandas as pd
import collections

import concoct.utils.dir_utils as dir_utils

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data", "scg_bins"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, 'nose_tmp_output')
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, 'extract_scg_bins')
SCRIPT_PATH = ospj(TEST_DIR_PATH, '..')

# Add script dir to python path to import functions
sys.path.append(SCRIPT_PATH)
from extract_scg_bins import get_approved_bins, sum_bases_in_bins, \
    get_winning_bins, write_approved_bins, main

CWD = os.getcwd()


class TestDnaDiff(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        dir_utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temporary output files"""
        #dir_utils.rm_rf(TMP_DIR_PATH)

    def test_get_approved_bins(self):
        """Test get_approved_bins"""
        df = get_approved_bins(ospj(DATA_PATH, "sample0_gt500_scg.tsv"),
                max_missing_scg=2, max_multicopy_scg=4)
        assert_equal(2, int(df.Cluster))

    def test_sum_bases_in_bins(self):
        """Test sum_bases_in_bins"""
        scg_tsv = ospj(DATA_PATH, "sample0_gt500_scg.tsv")
        b = sum_bases_in_bins(pd.read_csv(scg_tsv, sep="\t"),
                ospj(DATA_PATH, "sample0_gt500.fa"))
        assert_equal(12, b)
        df = get_approved_bins(ospj(DATA_PATH, "sample0_gt500_scg.tsv"),
                max_missing_scg=2, max_multicopy_scg=4)
        b = sum_bases_in_bins(df, ospj(DATA_PATH, "sample0_gt500.fa"))
        assert_equal(4, b)

    def test_get_winning_bins(self):
        """Test get_winning_bins"""
        scg_tsvs = [ospj(DATA_PATH, p) for p in ["sample0_gt300_scg.tsv",
            "sample0_gt500_scg.tsv"]]
        fasta_files = [ospj(DATA_PATH, p) for p in ["sample0_gt300.fa",
            "sample0_gt500.fa"]]
        winning_index, df = get_winning_bins(scg_tsvs, fasta_files,
                max_missing_scg=2, max_multicopy_scg=4)
        assert_equal(1, winning_index)
        winning_index, df = get_winning_bins(list(reversed(scg_tsvs)),
            list(reversed(fasta_files)), max_missing_scg=2,
            max_multicopy_scg=4)
        assert_equal(0, winning_index)

    def test_write_approved_bins(self):
        """Test write_approved_bins"""
        df = get_approved_bins(ospj(DATA_PATH, "sample0_gt500_scg.tsv"),
                max_missing_scg=2, max_multicopy_scg=4)
        assert_equal(2, int(df.Cluster))
        write_approved_bins(df, ospj(DATA_PATH, "sample0_gt500.fa"),
                TMP_BASENAME_DIR, "sample0_gt500")
        ok_(os.path.exists(ospj(TMP_BASENAME_DIR, "sample0_gt500_bin2.fa")))
        # make sure both have equal amount of records
        assert_equal(
            open(ospj(TMP_BASENAME_DIR, "sample0_gt500_bin2.fa")).read().count(">"),
            open(ospj(DATA_PATH, "sample0_gt500_bin2.fa")).read().count(">"))

    def test_find_best_per_group(self):
        fasta_files = [
            ospj(DATA_PATH, "sample0_gt300.fa"),
            ospj(DATA_PATH, "sample0_gt500.fa"),
        ]
        args = collections.namedtuple('Arguments', " ".join(["output_folder",
            "scg_tsvs", "fasta_files", "names", "max_missing_scg",
            "max_multicopy_scg", "groups"]))
        groupargs = args(
            output_folder=TMP_BASENAME_DIR,
            scg_tsvs=[os.path.splitext(f)[0] + "_scg.tsv" for f in fasta_files],
            fasta_files=fasta_files,
            names=[os.path.splitext(os.path.basename(f))[0] for f in fasta_files],
            max_missing_scg=2,
            max_multicopy_scg=4,
            groups=("gt300", "gt500")
        )
        main(groupargs)
