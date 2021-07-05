import sys
import os
from os.path import join as ospj
from nose.tools import ok_, assert_equal
import numpy as np

import concoct.utils.dir_utils as dir_utils

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data", "bins"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, 'nose_tmp_output')
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, 'dnadiff')
SCRIPT_PATH = ospj(TEST_DIR_PATH, '..')

# Add script dir to python path to import functions
sys.path.append(SCRIPT_PATH)
import dnadiff_dist_matrix

CWD = os.getcwd()


class TestDnaDiff(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        dir_utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temporary output files"""
        dir_utils.rm_rf(TMP_DIR_PATH)

    def test_run_dnadiff(self):
        """Test single dnadiff run"""
        dnadiff_dist_matrix.run_dnadiff(ospj(DATA_PATH,
            "sample0_gt1000_bin0.fa"), ospj(DATA_PATH,
                "sample0_gt1000_bin1.fa"), ospj(TMP_BASENAME_DIR, "out"))
        ok_(os.path.exists(ospj(TMP_BASENAME_DIR, "out.report")))

    def test_run_dnadiff_pairwise(self):
        """Test dnadiff pairwise on multiple bins"""
        names = ["bin{0}".format(i) for i in range(3)]
        dnadiff_dist_matrix.run_dnadiff_pairwise(
            [ospj(DATA_PATH, b) for b in
                ["sample0_gt1000_bin0.fa", "sample0_gt1000_bin1.fa",
                    "sample0_gt1000_bin2.fa"]],
            names, TMP_BASENAME_DIR)
        for f in [
            ospj(TMP_BASENAME_DIR,
                "{}_vs_{}".format(names[i], names[j]), "out.report")
                for i in range(len(names))
                for j in range(i + 1, len(names))]:
            ok_(os.path.exists(f))

    def test_parallel_run_dnadiff_pairwise(self):
        """Test dnadiff pairwise on multiple bins"""
        names = ["bin{0}".format(i) for i in range(3)]
        dnadiff_dist_matrix.parallel_run_dnadiff_pairwise(
            [ospj(DATA_PATH, b) for b in
                ["sample0_gt1000_bin0.fa", "sample0_gt1000_bin1.fa",
                    "sample0_gt1000_bin2.fa"]],
            names, TMP_BASENAME_DIR)
        for f in [
            ospj(TMP_BASENAME_DIR,
                "{}_vs_{}".format(names[i], names[j]), "out.report")
                for i in range(len(names))
                for j in range(i + 1, len(names))]:
            ok_(os.path.exists(f))

    def test_mummer_report_class(self):
        """Test mummer report class"""
        dnadiff_dist_matrix.run_dnadiff(ospj(DATA_PATH,
            "sample0_gt1000_bin0.fa"), ospj(DATA_PATH,
                "sample0_gt1000_bin1.fa"), ospj(TMP_BASENAME_DIR, "out"))
        ok_(os.path.exists(ospj(TMP_BASENAME_DIR, "out.report")))

        mumr = dnadiff_dist_matrix.MUMmerReport(ospj(TMP_BASENAME_DIR,
            "out.report"))
        assert_equal(mumr.tot_bases[0], 3213)
        assert_equal(mumr.tot_bases[1], 43514)
        assert_equal(mumr.aligned_bases[0], 0)
        assert_equal(mumr.aligned_bases[1], 0)

    def test_get_dist_matrix(self):
        """Test dnadiff pairwise on multiple bins and get the resulting
        distance matrix"""
        names = [
            "sample0_gt1000_bin0",
            "sample0_gt1000_bin10",
            "sample0_gt1000_bin11",
            "sample0_gt1000_bin1",
            "sample0_gt1000_bin2",
            "sample1_gt1000_bin51",
            "sample1_gt1000_bin67",
            "sample1_gt1000_bin6",
        ]
        files = [ospj(DATA_PATH, "{}.fa".format(n)) for n in names]
        dnadiff_dist_matrix.run_dnadiff_pairwise(files, names,
                TMP_BASENAME_DIR)
        for f in [
            ospj(TMP_BASENAME_DIR,
                "{}_vs_{}".format(names[i], names[j]), "out.report")
                for i in range(len(names))
                for j in range(i + 1, len(names))]:
            ok_(os.path.exists(f))

        matrix = dnadiff_dist_matrix.get_dist_matrix(TMP_BASENAME_DIR, names,
                50)
        matrix_exp = np.genfromtxt(ospj(DATA_PATH, "expected_dist_matrix.tsv"),
                delimiter="\t")
        np.testing.assert_almost_equal(matrix, matrix_exp, decimal=2)

    def test_plot_dist_matrix(self):
        """Plot a distance matrix"""
        names = [
            "sample0_gt1000_bin0",
            "sample0_gt1000_bin10",
            "sample0_gt1000_bin11",
            "sample0_gt1000_bin1",
            "sample0_gt1000_bin2",
            "sample1_gt1000_bin51",
            "sample1_gt1000_bin67",
            "sample1_gt1000_bin6",
        ]
        matrix = np.genfromtxt(ospj(DATA_PATH, "expected_dist_matrix.tsv"),
                delimiter="\t")
        heatmap = ospj(TMP_BASENAME_DIR, "hclust_heatmap.pdf")
        dendrogram = ospj(TMP_BASENAME_DIR, "hclust_dendrogram.pdf")
        clustering = ospj(TMP_BASENAME_DIR, "clustering.tsv")
        clustering_threshold = 0.05
        dnadiff_dist_matrix.plot_dist_matrix(matrix, names, heatmap, dendrogram, clustering_threshold, clustering)
        ok_(os.path.exists(heatmap))
        ok_(os.path.exists(dendrogram))
        ok_(os.path.exists(clustering))

    def test_plot_dist_matrix_88_bins(self):
        """Plot a distance matrix with 88 samples"""
        names = [
            "sample0_gt1000_bin{}".format(i) for i in range(88)
        ]
        matrix = np.genfromtxt(ospj(DATA_PATH, "expected_dist_matrix_88_bins.tsv"),
                delimiter="\t")
        heatmap = ospj(TMP_BASENAME_DIR, "hclust_heatmap.pdf")
        dendrogram = ospj(TMP_BASENAME_DIR, "hclust_dendrogram.pdf")
        clustering = ospj(TMP_BASENAME_DIR, "clustering.tsv")
        clustering_threshold = 0.05
        dnadiff_dist_matrix.plot_dist_matrix(matrix, names, heatmap, dendrogram, clustering_threshold, clustering)
        ok_(os.path.exists(heatmap))
        ok_(os.path.exists(dendrogram))
        ok_(os.path.exists(clustering))

    def test_write_fasta_names(self):
        names = [
            "sample0_gt1000_bin0",
            "sample0_gt1000_bin10",
            "sample0_gt1000_bin11",
            "sample0_gt1000_bin1",
            "sample0_gt1000_bin2",
            "sample1_gt1000_bin51",
            "sample1_gt1000_bin67",
            "sample1_gt1000_bin6",
        ]
        files = [ospj(DATA_PATH, "{}.fa".format(n)) for n in names]
        dnadiff_dist_matrix.write_fasta_names(names, files,
                ospj(TMP_BASENAME_DIR, "fasta_names.tsv"), "\t")
        ok_(os.path.exists(ospj(TMP_BASENAME_DIR, "fasta_names.tsv")))
