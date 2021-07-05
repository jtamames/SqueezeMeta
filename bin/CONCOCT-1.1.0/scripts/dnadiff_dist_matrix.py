#!/usr/bin/env python
"""
Output distance matrix between fasta files using dnadiff from MUMmer. Generates
dnadiff output files in folders:

output_folder/fastaname1_vs_fastaname2/
output_folder/fastaname1_vs_fastaname3/

etc

where fastaname for each fasta file can be supplied as an option to the script.
Otherwise they are just counted from 0 to len(fastafiles)

The distance between each bin is computed using the 1-to-1 alignments of the
report files (not M-to-M):

1 - AvgIdentity if min(AlignedBases) >= min_coverage. Otherwise distance is 1.
Or 0 to itself.

Resulting matrix is printed to stdout and to output_folder/dist_matrix.tsv. The
rows and columns of the matrix follow the order of the supplied fasta files. The
names given to each fasta file are also outputted to the file
output_folder/fasta_names.tsv

A hierarchical clustering of the distance using euclidean average linkage
clustering is plotted. This can be deactivated by using --skip_plot. The
resulting heatmap is in output_folder/hclust_heatmap.pdf or
output_folder/hclust_dendrogram.pdf and the resulting clustering is presented
in output_folder/clustering.tsv. The image extension can be changed.
"""
import argparse
import subprocess
import re
import os
import sys
from os.path import join as ospj
import numpy as np
from multiprocessing import Pool
import logging

from concoct.utils import dir_utils
from concoct.utils import check_dependencies


class CmdException(Exception):
    """Exception for a shell command that did not exit with returncode 0."""
    def __init__(self, cmd, stdout, stderr, returncode):
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode

    def __str__(self):
        return \
          "cmd:\n{}\n" \
          "returncode:\n{}\n" \
          "stdout:\n{}\n" \
          "stderr:\n{}\n".format(self.cmd, self.returncode, self.stdout,
          self.stderr)


class MUMmerReport(object):
    """Represents .report file from MUMmer's dnadiff. Stores TotalBases and
    AlignedBases (%) stats."""
    def __init__(self, report_file):
        self.report_file = report_file

        with open(report_file) as f:
            one_to_one_parsed = False

            for line in f:
                if line.startswith("TotalBases"):
                    self.tot_bases = [int(b) for b in line.split()[1:]]
                if line.startswith("AlignedBases"):
                    # store percentage
                    self.aligned_bases = [float(p) for p in
                            re.findall(r'\((.*?)\%\)', line)]
                if not one_to_one_parsed and line.startswith("AvgIdentity"):
                    self.avg_identity = [float(p) for p in line.split()[1:]]
                    one_to_one_parsed = True


def run_dnadiff(fasta1, fasta2, prefix):
    """Runs MUMmer's dnadiff"""
    cmd = "dnadiff -p {prefix} {f1} {f2}".format(f1=fasta1, f2=fasta2,
            prefix=prefix)
    with open("{}.cmd".format(prefix), "w") as f:
        f.write(cmd)
    proc = subprocess.Popen(cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, shell=True)
    stdout, stderr = proc.communicate()
    if proc.returncode > 0:
        raise CmdException(cmd, stdout, stderr, proc.returncode)


def run_dnadiff_star(args):
    """Converts given list to arguments for run_dnadiff. Useful for
    multiprocessing. http://stackoverflow.com/questions/5442910"""
    try:
        return run_dnadiff(*args)
    except CmdException as e:
        # Custom CmdException doesn't work well with multiprocessing so change
        # to regular Exception http://bugs.python.org/issue16558
        raise Exception


def run_dnadiff_pairwise(fasta_files, fasta_names, output_folder):
    """Runs MUMmer's dnadiff pairwise for given fasta_files. Uses fasta_names
    to organize output folders for dnadiff as fastaname1_vs_fastaname2."""
    assert len(fasta_files) == len(fasta_names)

    for i in range(len(fasta_files)):
        for j in range(i + 1, len(fasta_files)):
            out_dir = ospj(output_folder, "{fn1}_vs_{fn2}".format(
                fn1=fasta_names[i], fn2=fasta_names[j]))
            dir_utils.mkdir_p(out_dir)
            run_dnadiff(fasta_files[i], fasta_files[j], ospj(out_dir, "out"))


def parallel_run_dnadiff_pairwise(fasta_files, fasta_names, output_folder):
    """Runs MUMmer's dnadiff pairwise for given fasta_files using
    multiprocessing. Uses fasta_names to organize output folders for dnadiff as
    fastaname1_vs_fastaname2."""
    assert len(fasta_files) == len(fasta_names)

    pool = Pool()
    args = []
    for i in range(len(fasta_files)):
        for j in range(i + 1, len(fasta_files)):
            out_dir = ospj(output_folder, "{fn1}_vs_{fn2}".format(
                fn1=fasta_names[i], fn2=fasta_names[j]))
            dir_utils.mkdir_p(out_dir)
            args.append((fasta_files[i], fasta_files[j], ospj(out_dir, "out")))
    pool.map(run_dnadiff_star, args)
    pool.close()
    pool.join()


def get_dist_matrix(pairwise_folder, fasta_names, min_coverage):
    """Returns distance matrix from folder constructed with
    run_dnadiff_pairwise"""
    matrix = np.array(len(fasta_names) * [len(fasta_names) * [0.0]])
    for i in range(len(fasta_names)):
        for j in range(i + 1, len(fasta_names)):
            repfile = ospj(pairwise_folder, "{fn1}_vs_{fn2}".format(
                fn1=fasta_names[i], fn2=fasta_names[j]), "out.report")
            mumr = MUMmerReport(repfile)
            if min(mumr.aligned_bases) >= min_coverage:
                # take distance as 1 - AvgIdentity (same for both ref and qry)
                matrix[i][j] = 1.0 - (mumr.avg_identity[0] / 100.0)
            else:
                matrix[i][j] = 1.0
            matrix[j][i] = matrix[i][j]

    return matrix


def plot_dist_matrix(matrix, fasta_names, heatmap_out, dendrogram_out, cluster_threshold, clustering_out):
    """Cluster the distance matrix hierarchically and plot using seaborn.
    Average linkage method is used."""
    # Load required modules for plotting
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

    # Create
    pdm = pd.DataFrame(matrix, index=fasta_names, columns=fasta_names)

    # Create linkage clustering
    link = linkage(pdm, metric='euclidean', method='average')
    flat_clusters = fcluster(link, cluster_threshold, criterion='distance')

    clustering = pd.Series(dict(list(zip(fasta_names, flat_clusters))))
    clustering.to_csv(clustering_out, sep='\t')

    # Plot heatmap
    figsizex = max(10, len(fasta_names) / 4)
    clustergrid = sns.clustermap(pdm, col_linkage=link, row_linkage=link,
            figsize=(figsizex, figsizex))
    clustergrid.savefig(heatmap_out)

    # Plot dendrogram
    sns.set_style('white')
    figsizey = max(10, len(fasta_names) / 8)
    f, ax = plt.subplots(figsize=(figsizex, figsizey))
    dendrogram(link, labels=pdm.index, ax=ax)
    no_spine = {'left': True, 'bottom': True, 'right': True, 'top': True}
    sns.despine(**no_spine)
    plt.xticks(rotation=90)
    f.tight_layout()
    plt.savefig(dendrogram_out)


def write_fasta_names(fasta_names, fasta_files, output_file, separator):
    fnames_files = np.array([fasta_names, fasta_files])
    np.savetxt(output_file, fnames_files.transpose(), fmt="%s", delimiter=separator)


def parse_input():
    """Return input arguments using argparse"""
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output_folder", help="Output folder")
    parser.add_argument("fasta_files", nargs='+',
            help="fasta files to compare pairwise using MUMmer's dnadiff")
    parser.add_argument("--min_coverage", type=float, default=50,
            help="Minimum coverage of bin in percentage to calculate distance "
                 "otherwise distance is 1. Default is 50.")
    parser.add_argument("--fasta_names", default=None, help="File with names "
            "for fasta file, one line each. Could be sample names, bin names, "
            "genome names, whatever you want. The names are used when storing "
            "the MUMmer dnadiff results as in "
            "output_folder/fastaname1_vs_fastaname2/. The names are also used "
            "for the plots.")
    parser.add_argument("--plot_image_extension", default='pdf', help="Type of "
            "image to plotted e.g. pdf, png, svg.")
    parser.add_argument("--skip_dnadiff", action="store_true", help="Skips "
            "running MUMmer and uses output_folder as given input to "
            "calculate the distance matrix. Expects dnadiff output as "
            "output_folder/fastaname1_vs_fastaname2/out.report")
    parser.add_argument("--skip_matrix", action="store_true", help="Skips "
            "Calculating the distance matrix.")
    parser.add_argument("--skip_plot", action="store_true", help="Skips "
        "plotting the distance matrix. By default the distance matrix is "
        "clustered hierarchically using euclidean average linkage clustering. "
        "This step requires seaborn and scipy.")
    parser.add_argument("--cluster-threshold", type=float, default=0.05,
            help=("The maximum within cluster distance allowed."))
    args = parser.parse_args()
    # Get fasta names
    if args.fasta_names is not None:
        fasta_names = [s[:-1] for s in open(args.fasta_names).readlines()]
        if len(fasta_names) != len(args.fasta_files):
            raise Exception("Nr of names in fasta_names should be equal to nr "
                    "of given fasta_files")
    else:
        # Get basename from fasta files and see if those are unique
        fasta_names = [os.path.basename(f).split(".")[0] for f in
                args.fasta_files]
        if len(set(fasta_names)) != len(args.fasta_files):
            # assign integer id if non-unique fasta basenames
            fasta_names = list(range(len(args.fasta_files)))
    if args.skip_dnadiff and args.skip_matrix:
        raise Exception("If running dnadiff and calculating the distance "
                "matrix are both skipped, the program does not run any "
                "steps at all.")

    return args.output_folder, args.fasta_files, fasta_names, \
           args.min_coverage, args.skip_dnadiff, args.skip_matrix, \
           args.skip_plot, args.plot_image_extension, args.cluster_threshold


def verbose_check_dependencies(progs):
    for p in progs:
        path = check_dependencies.which(p)
        if path:
            logging.info("Using {}".format(path))
        else:
            raise Exception


def main(output_folder, fasta_files, fasta_names, min_coverage,
        skip_dnadiff=False, skip_matrix=False, skip_plot=False,
        plot_image_extension="pdf", cluster_threshold=0.05):
    """Output distance matrix between fasta files using MUMmer's dnadiff"""
    # create logger
    logging.basicConfig(
        stream=sys.stdout,
        level=logging.INFO,
        format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
    )
    logging.info("Checking dependencies")
    verbose_check_dependencies(["dnadiff"])
    if not skip_dnadiff:
        logging.info("Running dnadiff pairwise (this could take a while)")
        parallel_run_dnadiff_pairwise(fasta_files, fasta_names, output_folder)
    else:
        logging.info("Skipping dnadiff")
    if not skip_matrix:
        matrix = get_dist_matrix(output_folder, fasta_names, min_coverage)
        # save distance matrix to file
        logging.info("Writing distance matrix to "
                "{}".format(ospj(output_folder, "dist_matrix.tsv")))
        np.savetxt(ospj(output_folder, "dist_matrix.tsv"), matrix, fmt="%.2f",
                delimiter="\t")
    else:
        logging.info("Skipping matrix calculation")
    if not skip_plot:
        if skip_matrix:
            logging.info("Reading matrix from "
                    "{}".format(ospj(output_folder, "dist_matrix.tsv")))
            matrix = np.genfromtxt(ospj(output_folder, "dist_matrix.tsv"),
                    delimiter="\t")
        # plot distance matrix
        plot_dist_matrix(matrix, fasta_names,
            ospj(output_folder,
                "hclust_heatmap.{}".format(plot_image_extension)),
            ospj(output_folder,
                "hclust_dendrogram.{}".format(plot_image_extension)),
            cluster_threshold,
            ospj(output_folder,
                "clustering.tsv"))
    else:
        logging.info("Skipping plotting")
    logging.info("Writing fasta names of fasta files to "
            "{}".format(ospj(output_folder, "fasta_names.tsv")))
    write_fasta_names(fasta_names, fasta_files,
        ospj(output_folder, "fasta_names.tsv"), "\t")
    logging.info("Done")


if __name__ == "__main__":
    main(*parse_input())
