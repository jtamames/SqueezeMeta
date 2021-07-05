CONCOCT Scripts
===============
CONCOCT ships with some additional scripts which are very useful to e.g. create input files and to extract output fastas for concoct.
These scripts are:

    - ``cut_up_fasta.py``
    - ``concoct_coverage_table.py``
    - ``merge_cutup_clustering.py``
    - ``extract_fasta_bins.py``

The repository CONCOCT contains additional scripts in the ``CONCOCT/scripts`` directory which are not fully maintained.
They implement methods that we apply after binning with CONCOCT and it might be useful as a starting point or inspiration when creating your own scripts for downstream processing of the output files.
Out of these scripts, the ones documented here are:

    - ``dnadiff_dist_matrix.py``
    - ``extract_scg_bins.py`` [Deprecated]

Contents:

.. toctree::
    :maxdepth: 2

    cut_up_fasta
    concoct_coverage_table
    merge_cutup_clustering
    extract_fasta_bins
    dnadiff_dist_matrix
    extract_scg_bins
