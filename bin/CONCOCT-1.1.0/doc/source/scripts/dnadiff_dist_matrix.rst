======================
dnadiff_dist_matrix.py
======================

Usage
=====
The usage and help documentation of ``dnadiff_dist_matrix.py`` can be seen by
running ``pyhton dnadiff_dist_matrix -h``:

.. program-output:: cat ../../scripts/dnadiff_dist_matrix.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``dnadiff_dist_matrix`` on the test data::

    cd CONCOCT/scripts
    python dnadiff_dist_matrix.py test_dnadiff_out tests/test_data/bins/sample*.fa

This results in the following output files in the folder ``test_dnadiff_out/``:

    - ``dist_matrix.stv`` The distance matrix
    - ``fasta_names.tsv`` The names given to each bin (or fasta file)
    - ``clustering.tsv`` This file will give a cluster assignment for each bin (or fasta file)
    - :download:`hcust_dendrogram.pdf <../_static/scripts/dna_diff_dist_matrix/hclust_dendrogram.pdf>`
      Dendrogram of the clustering (click for example)
    - :download:`hcust_heatmap.pdf <../_static/scripts/dna_diff_dist_matrix/hclust_heatmap.pdf>`
      Heatmap of the clustering (click for example)

Then there is also for each pairwise ``dnadiff`` alignment the following output
files in a subfolder ``fastaname1_vs_fastaname2/``::

    out.1coords
    out.1delta
    out.cmd
    out.delta
    out.mcoords
    out.mdelta
    out.qdiff
    out.rdiff
    out.report
    out.snps
    out.unqry
    out.unref

See MUMmer's own manual for an explanation of each file with ``dnadiff --help``.
