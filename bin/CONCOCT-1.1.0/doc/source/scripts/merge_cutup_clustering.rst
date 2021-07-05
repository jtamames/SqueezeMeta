=========================
merge_cutup_clustering.py
=========================

Usage
=====
The usage and help documentation of ``merge_cutup_clustering.py`` can be seen by
running ``merge_cutup_clustering.py -h``:

.. program-output:: cat ../../scripts/merge_cutup_clustering.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``merge_cutup_clustering.py``::

    merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

This merges the clustering ``clustering_gt1000.csv`` created by concoct by looking at cluster assignments per contig part and assigning a concensus cluster for the original contig.
The output clustering_merged.csv contains a header line and contig_id and cluster_id per line, separated by a comma.
