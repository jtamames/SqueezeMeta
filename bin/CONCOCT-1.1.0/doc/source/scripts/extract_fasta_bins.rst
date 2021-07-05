======================
extract_fasta_bins.py
======================

Usage
=====
The usage and help documentation of ``extract_fasta_bins.py`` can be seen by
running ``extract_fasta_bins.py -h``:

.. program-output:: cat ../../scripts/extract_fasta_bins.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``extract_fasta_bins.py``::

    mkdir concoct_output/fasta_bins
    extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins

This creates a fasta file for each cluster assigned by concoct.
The clusters assigned need not to be complete or uncontaminated and should be investigated closer with e.g. CheckM.
