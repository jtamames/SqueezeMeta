=========================
concoct_coverage_table.py
=========================

Usage
=====
The usage and help documentation of ``concoct_coverage_table.py`` can be seen by
running ``concoct_coverage_table.py -h``:

.. program-output:: cat ../../scripts/concoct_coverage_table.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``concoct_coverage_table.py``::

    concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv

This creates a coverage table suitable as input for concoct as the `coverage_file` parameter.
The ``contigs_10K.bed`` file is created from the ``cut_up_fasta.py`` script and the ``bam``-files needs to be sorted and indexed.
