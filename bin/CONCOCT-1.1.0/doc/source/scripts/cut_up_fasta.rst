======================
cut_up_fasta.py
======================

Usage
=====
The usage and help documentation of ``cut_up_fasta.py`` can be seen by
running ``cut_up_fasta.py -h``:

.. program-output:: cat ../../scripts/cut_up_fasta.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``cut_up_fasta.py``::

    cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

This creates a fasta file and a BED file.
The fasta file ``contigs_10K.fa`` contains the original contigs cut up into parts of length exactly 10K, except for the last contig part which is between 10K and 20K long.
The BED file ``contigs_10K.bed`` contains a list of the contig parts created with coordinates in the original contigs.
