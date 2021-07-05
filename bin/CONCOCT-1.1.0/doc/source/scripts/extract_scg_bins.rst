================================
[Deprecated] extract_scg_bins.py
================================

Usage
=====
The usage and help documentation of ``extract_scg_bins.py`` can be seen by
running ``pyhton extract_scg_bins -h``:

.. program-output:: cat ../../scripts/extract_scg_bins.py | sed 's/import argparse/import argparse, conf/' | python - --help
   :shell:

Example
=======
An example of how to run ``extract_scg_bins`` on the test data::

    cd CONCOCT/scripts/tests/test_data
    python extract_scg_bins.py \
        --output_folder test_extract_scg_bins_out \
        --scg_tsvs tests/test_data/scg_bins/sample0_gt300_scg.tsv \
                   tests/test_data/scg_bins/sample0_gt500_scg.tsv \
        --fasta_files tests/test_data/scg_bins/sample0_gt300.fa \
                      tests/test_data/scg_bins/sample0_gt500.fa \
        --names sample0_gt300 sample0_gt500 \
        --max_missing_scg 2 --max_multicopy_scg 4 \
        --groups gt300 gt500

This results in the following output files in the folder ``test_extraxt_scg_bins_out/``::

    $ ls test_extract_scg_bins_out/
    sample0_gt300_bin2.fa  sample0_gt500_bin2.fa

Only bin2 satisfies the given criteria for both binnings. If we want to get the
best binning of the two, one can remove the ``--groups`` parameter (or give
them the same group id). That would only output ``sample0_gt500_bin2.fa``,
because the sum of bases in the approved bins of ``sample0_gt500`` is higher
than that of ``sample0_gt300``.
