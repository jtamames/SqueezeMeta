***********************
remove_contigs_from_bin
***********************

======================= ===============
remove_contigs_from_bin R Documentation
======================= ===============

Remove contigs from a given bin
-------------------------------

Description
~~~~~~~~~~~

Remove contigs from a given bin

Usage
~~~~~

.. code:: R

   remove_contigs_from_bin(SQM, bin, contigs)

Arguments
~~~~~~~~~

+-------------+--------------------------------------------------------+
| ``SQM``     | A SQM object.                                          |
+-------------+--------------------------------------------------------+
| ``bin``     | character. Name of the bin from which the contigs will |
|             | be removed.                                            |
+-------------+--------------------------------------------------------+
| ``contigs`` | character. Vector with the names of the contigs that   |
|             | will be removed from the new bin.                      |
+-------------+--------------------------------------------------------+

Value
~~~~~

SQM object with the new binning information, including recalculated bin
statistics if possible.

See Also
~~~~~~~~

``find_redundant_contigs``, ``create_bin``
