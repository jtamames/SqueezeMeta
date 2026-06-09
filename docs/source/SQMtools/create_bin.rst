**********
create_bin
**********

========== ===============
create_bin R Documentation
========== ===============

Create a bin from a vector of contigs
-------------------------------------

Description
~~~~~~~~~~~

Create a bin from a vector of contigs

Usage
~~~~~

.. code:: R

   create_bin(SQM, bin, contigs, delete_overlapping_bins = FALSE)

Arguments
~~~~~~~~~

+-----------------------------+----------------------------------------+
| ``SQM``                     | A SQM object.                          |
+-----------------------------+----------------------------------------+
| ``bin``                     | character. Name of the bin to be       |
|                             | created.                               |
+-----------------------------+----------------------------------------+
| ``contigs``                 | character. Vector with the names of    |
|                             | the contigs that will be included in   |
|                             | the new bin.                           |
+-----------------------------+----------------------------------------+
| ``delete_overlapping_bins`` | logical. If ``TRUE``, bins that        |
|                             | originally contained any of the        |
|                             | provided contigs will be removed from  |
|                             | the object. Default ``FALSE``.         |
+-----------------------------+----------------------------------------+

Value
~~~~~

SQM object with the new binning information, including recalculated bin
statistics if possible.

See Also
~~~~~~~~

``find_redundant_contigs``, ``remove_contigs_from_bin``
