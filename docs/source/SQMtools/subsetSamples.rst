*************
subsetSamples
*************

============= ===============
subsetSamples R Documentation
============= ===============

Filter results by sample
------------------------

Description
~~~~~~~~~~~

Create a SQM or SQMlite object containing only samples specified by the
user, and the ORFs, contigs, bins, taxa and functions present in those
samples.

Usage
~~~~~

.. code:: R

   subsetSamples(SQM, samples, remove_missing = TRUE, new_sample_names = NULL)

Arguments
~~~~~~~~~

+----------------------+-----------------------------------------------+
| ``SQM``              | SQM or SQMlite object to be subsetted.        |
+----------------------+-----------------------------------------------+
| ``samples``          | character. Samples to be included in the      |
|                      | subset.                                       |
+----------------------+-----------------------------------------------+
| ``remove_missing``   | bool. If ``TRUE``, ORFs, contigs, bins, taxa  |
|                      | and functions absent from the selected        |
|                      | samples will be removed from the subsetted    |
|                      | object (default ``TRUE``).                    |
+----------------------+-----------------------------------------------+
| ``new_sample_names`` | character. New sample names to be included in |
|                      | the subset (default ``NULL``).                |
+----------------------+-----------------------------------------------+

Value
~~~~~

SQM or SQMlite object containing only the requested samples.

See Also
~~~~~~~~

``subsetTax``, ``subsetFun``, ``subsetORFs``, ``combineSQM``. The most
abundant items of a particular table contained in a SQM object can be
selected with ``mostAbundant``.
