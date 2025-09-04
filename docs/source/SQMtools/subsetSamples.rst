*************
subsetSamples
*************

.. container::

   ============= ===============
   subsetSamples R Documentation
   ============= ===============

   .. rubric:: Filter results by sample
      :name: subsetSamples

   .. rubric:: Description
      :name: description

   Create a SQM or SQMlite object containing only samples specified by
   the user, and the ORFs, contigs, bins, taxa and functions present in
   those samples.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetSamples(SQM, samples, remove_missing = TRUE, new_sample_names = NULL)

   .. rubric:: Arguments
      :name: arguments

   +----------------------+----------------------------------------------+
   | ``SQM``              | SQM or SQMlite object to be subsetted.       |
   +----------------------+----------------------------------------------+
   | ``samples``          | character. Samples to be included in the     |
   |                      | subset.                                      |
   +----------------------+----------------------------------------------+
   | ``remove_missing``   | bool. If ``TRUE``, ORFs, contigs, bins, taxa |
   |                      | and functions absent from the selected       |
   |                      | samples will be removed from the subsetted   |
   |                      | object (default ``TRUE``).                   |
   +----------------------+----------------------------------------------+
   | ``new_sample_names`` | character. New sample names to be included   |
   |                      | in the subset (default ``NULL``).            |
   +----------------------+----------------------------------------------+

   .. rubric:: Value
      :name: value

   SQM or SQMlite object containing only the requested samples.

   .. rubric:: See Also
      :name: see-also

   ``subsetTax``, ``subsetFun``, ``subsetORFs``, ``combineSQM``. The
   most abundant items of a particular table contained in a SQM object
   can be selected with ``mostAbundant``.
