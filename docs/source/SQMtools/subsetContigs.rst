*************
subsetContigs
*************

.. container::

   ============= ===============
   subsetContigs R Documentation
   ============= ===============

   .. rubric:: Select contigs
      :name: subsetContigs

   .. rubric:: Description
      :name: description

   Create a SQM object containing only the requested contigs, the ORFs
   contained in them and the bins that contain them.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetContigs(
        SQM,
        contigs,
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = FALSE,
        rescale_copy_number = FALSE,
        recalculate_bin_stats = TRUE,
        allow_empty = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``SQM``                          | SQM object to be subsetted.      |
   +----------------------------------+----------------------------------+
   | ``contigs``                      | character. Vector of contigs to  |
   |                                  | be selected.                     |
   +----------------------------------+----------------------------------+
   | ``trusted_functions_only``       | logical. If ``TRUE``, only       |
   |                                  | highly trusted functional        |
   |                                  | annotations (best hit + best     |
   |                                  | average) will be considered when |
   |                                  | generating aggregated function   |
   |                                  | tables. If ``FALSE``, best hit   |
   |                                  | annotations will be used         |
   |                                  | (default ``FALSE``).             |
   +----------------------------------+----------------------------------+
   | `                                | logical. If ``FALSE``, ORFs with |
   | `ignore_unclassified_functions`` | no functional classification     |
   |                                  | will be aggregated together into |
   |                                  | an "Unclassified" category. If   |
   |                                  | ``TRUE``, they will be ignored   |
   |                                  | (default ``FALSE``).             |
   +----------------------------------+----------------------------------+
   | ``rescale_tpm``                  | logical. If ``TRUE``, TPMs for   |
   |                                  | KEGGs, COGs, and PFAMs will be   |
   |                                  | recalculated (so that the TPMs   |
   |                                  | in the subset actually add up to |
   |                                  | 1 million). Otherwise,           |
   |                                  | per-function TPMs will be        |
   |                                  | calculated by aggregating the    |
   |                                  | TPMs of the ORFs annotated with  |
   |                                  | that function, and will thus     |
   |                                  | keep the scaling present in the  |
   |                                  | parent object (default           |
   |                                  | ``FALSE``).                      |
   +----------------------------------+----------------------------------+
   | ``rescale_copy_number``          | logical. If ``TRUE``, copy       |
   |                                  | numbers with be recalculated     |
   |                                  | using the RecA/RadA coverages in |
   |                                  | the subset. Otherwise, RecA/RadA |
   |                                  | coverages will be taken from the |
   |                                  | parent object. By default it is  |
   |                                  | set to ``FALSE``, which means    |
   |                                  | that the returned copy numbers   |
   |                                  | for each function will represent |
   |                                  | the average copy number of that  |
   |                                  | function per genome in the       |
   |                                  | parent object.                   |
   +----------------------------------+----------------------------------+
   | ``recalculate_bin_stats``        | logical. If ``TRUE``, bin stats  |
   |                                  | and taxonomy are recalculated    |
   |                                  | based on the contigs present in  |
   |                                  | the subsetted object (default    |
   |                                  | ``TRUE``).                       |
   +----------------------------------+----------------------------------+
   | ``allow_empty``                  | (internal use only).             |
   +----------------------------------+----------------------------------+

   .. rubric:: Value
      :name: value

   SQM object containing only the selected contigs.

   .. rubric:: See Also
      :name: see-also

   ``subsetORFs``

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Which contigs have a GC content below 40?
      lowGCcontigNames = rownames(Hadza$contigs$table[Hadza$contigs$table[,"GC perc"]<40,])
      lowGCcontigs = subsetContigs(Hadza, lowGCcontigNames)
      hist(lowGCcontigs$contigs$table[,"GC perc"])
