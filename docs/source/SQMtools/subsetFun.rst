*********
subsetFun
*********

.. container::

   ========= ===============
   subsetFun R Documentation
   ========= ===============

   .. rubric:: Filter results by function
      :name: subsetFun

   .. rubric:: Description
      :name: description

   Create a SQM or SQMbunch object containing only the ORFs with a given
   function, and the contigs and bins that contain them.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetFun(
        SQM,
        fun,
        columns = NULL,
        ignore_case = TRUE,
        fixed = FALSE,
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = FALSE,
        rescale_copy_number = FALSE,
        recalculate_bin_stats = FALSE,
        allow_empty = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``SQM``                          | SQM or SQMbunch object to be     |
   |                                  | subsetted.                       |
   +----------------------------------+----------------------------------+
   | ``fun``                          | character. Pattern to search for |
   |                                  | in the different functional      |
   |                                  | classifications.                 |
   +----------------------------------+----------------------------------+
   | ``columns``                      | character. Restrict the search   |
   |                                  | to the provided column names     |
   |                                  | from ``SQM$orfs$table``. If not  |
   |                                  | provided the search will be      |
   |                                  | performed in all the columns     |
   |                                  | containing functional            |
   |                                  | information (default ``NULL``).  |
   +----------------------------------+----------------------------------+
   | ``ignore_case``                  | logical Make pattern matching    |
   |                                  | case-insensitive (default        |
   |                                  | ``TRUE``).                       |
   +----------------------------------+----------------------------------+
   | ``fixed``                        | logical. If ``TRUE``, pattern is |
   |                                  | a string to be matched as is. If |
   |                                  | ``FALSE`` the pattern is treated |
   |                                  | as a regular expression (default |
   |                                  | ``FALSE``).                      |
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
   |                                  | using the median single-copy     |
   |                                  | gene coverages in the subset.    |
   |                                  | Otherwise, single-copy gene      |
   |                                  | coverages will be taken from the |
   |                                  | parent object. By default it is  |
   |                                  | set to ``FALSE``, which means    |
   |                                  | that the returned copy numbers   |
   |                                  | for each function will represent |
   |                                  | the average copy number of that  |
   |                                  | function per genome in the       |
   |                                  | parent object.                   |
   +----------------------------------+----------------------------------+
   | ``recalculate_bin_stats``        | logical. If ``TRUE``, bin        |
   |                                  | abundance, quality and taxonomy  |
   |                                  | are recalculated based on the    |
   |                                  | contigs present in the subsetted |
   |                                  | object (default ``FALSE``).      |
   +----------------------------------+----------------------------------+
   | ``allow_empty``                  | (internal use only).             |
   +----------------------------------+----------------------------------+

   .. rubric:: Value
      :name: value

   SQM or SQMbunch object containing only the requested function.

   .. rubric:: See Also
      :name: see-also

   ``subsetTax``, ``subsetORFs``, ``subsetSamples``, ``combineSQM``. The
   most abundant items of a particular table contained in a SQM object
   can be selected with ``mostAbundant``.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      Hadza.iron = subsetFun(Hadza, "iron")
      Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
      # Search for multiple patterns using regular expressions
      Hadza.twoKOs = subsetFun(Hadza, "K00812|K00813", fixed=FALSE)
