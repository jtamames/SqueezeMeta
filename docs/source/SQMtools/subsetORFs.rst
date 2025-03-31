**********
subsetORFs
**********

.. container::

   ========== ===============
   subsetORFs R Documentation
   ========== ===============

   .. rubric:: Select ORFs
      :name: subsetORFs

   .. rubric:: Description
      :name: description

   Create a SQM object containing only the requested ORFs, and the
   contigs and bins that contain them. Internally, all the other subset
   functions in this package end up calling ``subsetORFs`` to do the
   work for them.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetORFs(
        SQM,
        orfs,
        tax_source = "orfs",
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = FALSE,
        rescale_copy_number = FALSE,
        recalculate_bin_stats = TRUE,
        contigs_override = NULL,
        allow_empty = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``SQM``                          | SQM object to be subsetted.      |
   +----------------------------------+----------------------------------+
   | ``orfs``                         | character. Vector of ORFs to be  |
   |                                  | selected.                        |
   +----------------------------------+----------------------------------+
   | ``tax_source``                   | character. Features used for     |
   |                                  | calculating aggregated           |
   |                                  | abundances at the different      |
   |                                  | taxonomic ranks. Either          |
   |                                  | ``"orfs"`` or ``"contigs"``      |
   |                                  | (default ``"orfs"``).            |
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
   | ``contigs_override``             | character. Optional vector of    |
   |                                  | contigs to be included in the    |
   |                                  | subsetted object.                |
   +----------------------------------+----------------------------------+
   | ``allow_empty``                  | (internal use only).             |
   +----------------------------------+----------------------------------+

   .. rubric:: Value
      :name: value

   SQM object containing the requested ORFs.

   .. rubric:: A note on contig/bins subsetting
      :name: a-note-on-contigbins-subsetting

   While this function selects the contigs and bins that contain the
   desired orfs, it DOES NOT recalculate contig abundance and statistics
   based on the selected ORFs only. This means that the abundances
   presented in tables such as ``SQM$contig$abund`` will still refer to
   the complete contigs, regardless of whether only a fraction of their
   ORFs are actually present in the returned SQM object. This is also
   true for the statistics presented in ``SQM$contigs$table``. Bin
   statistics may be recalculated if ``rescale_copy_number`` is set to
   ``TRUE``, but recalculation will be based on contigs, not ORFs.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Select the 100 most abundant ORFs in our dataset.
      mostAbundantORFnames = names(sort(rowSums(Hadza$orfs$tpm), decreasing=TRUE))[1:100]
      mostAbundantORFs = subsetORFs(Hadza, mostAbundantORFnames)
