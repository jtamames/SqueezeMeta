**********
subsetBins
**********

.. container::

   ========== ===============
   subsetBins R Documentation
   ========== ===============

   .. rubric:: Create a SQM object containing only the requested bins,
      and the contigs and ORFs contained in them.
      :name: subsetBins

   .. rubric:: Description
      :name: description

   Create a SQM object containing only the requested bins, and the
   contigs and ORFs contained in them.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetBins(
        SQM,
        bins = NULL,
        rank = NULL,
        tax = NULL,
        min_completeness = NULL,
        max_contamination = NULL,
        tax_source = "bins",
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = TRUE,
        rescale_copy_number = TRUE,
        allow_empty = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``SQM``                          | SQM object to be subsetted.      |
   +----------------------------------+----------------------------------+
   | ``bins``                         | character. Vector of bins to be  |
   |                                  | selected. If provided, will      |
   |                                  | override ``rank``, ``tax``,      |
   |                                  | ``min_completeness`` and         |
   |                                  | ``max_contamination``.           |
   +----------------------------------+----------------------------------+
   | ``rank``                         | character. The taxonomic rank    |
   |                                  | from which to select the desired |
   |                                  | taxa (``superkingdom``,          |
   |                                  | ``phylum``, ``class``,           |
   |                                  | ``order``, ``family``,           |
   |                                  | ``genus``, ``species``)          |
   +----------------------------------+----------------------------------+
   | ``tax``                          | character. A taxon or vector of  |
   |                                  | taxa to be selected.             |
   +----------------------------------+----------------------------------+
   | ``min_completeness``             | numeric. Discard bins with       |
   |                                  | completeness lower than this     |
   |                                  | value (default ``NULL``).        |
   +----------------------------------+----------------------------------+
   | ``max_contamination``            | numeric. Discard bins with       |
   |                                  | contamination higher than this   |
   |                                  | value (default ``NULL``).        |
   +----------------------------------+----------------------------------+
   | ``tax_source``                   | character, source data used for  |
   |                                  | taxonomic subsetting (if         |
   |                                  | ``rank`` and ``tax`` are         |
   |                                  | provided) and for the aggregate  |
   |                                  | taxonomy tables present in       |
   |                                  | ``SQM$taxa``, either ``"orfs"``, |
   |                                  | ``"contigs"``, ``"bins"`` (GTDB  |
   |                                  | bin taxonomy if available, SQM   |
   |                                  | bin taxonomy otherwise),         |
   |                                  | ``"bins_gtdb"`` (GTDB bin        |
   |                                  | taxonomy) or ``"bins_sqm"`` (SQM |
   |                                  | bin taxonomy). If using          |
   |                                  | ``bins_gtdb``, note that GTDB    |
   |                                  | taxonomy may differ from the     |
   |                                  | NCBI taxonomy used throughout    |
   |                                  | the rest of SqueezeMeta. Default |
   |                                  | ``"bins"``.                      |
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
   |                                  | parent object. By default it is  |
   |                                  | set to ``TRUE``, which means     |
   |                                  | that the returned TPMs will be   |
   |                                  | scaled *by million of reads of   |
   |                                  | the selected bins*.              |
   +----------------------------------+----------------------------------+
   | ``rescale_copy_number``          | logical. If ``TRUE``, copy       |
   |                                  | numbers with be recalculated     |
   |                                  | using the median single-copy     |
   |                                  | gene coverages in the subset.    |
   |                                  | Otherwise, single-copy gene      |
   |                                  | coverages will be taken from the |
   |                                  | parent object. By default it is  |
   |                                  | set to ``TRUE``, which means     |
   |                                  | that the returned copy numbers   |
   |                                  | for each function will represent |
   |                                  | the average copy number of that  |
   |                                  | function *per genome of the      |
   |                                  | selected taxon*.                 |
   +----------------------------------+----------------------------------+
   | ``allow_empty``                  | (internal use only).             |
   +----------------------------------+----------------------------------+

   .. rubric:: Value
      :name: value

   SQM object containing only the requested bins.

   .. rubric:: See Also
      :name: see-also

   ``subsetContigs``, ``subsetORFs``

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Which are the most complete bins?
      topBinNames = rownames(Hadza$bins$table)[order(Hadza$bins$table[,"Completeness"],
                                               decreasing=TRUE)][1:2]
      # Subset with the most complete bin.
      topBin = subsetBins(Hadza, topBinNames[1])

      # Subset with all the bins over 90% completeness
      over90 = subsetBins(Hadza, min_completeness = 90)

      # Subset with bins from the Phascolarctobacterium genus using SqueezeMeta's taxonomy
      phasco = subsetBins(Hadza, tax_source = "bins", rank = "genus", tax = "Phascolarctobacterium")

      # Subset with binsfrom the Bacteroidota phylum using GTDB taxonomy
      bact = subsetBins(Hadza, tax_source = "bins_gtdb", rank = "phylum", tax = "p__Bacteroidota")
