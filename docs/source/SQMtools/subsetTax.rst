*********
subsetTax
*********

.. container::

   ========= ===============
   subsetTax R Documentation
   ========= ===============

   .. rubric:: Filter results by taxonomy
      :name: subsetTax

   .. rubric:: Description
      :name: description

   Create a SQM or SQMbunch object containing only the contigs with a
   given consensus taxonomy, the ORFs contained in them and the bins
   that contain them.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      subsetTax(
        SQM,
        rank,
        tax,
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = TRUE,
        rescale_copy_number = TRUE,
        recalculate_bin_stats = TRUE,
        allow_empty = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``SQM``                          | SQM object to be subsetted.      |
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
   |                                  | the selected taxon*.             |
   +----------------------------------+----------------------------------+
   | ``rescale_copy_number``          | logical. If ``TRUE``, copy       |
   |                                  | numbers with be recalculated     |
   |                                  | using the RecA/RadA coverages in |
   |                                  | the subset. Otherwise, RecA/RadA |
   |                                  | coverages will be taken from the |
   |                                  | parent object. By default it is  |
   |                                  | set to ``TRUE``, which means     |
   |                                  | that the returned copy numbers   |
   |                                  | for each function will represent |
   |                                  | the average copy number of that  |
   |                                  | function *per genome of the      |
   |                                  | selected taxon*.                 |
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

   SQM or SQMbunch object containing only the requested taxon.

   .. rubric:: See Also
      :name: see-also

   ``subsetFun``, ``subsetContigs``, ``subsetSamples``, ``combineSQM``.
   The most abundant items of a particular table contained in a SQM
   object can be selected with ``mostAbundant``.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      Hadza.Prevotella = subsetTax(Hadza, "genus", "Prevotella")
      Hadza.Bacteroidota = subsetTax(Hadza, "phylum", "Bacteroidota")
