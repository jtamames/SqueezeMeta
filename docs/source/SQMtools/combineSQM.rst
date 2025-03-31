**********
combineSQM
**********

.. container::

   ========== ===============
   combineSQM R Documentation
   ========== ===============

   .. rubric:: Combine several SQM objects
      :name: combineSQM

   .. rubric:: Description
      :name: description

   Combine an arbitrary number of SQM objects into a single SQM object
   (if the input objects contain the same samples, i.e. they come from
   the same SqueezeMeta run) or a single SQMbunch object. For combining
   results from sqm_reads.pl or sqm_longreads.pl please check
   ``combineSQMlite``. The parameters below (other than ...) will take
   only effect if the input objects contain the same samples. Otherwise
   the input objects will be taken as they are, with no recalculation of
   taxonomy, function or rescaling,

   .. rubric:: Usage
      :name: usage

   .. code:: R

      combineSQM(
        ...,
        tax_source = "orfs",
        trusted_functions_only = FALSE,
        ignore_unclassified_functions = FALSE,
        rescale_tpm = TRUE,
        rescale_copy_number = TRUE,
        recalculate_bin_stats = TRUE
      )

   .. rubric:: Arguments
      :name: arguments

   +----------------------------------+----------------------------------+
   | ``...``                          | an arbitrary number of SQM       |
   |                                  | objects. Alternatively, a single |
   |                                  | list containing an arbitrary     |
   |                                  | number of SQM objects.           |
   +----------------------------------+----------------------------------+
   | ``tax_source``                   | character. Features used for     |
   |                                  | calculating aggregated           |
   |                                  | abundances at the different      |
   |                                  | taxonomic ranks. Either          |
   |                                  | ``"orfs"`` or ``"contigs"``      |
   |                                  | (default ``"orfs"``). If the     |
   |                                  | objects being combined contain a |
   |                                  | subset of taxa or bins, this     |
   |                                  | parameter can be set to          |
   |                                  | ``TRUE``.                        |
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
   |                                  | ``TRUE``).                       |
   +----------------------------------+----------------------------------+
   | ``rescale_copy_number``          | logical. If ``TRUE``, copy       |
   |                                  | numbers with be recalculated     |
   |                                  | using the RecA/RadA coverages in |
   |                                  | the subset. Otherwise, RecA/RadA |
   |                                  | coverages will be taken from the |
   |                                  | parent object with the highest   |
   |                                  | RecA/RadA coverages. By default  |
   |                                  | it is set to ``TRUE``, which     |
   |                                  | means that the returned copy     |
   |                                  | numbers will represent the       |
   |                                  | average copy number per function |
   |                                  | *in the genomes of the selected  |
   |                                  | bins or contigs*. If any SQM     |
   |                                  | objects that are being combined  |
   |                                  | contain a functional subset      |
   |                                  | rather than a contig/bins        |
   |                                  | subset, this parameter should be |
   |                                  | set to ``FALSE``.                |
   +----------------------------------+----------------------------------+
   | ``recalculate_bin_stats``        | logical. If ``TRUE``, bin stats  |
   |                                  | and taxonomy are recalculated    |
   |                                  | based on the contigs present in  |
   |                                  | the subsetted object (default    |
   |                                  | ``TRUE``).                       |
   +----------------------------------+----------------------------------+

   .. rubric:: Value
      :name: value

   A SQM or SQMbunch object

   .. rubric:: See Also
      :name: see-also

   ``subsetFun``, ``subsetTax``, ``combineSQMlite``

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      # Select Carbohydrate metabolism ORFs in Bacteroidota,
      #  and Amino acid metabolism ORFs in Proteobacteria
      bact = subsetTax(Hadza, "phylum", "Bacteroidota")
      bact.carb = subsetFun(bact, "Carbohydrate metabolism")
      baci = subsetTax(Hadza, "phylum", "Bacillota")
      baci.amins = subsetFun(baci, "Amino acid metabolism")
      bact.carb_proteo.amins = combineSQM(bact.carb, baci.amins, rescale_copy_number=FALSE)
