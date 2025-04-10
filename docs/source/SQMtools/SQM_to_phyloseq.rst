***************
SQM_to_phyloseq
***************

.. container::

   =============== ===============
   SQM_to_phyloseq R Documentation
   =============== ===============

   .. rubric:: Convert a SQM object into a phyloseq object from the
      *phyloseq* package
      :name: SQM_to_phyloseq

   .. rubric:: Description
      :name: description

   This function will convert the selected features from a SQM object
   into a phyloseq object from the
   `phyloseq <https://joey711.github.io/phyloseq/>`__ package. When
   possible, it will also include the taxonomy of the included features.
   Optionally, it accepts a meta table that will be passed as provided
   to the ``phyloseq`` object constructor.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      SQM_to_phyloseq(
        SQM,
        features = "genus",
        count = "abund",
        md = NULL,
        nocds = "treat_separately",
        no_partial_classifications = FALSE,
        ignore_unclassified = FALSE,
        ignore_unmapped = FALSE,
        bin_tax_source = "SQM",
        include_seqs = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +--------------------------------+------------------------------------+
   | ``SQM``                        | A SQM, SQMbunch or SQMlite object. |
   +--------------------------------+------------------------------------+
   | ``features``                   | character. Either ``"orfs"``,      |
   |                                | ``"contigs"``, ``"bins"``, any     |
   |                                | taxonomic rank included in         |
   |                                | ``SQM$taxa`` or any functional     |
   |                                | classication included in           |
   |                                | ``SQM$functions`` (default         |
   |                                | ``"tax"``). Note that a given      |
   |                                | feature type might not be          |
   |                                | available in this objects (e.g.    |
   |                                | ``"contigs"`` in SQMlite objects   |
   |                                | originating from a SQM reads       |
   |                                | project).                          |
   +--------------------------------+------------------------------------+
   | ``count``                      | character. Either ``"abund"`` for  |
   |                                | raw abundances, ``"percent"`` for  |
   |                                | percentages, ``"bases"`` for raw   |
   |                                | base counts, ``"cov"`` for         |
   |                                | coverages, code"cpm" for coverages |
   |                                | per million reads, ``"tpm"`` for   |
   |                                | TPM normalized values or           |
   |                                | ``"copy_number"`` for copy numbers |
   |                                | (default ``"abund"``). Note that a |
   |                                | given count type might not         |
   |                                | available in this object (e.g. TPM |
   |                                | or copy number in SQMlite objects  |
   |                                | originating from a SQM reads       |
   |                                | project).                          |
   +--------------------------------+------------------------------------+
   | ``md``                         | data.frame. A optional data.frame  |
   |                                | containing metadata for the        |
   |                                | samples in the SQM object.         |
   +--------------------------------+------------------------------------+
   | ``nocds``                      | character. Either                  |
   |                                | ``"treat_separately"`` to treat    |
   |                                | features annotated as No CDS       |
   |                                | separately,                        |
   |                                | ``"treat_as_unclassified"`` to     |
   |                                | treat them as Unclassified or      |
   |                                | ``"ignore"`` to ignore them in the |
   |                                | output (default                    |
   |                                | ``"treat_separately"``).           |
   +--------------------------------+------------------------------------+
   | ``no_partial_classifications`` | logical. When ``features`` is a    |
   |                                | taxonomic rank, treat features not |
   |                                | fully classified at the requested  |
   |                                | level (e.g. "Unclassified          |
   |                                | bacteroidota" at the class level   |
   |                                | or below) as fully unclassified.   |
   |                                | This takes effect before           |
   |                                | ``ignore_unclassified``, so if     |
   |                                | both are ``TRUE`` the plot will    |
   |                                | only contain features that were    |
   |                                | fully classified at the requested  |
   |                                | level (default ``FALSE``).         |
   +--------------------------------+------------------------------------+
   | ``ignore_unclassified``        | logical. When ``features`` is a    |
   |                                | taxonomic rank or functional       |
   |                                | category, don't include            |
   |                                | unclassified reads in the output   |
   |                                | (default ``FALSE``).               |
   +--------------------------------+------------------------------------+
   | ``ignore_unmapped``            | logical. Don't include unmapped    |
   |                                | reads in the output (default       |
   |                                | ``FALSE``).                        |
   +--------------------------------+------------------------------------+
   | ``bin_tax_source``             | character. Source of taxonomy when |
   |                                | ``features = "bins"``, either      |
   |                                | ``"SQM"`` of ``"gtdb"`` (default   |
   |                                | ``"gtdb"``).                       |
   +--------------------------------+------------------------------------+
   | ``include_seqs``               | logical. Whether to include        |
   |                                | sequences or not if creating a     |
   |                                | microtable from ORFs or contigs    |
   |                                | (default ``FALSE``).               |
   +--------------------------------+------------------------------------+

   .. rubric:: Value
      :name: value

   A phyloseq object.

   .. rubric:: See Also
      :name: see-also

   ``SQM_to_microeco`` for exporting a SQM/SQMlite/SQM object as a
   microtable object.
