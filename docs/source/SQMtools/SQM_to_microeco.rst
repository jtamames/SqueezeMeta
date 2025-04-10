***************
SQM_to_microeco
***************

.. container::

   =============== ===============
   SQM_to_microeco R Documentation
   =============== ===============

   .. rubric:: Convert a SQM object into a microtable object from the
      *microeco* package
      :name: SQM_to_microeco

   .. rubric:: Description
      :name: description

   This function will convert the selected features from a SQM object
   into an object of the ``microtable`` class from the
   `microeco <https://chiliubio.github.io/microeco>`__ package. When
   possible, it will also include the taxonomy of the included features.
   Optionally, it accepts a meta table that will be passed as provided
   to ``microtable$new``.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      SQM_to_microeco(
        SQM,
        features = "genus",
        count = "abund",
        md = NULL,
        bin_tax_source = "SQM",
        include_seqs = FALSE
      )

   .. rubric:: Arguments
      :name: arguments

   +--------------------+------------------------------------------------+
   | ``SQM``            | A SQM, SQMbunch or SQMlite object.             |
   +--------------------+------------------------------------------------+
   | ``features``       | character. Either ``"orfs"``, ``"contigs"``,   |
   |                    | ``"bins"``, any taxonomic rank included in     |
   |                    | ``SQM$taxa`` or any functional classication    |
   |                    | included in ``SQM$functions`` (default         |
   |                    | ``"tax"``). Note that a given feature type     |
   |                    | might not be available in this objects (e.g.   |
   |                    | ``"contigs"`` in SQMlite objects originating   |
   |                    | from a SQM reads project).                     |
   +--------------------+------------------------------------------------+
   | ``count``          | character. Either ``"abund"`` for raw          |
   |                    | abundances, ``"percent"`` for percentages,     |
   |                    | ``"bases"`` for raw base counts, ``"cov"`` for |
   |                    | coverages, code"cpm" for coverages per million |
   |                    | reads, ``"tpm"`` for TPM normalized values or  |
   |                    | ``"copy_number"`` for copy numbers (default    |
   |                    | ``"abund"``). Note that a given count type     |
   |                    | might not available in this object (e.g. TPM   |
   |                    | or copy number in SQMlite objects originating  |
   |                    | from a SQM reads project).                     |
   +--------------------+------------------------------------------------+
   | ``md``             | data.frame. A optional data.frame containing   |
   |                    | metadata for the samples in the SQM object.    |
   +--------------------+------------------------------------------------+
   | ``bin_tax_source`` | character. Source of taxonomy when             |
   |                    | ``features = "bins"``, either ``"SQM"`` of     |
   |                    | ``"gtdb"`` (default ``"gtdb"``).               |
   +--------------------+------------------------------------------------+
   | ``include_seqs``   | logical. Whether to include sequences or not   |
   |                    | if creating a microtable from contigs (default |
   |                    | ``FALSE``).                                    |
   +--------------------+------------------------------------------------+

   .. rubric:: Value
      :name: value

   A ``microtable``.

   .. rubric:: See Also
      :name: see-also

   ``SQM_to_phyloseq`` for exporting a SQM/SQMlite/SQM object as a
   phyloseq object.
