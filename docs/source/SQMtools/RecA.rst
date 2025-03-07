****
RecA
****

.. container::

   ==== ===============
   RecA R Documentation
   ==== ===============

   .. rubric:: RecA/RadA recombinase
      :name: RecA

   .. rubric:: Description
      :name: description

   The recombination protein RecA/RadA is essential for the repair and
   maintenance of DNA, and has homologs in every bacteria and archaea.
   By dividing the coverage of functions by the coverage of RecA,
   abundances can be transformed into copy numbers, which can be used to
   compare functional profiles in samples with different sequencing
   depths. RecA-derived copy numbers are available in the SQM object
   (``SQM$functions$<annotation_type>$copy_number``).

   .. rubric:: Usage
      :name: usage

   .. code:: R

      data(RecA)

   .. rubric:: Format
      :name: format

   Character vector with the COG identifier for RecA/RadA.

   .. rubric:: Source
      :name: source

   `EggNOG
   Database <http://eggnogdb.embl.de/#/app/results?seqid=P0A7G6&target_nogs=COG0468#COG0468_datamenu>`__.

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      data(RecA)
      ### Let's calculate the average copy number of each function in our samples.
      # We do it for COG annotations here, but we could also do it for KEGG or PFAMs.
      COG.coverage = Hadza$functions$COG$cov
      COG.copynumber = t(t(COG.coverage) / COG.coverage[RecA,]) # Sample-wise division by RecA coverage.
