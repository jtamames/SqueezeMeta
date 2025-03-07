*****
MGKOs
*****

.. container::

   ===== ===============
   MGKOs R Documentation
   ===== ===============

   .. rubric:: Single Copy Phylogenetic Marker Genes from Sunagawa's
      group (KOs)
      :name: MGKOs

   .. rubric:: Description
      :name: description

   Lists of Single Copy Phylogenetic Marker Genes. These are useful for
   transforming coverages or tpms into copy numbers. This is an
   alternative way of normalizing data in order to be able to compare
   functional profiles in samples with different sequencing depths.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      data(MGKOs)

   .. rubric:: Format
      :name: format

   Character vector with the KEGG identifiers for 10 Single Copy
   Phylogenetic Marker Genes.

   .. rubric:: References
      :name: references

   Salazar, G *et al.* (2019). Gene Expression Changes and Community
   Turnover Differentially Shape the Global Ocean Metatranscriptome
   *Cell* **179**:1068-1083.
   (`PubMed <https://pubmed.ncbi.nlm.nih.gov/31730850/>`__).

   .. rubric:: See Also
      :name: see-also

   ``MGOGs`` for an equivalent list using OGs instead of KOs; ``USiCGs``
   for an alternative set of single copy genes, and for examples on how
   to generate copy numbers.
