*****
MGKOs
*****

===== ===============
MGKOs R Documentation
===== ===============

Single Copy Phylogenetic Marker Genes from Sunagawa's group (KOs)
-----------------------------------------------------------------

Description
~~~~~~~~~~~

Lists of Single Copy Phylogenetic Marker Genes. These are useful for
transforming coverages or tpms into copy numbers. This is an alternative
way of normalizing data in order to be able to compare functional
profiles in samples with different sequencing depths.

Usage
~~~~~

.. code:: R

   data(MGKOs)

Format
~~~~~~

Character vector with the KEGG identifiers for 10 Single Copy
Phylogenetic Marker Genes.

References
~~~~~~~~~~

Salazar, G *et al.* (2019). Gene Expression Changes and Community
Turnover Differentially Shape the Global Ocean Metatranscriptome *Cell*
**179**:1068-1083.
(`PubMed <https://pubmed.ncbi.nlm.nih.gov/31730850/>`__).

See Also
~~~~~~~~

``MGOGs`` for an equivalent list using OGs instead of KOs; ``USiCGs``
for an alternative set of single copy genes, and for examples on how to
generate copy numbers.
