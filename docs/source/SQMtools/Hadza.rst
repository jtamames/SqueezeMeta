*****
Hadza
*****

.. container::

   ===== ===============
   Hadza R Documentation
   ===== ===============

   .. rubric:: Hadza hunter-gatherer gut metagenomes
      :name: Hadza

   .. rubric:: Description
      :name: description

   Subset of two bins (and the associated contigs and genes) generated
   by running SqueezeMeta on two gut metagenomic samples obtained from
   two hunter-gatherers of the Hadza ethnic group.

   .. rubric:: Usage
      :name: usage

   .. code:: R

      data(Hadza)

   .. rubric:: Format
      :name: format

   A SQM object; see ``loadSQM``.

   .. rubric:: Source
      :name: source

   `SRR1927149 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR1927149>`__,
   `SRR1929485 <https://www.ncbi.nlm.nih.gov/sra/?term=SRR1929485>`__.

   .. rubric:: References
      :name: references

   Rampelli *et al.*, 2015. Metagenome Sequencing of the Hadza
   Hunter-Gatherer Gut Microbiota. *Curr. biol.* **25**:1682-93
   (`PubMed <https://pubmed.ncbi.nlm.nih.gov/25981789/>`__).

   .. rubric:: Examples
      :name: examples

   .. code:: R

      data(Hadza)
      plotTaxonomy(Hadza, "genus", rescale=TRUE)
      plotFunctions(Hadza, "COG")
