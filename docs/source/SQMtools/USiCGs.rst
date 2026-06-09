******
USiCGs
******

====== ===============
USiCGs R Documentation
====== ===============

Universal Single-Copy Genes
---------------------------

Description
~~~~~~~~~~~

Lists of Universal Single Copy Genes for Bacteria and Archaea. These are
useful for transforming coverages or tpms into copy numbers. This is an
alternative way of normalizing data in order to be able to compare
functional profiles in samples with different sequencing depths.

Usage
~~~~~

.. code:: R

   data(USiCGs)

Format
~~~~~~

Character vector with the KEGG identifiers for 15 Universal Single Copy
Genes.

Source
~~~~~~

`Carr et al., 2013. Table
S1 <https://pubmed.ncbi.nlm.nih.gov/24146609/>`__.

References
~~~~~~~~~~

Carr, Shen-Orr & Borenstein (2013). Reconstructing the Genomic Content
of Microbiome Taxa through Shotgun Metagenomic Deconvolution *PLoS
Comput. Biol.* **9**:e1003292.
(`PubMed <https://pubmed.ncbi.nlm.nih.gov/24146609/>`__).

See Also
~~~~~~~~

``MGOGs`` and ``MGKOs`` for an alternative set of single copy genes, and
for examples on how to generate copy numbers.

Examples
~~~~~~~~

.. code:: R

   data(Hadza)
   data(USiCGs)
   ### Let's look at the Universal Single Copy Gene distribution in our samples.
   KEGG.tpm = Hadza$functions$KEGG$tpm
   all(USiCGs %in% rownames(KEGG.tpm)) # Are all the USiCGs present in our dataset?
   # Plot a boxplot of USiCGs tpms and calculate median USiCGs tpm.
   # This looks weird in the test dataset because it contains only a small subset of the metagenomes.
   # In a set of complete metagenomes USiCGs should have fairly similar TPM averages
   # and low dispersion across samples.
   boxplot(t(KEGG.tpm[USiCGs,]), names=USiCGs, ylab="TPM", col="slateblue2")
    
   ### Now let's calculate the average copy numbers of each function.
   # We do it for KEGG annotations here, but we could also do it for COGs or PFAMs.
   USiCGs.cov = apply(Hadza$functions$KEGG$cov[USiCGs,], 2, median)
   # Sample-wise division by the median USiCG coverage.
   KEGG.copynumber = t(t(Hadza$functions$KEGG$cov) / USiCGs.cov)
