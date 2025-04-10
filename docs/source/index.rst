.. image:: ../resources/logo.svg
  :width: 20%
  :align: right
  :alt: SqueezeMeta logo

Welcome to SqueezeMeta's documentation!
=======================================

**SqueezeMeta** is a fully automatic pipeline for
metagenomics/metatranscriptomics, covering all steps of the analysis.
SqueezeMeta includes multi-metagenome support allowing the co-assembly
of related metagenomes and the retrieval of individual metagenome-assembled genomes (MAGs)
via binning procedures. Thus, SqueezeMeta features several characteristics:

1) Several :ref:`assembly and co-assembly algorithms and strategies <Assembly Strategy>` for short and long reads
2) Several binning algorithms for the recovery of metagenome-assembled genomes (MAGs)
3) Taxonomic annotation, functional annotation and quantification of genes, contigs, and bins
4) Support for the annotation and quantification of :ref:`pre-existing assemblies or collections of genomes <extassembly>`
5) Support for :ref:`de-novo metatranscriptome assembly <metatranscriptomics>` and :ref:`hybrid metagenomics/metatranscriptomics projects <metag metat>`
6) Support for the :ref:`annotation of unassembled shotgun metagenomic reads <alt modes short>`
7) An :doc:`R package <SQMtools>` to easily explore your results, including bindings for `microeco <https://chiliubio.github.io/microeco/>`_ and `phyloseq <https://joey711.github.io/phyloseq/>`_

.. note::
  Check out the :doc:`use_cases` section for more information.

SqueezeMeta uses a combination of custom scripts and external
software packages for the different steps of the analysis:

1)  Assembly
2)  RNA prediction and classification
3)  ORF (CDS) prediction
4)  Homology searching against taxonomic and functional databases
5)  Hmmer searching against Pfam database
6)  Taxonomic assignment of genes
7)  Functional assignment of genes (OPTIONAL)
8)  Blastx on parts of the contigs with no gene prediction or no hits
9)  Taxonomic assignment of contigs, and check for taxonomic disparities
10) Coverage and abundance estimation for genes and contigs
11) Estimation of taxa abundances
12) Estimation of function abundances
13) Merging of previous results to obtain the ORF table
14) Binning with different methods
15) Binning integration with DAS tool
16) Taxonomic assignment of bins, and check for taxonomic disparities
17) Checking of bins with CheckM2 (and optionally classify them with
    GTDB-Tk)
18) Merging of previous results to obtain the bin table
19) Merging of previous results to obtain the contig table
20) Prediction of kegg and metacyc patwhays for each bin
21) Final statistics for the run
22) Generation of tables with aggregated taxonomic and functional
    profiles

Detailed information about the different steps of the pipeline can be
found in the :doc:`scripts` section.


Contents
--------

.. toctree::
   :maxdepth: 2

   use_cases
   installation
   execution
   adv_annotation
   scripts
   alt_modes
   SQMtools
   new_binners
   utils
   alg_details
