.. CONCOCT documentation master file, created by
   sphinx-quickstart on Thu Aug  1 11:22:50 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CONCOCT's documentation
=======================
CONCOCT "bins" metagenomic contigs. Metagenomic binning is the process of clustering sequences into clusters corresponding to operational taxonomic units of some level.

Features
--------
CONCOCT does unsupervised binning of metagenomic contigs by using nucleotide composition - kmer frequencies - and coverage data for multiple samples.
CONCOCT can accurately (up to species level) bin metagenomic contigs. For optimal performance:

- Map several samples against your assembled contigs.
- Cut longer contigs into 10 - 20 kb pieces prior to mapping.
- Evaluate your bins using single copy genes.

Installation
------------
For a comprehensive guide on how to install CONCOCT and all its dependencies, see :doc:`installation`.

Contribute
----------

- Issue Tracker: `github <https://github.com/BinPro/CONCOCT/issues>`__
- Source Code: `github <https://github.com/BinPro/CONCOCT>`__

Support
-------
If you are having issues, please let us know.
We have a discussion thread on gitter:

.. image:: https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg?style=flat-square
    :alt: Join the chat at gitter.im/BinPro/CONCOCT
    :target: https://gitter.im/BinPro/CONCOCT

Known Issues
------------

 - Contig names consisting of digits only are not allowed. Please rename your contigs in both the fasta and the coverage table before proceeding.
 - Contig sequences can only contain letters A,C,G or T. For example Ns are currently not allowed.
 - Contigs need to be cut up prior to binning. This is covered in the :doc:`usage` page.

For a more up to date list of reported issues, check the issue tracker:
https://github.com/BinPro/CONCOCT/issues

Licence
-------
FreeBSD

Contents:
---------

.. toctree::
   :maxdepth: 2

   self
   installation
   usage
   cmd_options
   scripts/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
