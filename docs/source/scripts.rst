*************************************
Scripts, output files and file format
*************************************

.. note::
    Most of the information contained in the output files listed below can be more easily explored through :doc:`SQMtools`.


Upon startup, the pipeline will initially create the following files:

- ``<project>/creator.txt``: text file containing the SqueezeMeta version when the project was created
- ``<project>/SqueezeMeta_conf.pl``: run configuration. You will need to edit this file and :ref:`restart <restart>` if you want to change parameters mid-run
- ``<project>/parameters.pl``: additional parameters
- ``<project>/methods.txt``: text file containing names a citations of the software that is called during the run
- ``<project>/syslog``: text file containing the commands being called by the pipeline and their STDOUT/STDERR outputs, useful for debugging

Step 1: Assembly
================

**Script:** *01.run_assembly.pl*

Files produced
--------------
- ``<project>/results/01.<project>.fasta``: FASTA file containing the contigs resulting from the assembly 
- ``<project>/intermediate/01.<project>.lon``: Length of the contigs
- ``<project>/intermediate/01.<project>.stats``: Some statistics on the assembly (N50, N90, number of reads)

.. note::
  
  The merged/seqmerge modes will also produce a .fasta and a .lon file for each sample)


Step 2: RNA finding
===================

**Script:** *02.run_barrnap.pl*

Files produced
--------------
- ``<project>/results/02.<project>.rnas``: FASTA file containing all rRNAs and tRNAs found in the assembly
- ``<project>/results/02.<project>.16S``: Assignment (RDP classifier) for the 16S rRNAs sequences
- ``<project>/intermediate/02.<project>.maskedrna.fasta``: Fasta file containing the contigs resulting from the assembly, masking the positions where a rRNA/tRNA was found.

Step 3: Gene prediction
=======================

**Script:** *03.run_prodigal.pl*

Files produced
--------------
- ``<project>/results/03.<project>.fna``: Nucleotide sequences for predicted ORF
- ``<project>/results/03.<project>.faa``: Aminoacid sequences for predicted ORF
- ``<project>/results/03.<project>.gff``: Features and position in contigs for each of the predicted genes (this file will be moved to the ``intermediate`` directory if the ``-D`` option is selected)


Step 4: Homology searching against taxonomic (nr) and functional (COG, KEGG) databases
======================================================================================

**Script:** *04.rundiamond.pl*

Files produced
--------------
- ``<project>/intermediate/04.<project>.nr.diamond``: result of the homology search against the nr database
- ``<project>/intermediate/04.<project>.kegg.diamond``: result of the homology search against the KEGG database
- ``<project>/intermediate/04.<project>.eggnog.diamond``: result of the homology search against the eggNOG database
- ``<project>/intermediate/DB_BUILD_DATE``: date at which the SqueezeMeta database was originally created

.. note::

  If additional databases were provided using the ``-extdb`` option, this script will create additional diamond result files for each database

Step 5: HMM search for Pfam database
====================================

**Script:** *05.run_hmmer.pl*

Files produced
--------------
- ``<project>/intermediate/05.<project>.pfam.hmm``: results of the HMM search against the Pfam database


.. _lca script:
Step 6: Taxonomic assignment
============================

**Script:** *06.lca.pl*

Files produced
--------------
- ``<project>/results/06.<project>.fun3.tax.wranks``: taxonomic assignments for each ORF, including taxonomic ranks
- ``<project>/results/06.<project>.fun3.tax.noidfilter.wranks``: same as above, but the assignment is done without considering identity filters (see :ref:`lca`)

.. note::
  These files will be moved to the ``intermediate`` directory if the ``-D`` option is selected

.. _fun3 script:
Step 7: Functional assignment
=============================

**Script:** *07.fun3assign.pl*

Files produced
--------------
- ``<project>/results/07.<project>.fun3.cog``: PFAM functional assignment for each ORF
- ``<project>/results/07.<project>.fun3.kegg``: PFAM functional assignment for each ORF

Format of these files:

- Column 1: Name of the ORF
- Column 2: Best hit assignment
- Column 3: Best average assignment (see :ref:`fun3`)

.. note::
  - These files will be moved to the ``intermediate`` directory if the ``-D`` option is selected
  - If additional databases were provided using the ``-optdb`` option, this script will create additional result files for each database

- ``<project>/results/07.<project>.pfam``: PFAM functional assignment for each ORF

Step 8: Blastx on parts of the contigs without gene prediction or without hits
==============================================================================

**Script:** *08.blastx.pl*

This script will only be executed if the ``-D`` option was selected.

Files produced
--------------

- ``<project>/results/08.<project>.gff``: features and position in contigs for each of the Prodigal and BlastX ORFs Blastx 
- ``<project>/results/08.<project>.fun3.tax.wranks``: taxonomic assignment for the mix of Prodigal and BlastX ORFs, including taxonomic ranks
- ``<project>/results/08.<project>.fun3.tax.noidfilter.wranks``: same as above, but the assignment is done without considering identity filters (see :ref:`lca`)
- ``<project>/results/08.<project>.fun3.cog``: COG functional assignment for the mix of Prodigal and BlastX ORFs
- ``<project>/results/08.<project>.fun3.kegg``: KEGG functional assignment for the mix of Prodigal and BlastX ORFs 
- ``<project>/intermediate/blastx.fna``: nucleotide sequences for BlastX ORFs 

.. note::
  If additional databases were provided using the ``-optdb`` option, this script will create additional result files for each database

Step 9: Taxonomic assignment of contigs
=======================================

**Script:** *09.summarycontigs3.pl*

Files produced
--------------
- ``<project>/intermediate/09.<project>.contiglog``: consensus taxonomic assignment for the contigs (see :ref:`consensus tax`)

Format of the file:

- Column 1: name of the contig
- Column 2: taxonomic assignment, with ranks
- Column 3: lower rank of the assignment
- Column 4: disparity value (see :ref:`disparity`)
- Column 5: number of genes in the contig

.. _mappingstat:
Step 10: Mapping of reads to contigs and calculation of abundance measures
==========================================================================

**Script:** *10.mapsamples.pl*

Files produced
--------------
- ``<project>/results/10.<project>.mappingstat``:

- ``<project>/intermediate/10.<project>.mapcount``: several measurements regarding mapping of reads to ORFs

Format of the file:

    - Column 1: ORF name
    - Column 2: ORF length (nucleotides)
    - Column 3: number of reads mapped to that ORF
    - Column 4: number of bases mapped to that ORF
    - Column 5: RPKM value for the ORF
    - Column 6: coverage value for the ORF (Bases mapped / ORF length)
    - Column 7: TPM value for the ORF
    - Column 8: sample to which these abundance values correspond

- ``<project>/intermediate/10.<project>.contigcov``: several measurements regarding mapping of reads to contigs

Format of the file:

    - Column 1: ORF name
    - Column 2: coverage value for the contig
    - Column 3: RPKM value for the contig
    - Column 4: TPM value for the contig
    - Column 5: contig length (nucleotides)
    - Column 6: number of reads mapped to that contig
    - Column 7: number of bases mapped to that contig
    - Column 8: sample to which these abundance values correspond


Step 11: Calculation of the abundance of all taxa
=================================================

**Script:** *11.mcount.pl*

Files produced
--------------
- ``<project>/results/11.<project>.mcount``

Format of the file:

    - Column 1: taxonomic rank for the taxon
    - Column 2: taxon
    - Column 3: accumulated contig size: Sum of the length of all contigs for that taxon
    - Column 4 (and all even columns from this one): number of reads mapping to the taxon in the corresponding sample
    - Column 5 (and all odd columns from this one): number of bases mapping to the taxon in the corresponding sample

.. _funcover:
Step 12: Calculation of the abundance of all functions
======================================================

**Script:** *12.funcover.pl*

Files produced
--------------

- ``<project>/ext_tables/12.<project>.cog.stamp``: COG function table for `STAMP <http://kiwi.cs.dal.ca/Software/STAMP>`_

    - Column 1: functional class for the COG
    - Column 2: COG ID and function name
    - Column 3 and above: abundance of reads for that COG in the corresponding sample

- ``<project>/ext_tables/12.<project>.kegg.stamp``: KEGG function table for `STAMP <http://kiwi.cs.dal.ca/Software/STAMP>`_

    - Column 1: KEGG ID and function name
    - Column 2 and above: abundance of reads for that KEGG in the corresponding sample

- ``<project>/results/12.<project>.cog.funcover``: Several measurements of the abundance and distribution of each COG	

    - Column 1: COG ID
    - Column 2: sample name
    - Column 3: number of different ORFs of this function in the corresponding sample (copy number)
    - Column 4: sum of the length of all ORFs of this function in the corresponding sample (Total length)
    - Column 5: sum of the bases mapped to all ORFs of this function in the corresponding sample (Total bases)
    - Column 6: coverage of the function (Total bases / Total length)
    - Column 7: TPM value for the function
    - Column 9: number of the different taxa per rank (k: kingdom, p: phylum; c: class; o: order; f: family; g: genus; s: species) in which this COG has been found
    - Column 10: function of the COG

- ``<project>/results/12.<project>.kegg.funcover``: several measurements of the abundance and distribution of each KEGG. This has the same format as the ``cog.funcover`` file but replacing COGs by KEGGs. Additionally, the function of the KEGG will be present in column 11, while column 10 will contain the name of the KEGG

.. note::
  If additional databases were provided using the ``-extdb`` option, this script will create additional result files for each database

.. _ORF table:
Step 13: Creation of the ORF table
==================================

**Script:** *13.mergeannot2.pl*

Files produced
--------------
- ``<project>/results/13.<project>.orftable``
    - Column 1: ORF name
    - Column 2: Contig name
    - Column 3: molecule (CDS or RNA)
    - Column 4: method of ORF prediction (prodigal, barrnap, blastx)
    - Column 5: ORF length (nucleotides)
    - Column 6: ORF length (amino acids)
    - Column 7: GC percentage for the ORF
    - Column 8: Gene name
    - Column 9: Taxonomy for the ORF
    - Column 10: KEGG ID for the ORF (If a ``*`` sign is shown here, it means that the functional assignment was done by both best hit and best average scores, therefore is more reliable. Otherwise, the assignment was done using just the best hit, but there is evidence of a conflicting annotation)
    - Column 11: KEGG function
    - Column 12: KEGG functional class
    - Column 13: COG ID for the ORF (If a * sign is shown here, it means that the functional assignment was done by both best hit and best average scores, therefore is more reliable. Otherwise, the assignment was done using just the best hit, but there is evidence of a conflicting annotation)
    - Column 14: COG function
    - Column 15: COG functional class
    - Column 16: function in the external database provided
    - Column 17: Pfam annotation
    - Column 18 and beyond: TPM, coverage, read count and base count for the ORF in the different samples

.. note::                                                                                                                              If additional databases were provided using the ``-extdb`` option, functions and functional classes will be shown for each of them after column 15

Step 14: Binning
================

**Script:** *14.runbinning.pl*

Files produced
--------------

- ``<project>/intermediate/binners/maxbin``: directory containing fasta files with the contigs assigned to each bin by MaxBin (if selected)
- ``<project>/intermediate/binners/metabat``: directory containing fasta files with the contigs assigned to each bin by MetaBAT 2 (if selected)
- ``<project>/intermediate/binners/concoct``: directory containing fasta files with the contigs assigned to each bin by CONCOCT (if selected)


Step 15: Merging bins with DAS Tool
===================================

**Script:** *15.dastool.pl*

Files produced
--------------
- ``<project>/results/bins``: directory containing fasta files with the contigs associated to each bin after integrating the results for all binners with DAS Tool. If only one binner was selected, DAS Tool will not be run and the directory will instead contain the results for that binner


Step 16: Taxonomic assignment of bins
=====================================

**Script:** *16.addtax2.pl*

Files produced
--------------
- One taxonomy file for each fasta in the ``<project>/results/bins`` directory
- ``<project>/intermediate/16.<project>.bintax``: consensus taxonomic assignment for the bins (see :ref:`consensus tax`)
    - Column 1: binning method
    - Column 2: name of the bin
    - Column 3: taxonomic assignment for the bin, with ranks
    - Column 4: size of the bin (accumulated sum of contig lengths)
    - Column 5: disparity of the bin (see :ref:`disparity`)

.. note::
  Note that the taxonomy generated here is the consensus from the individual taxonomic assignments for each contig in the bin, not a GTDB-Tk taxonomy (which would be more precise). That can be achieved by adding the `--gtdbtk` flag, and is obtained during :ref:`bin annot` 

.. _bin annot:
Step 17: Running CheckM2 and optionally GTDB-Tk on bins
=======================================================

**Script:** *17.checkbins.pl*

Files produced
--------------
- ``<project>/intermediate/17.<project>.checkM``: Raw output from CheckM2
- If ``--gtdbtk`` is specified when running SqueezeMeta, also:
    - ``<project>/intermediate/17.<project>.gtdbtk``: GTDB-Tk output for archaeal and bacterial bins combined

Step 18: Creation of the bin table
==================================

**Script:** *18.getbins.pl*

Files produced
--------------

- ``<project>/intermediate/18.<project>.bincov``: coverage and TPM values for each bin
    - Column 1: bin name
    - Column 2: binning method
    - Column 3: coverage of the bin in the corresponding sample (Sum of bases from reads in the sample mapped to contigs in the bin / Sum of length of contigs in the bin)
    - Column 4: TPM for the bin in the corresponding sample (Sum of reads from the corresponding sample mapping to contigs in the bin x 10^6 /  Sum of length of contigs in the bin x Total number of reads)
    - Column 5: sample name

- ``<project>/results/18.<project>.bintable``: compilation of all data for bins
    - Column 1: bin name
    - Column 2: binning method
    - Column 3: taxonomic annotation (from the annotations of the contigs)
    - Column 4: taxonomy for the 16S rRNAs if the bin (if any)
    - Column 5: bin size (sum of length of the contigs)
    - Column 6: GC percentage for the bin
    - Column 7: number of contigs in the bin
    - Column 8: disparity of the bin
    - Column 9: completeness of the bin (CheckM2)
    - Column 10: contamination of the bin (CheckM2)
    - Column 11: strain heterogeneity of the bin (checkM)
    - Column 12 and beyond: coverage and TPM values for the bin in each sample.

.. note::                                                                                                                              If GTDB-Tk was run to classify the bins by adding the ``-gtdbtk`` option, an additional column named ``Tax GTDB-Tk`` will be present after column 4 in the file ``<project>/results/18.<project>.bintable``


Step 19: Creation of the contig table
=====================================

**Script:** *19.getcontigs.pl*

Files produced
--------------

- ``<project>/intermediate/19.<project>.contigsinbins``: list of contigs and corresponding bins

- ``<project>/results/19.<project>.contigtable``: compilation of data for contigs
    - Column 1: contig name
    - Column 2: taxonomic annotation for the contig (from the annotations of the ORFs)
    - Column 3: disparity of the contig
    - Column 4: GC percentage for the contig
    - Column 5: contig length
    - Column 6: number of genes in the contig
    - Column 7: bin to which the contig belong (if any)
    - Column 8 and beyond: values of coverage, TPM and number of mapped reads for the contig in each sample

Step 20: Prediction of pathway presence in bins using MinPath
=============================================================

**Script:** *20. minpath.pl*

Files produced
--------------

- ``<project>/results/20.<project>.kegg.pathways``: prediction of KEGG pathways in bins
    - Column 1: bin name
    - Column 2: taxonomic annotation for the bin
    - Column 3: number of KEGG pathways found
    - Column 4 and beyond: NF indicates that the pathway was not predicted. A number shows that the pathway was predicted to be present, and correspond to the number of enzymes of that pathway that were found.

- ``<project>/results/20.<project>.metacyc.pathways``: prediction of Metacyc pathways in bins. Format is similar as for the file above

.. _stats:
Step 21: Final statistics for the run
=====================================

**Script:** *21.stats.pl*

Files produced
--------------
- ``<project>/results/21.<project>.stats``: several statistics regarding ORFs, contigs and bins 

.. _sqm2tables in pipeline:
Step 22: Calculation of summary tables for the project
======================================================

**Script:** *sqm2tables.py*

Files produced
--------------
This script is executed with default parameters at the end of a SqueezeMeta run, and its results are placed in the ``<project>/results/tables`` directory. You may still want to run it on your own if you want to use non-default parameters. A list of output files can be found :ref:`here <sqm2tables output>`. This script is executed with default parameters at the end of a SqueezeMeta run, and its results are placed in the ``<project>/results/tables`` directory. You may still want to run it on your own if you want to use non-default parameters. A list of output files can be found :ref:`here <sqm2tables output>`.
