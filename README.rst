.. image:: https://github.com/jtamames/SqueezeMeta/blob/images/logo.svg
  :width: 20%
  :align: right
  :alt: SqueezeMeta logo

************************************************************************
SqueezeMeta: a fully automated metagenomics pipeline, from reads to bins
************************************************************************

**SqueezeMeta** is a fully automatic pipeline for
metagenomics/metatranscriptomics, covering all steps of the analysis.
SqueezeMeta includes multi-metagenome support allowing the co-assembly
of related metagenomes and the retrieval of individual metagenome-assembled genomes (MAGs)
via binning procedures. Thus, SqueezeMeta features several characteristics:

1) Several assembly and co-assembly algorithms and strategies for short and long reads
2) Several binning algorithms for the recovery of metagenome-assembled genomes (MAGs)
3) Taxonomic annotation, functional annotation and quantification of genes, contigs, and bins
4) Support for the annotation and quantification of pre-existing assemblies or collections of genomes
5) Support for de-novo metatranscriptome assembly and hybrid metagenomics/metatranscriptomics projects
6) Support for the annotation of unassembled shotgun metagenomic reads
7) An R package to easily explore your results, including bindings for `microeco <https://chiliubio.github.io/microeco/>`_ and `phyloseq <https://joey711.github.io/phyloseq/>`_

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
found in the `documentation <https://squeezemeta.readthedocs.io>`_.

*************
Documentation
*************

- The documentation for SqueezeMeta and SQMtools is available in `ReadTheDocs <https://squeezemeta.readthedocs.io>`_.
- The `wiki <https://github.com/jtamames/SqueezeMeta/wiki>`_ contains extra examples on how to use certain features of SqueezeMeta/SQMtools.
- You can also check the SqueezeMeta paper `here <https://www.frontiersin.org/articles/10.3389/fmicb.2018.03349/full>`_ and a second paper on how to analyse the output of SqueezeMeta `here <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03703-2>`_


************
Installation
************

SqueezeMeta is intended to be run in a x86-64 Linux OS (tested in Ubuntu
and CentOS). The easiest way to install it is by using conda. The
default conda solver might however be slow solving the dependencies, so
it’s better to first set up the libmamba solver with

::

   conda update -n base conda # if your conda version is below 22.11
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba

and then use conda to install SqueezeMeta

``conda create -n SqueezeMeta -c conda-forge -c bioconda -c fpusan squeezemeta --no-channel-priority --override-channels``

If you change ``squeezemeta`` to ``squeezemeta-dev`` you will instead
get the latest development version. This will contain additional bugfixes
and features, but potentially also new bugs, as it will not have been
tested as thoroughly as the stable version.

(If the environment does not solve and you get a message saying that
``__cuda`` is missing in your system, try adding ``CONDA_OVERRIDE_CUDA=12.4``
before the installation command: ``CONDA_OVERRIDE_CUDA=12.4 conda create ...``)

The command above will create a new conda environment named SqueezeMeta,
which must then be activated.

``conda activate SqueezeMeta``

When using conda, all the scripts from the SqueezeMeta distribution will
be available on ``$PATH``.

Alternatively, you can download the latest release from the GitHub
repository and uncompress the tarball in a suitable directory. The
tarball includes the SqueezeMeta scripts as well as the `third-party
software <https://squeezemeta.readthedocs.io/en/stable/installation.html#vendored-tools>`_
redistributed with SqueezeMeta. Note that, you may need to provide
additional dependencies, and potentially recompile some
binaries from source in order for the manual install to work.
The conda method is now the recommended way to install SqueezeMeta,
and we will not prioritize support to issues regarding manual installation.

The ``test_install.pl`` script can be run in order to check whether the
required dependencies are available in your environment.

``/path/to/SqueezeMeta/utils/install_utils/test_install.pl``

Downloading or building databases
=================================

SqueezeMeta uses several databases. GenBank nr for taxonomic assignment,
and eggnog, KEGG and Pfam for functional assignment. The script
*download_databases.pl* can be run to download a pre-formatted version
of all the databases required by SqueezeMeta.

``/path/to/SqueezeMeta/utils/install_utils/download_databases.pl /download/path``

, where ``/download/path`` is the destination folder. This is the
recommended option, but the files are hosted in our institutional
server, which can at times be unreachable.

Alternatively, the script ``make_databases.pl`` can be run to download
from source and format the latest version of the databases.

``/path/to/SqueezeMeta/utils/install_utils/make_databases.pl /download/path``

Generally, ``download_databases.pl`` is the safest choice for getting
your databases set up. When running ``make_databases.pl``, data download
(e.g. from the NCBI server) can be interrupted, leading to a corrupted
database. Always run ``test_install.pl`` to check that the database was
properly created. Otherwise, you can try re-running
``make_databases.pl``, or just run ``download_databases.pl`` instead.

The databases occupy 470Gb, but we recommend having at least 700Gb free
disk space during the building process.

Two directories will be generated after running either
``make_databases.pl`` or ``download_databases.pl``.

- ``/download/path/db``, which contains the actua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l databases.
- ``/download/path/test``, which contains data for a test run of SqueezeMeta.

If the SqueezeMeta databases are already built in another location in
the system, a different copy of SqueezeMeta can be configured to use
them with

``/path/to/SqueezeMeta/utils/install_utils/configure_nodb.pl /path/to/db``

, where ``/path/to/db`` is the route to the ``db`` folder that was
generated by either ``make_databases.pl`` or ``download_databases.pl``.

After configuring the databases, the ``test_install.pl`` can be run in
order to check that SqueezeMeta is ready to work (see previous section).

Testing SqueezeMeta
===================

The *download_databases.pl* and *make_databases.pl* scripts also
download two datasets for testing that the program is running correctly.
Assuming either was run with the directory ``/download/path`` as its
target the test run can be executed with

| ``cd </download/path/test>``
| ``SqueezeMeta.pl -m coassembly -p Hadza -s test.mock.samples -f raw``

Alternatively, ``-m sequential`` or ``-m merged`` can be used.

In addition to this mock dataset, we also provide two real metagenomes.
A test run on those can be executed with

``SqueezeMeta.pl -m coassembly -p Hadza -s test.samples -f raw``

Updating SqueezeMeta
====================

Assuming your databases are not inside the SqueezeMeta directory, just
remove it, download the new version and configure it with

``/path/to/SqueezeMeta/utils/install_utils/configure_nodb.pl /path/to/db``

********************
Usage considerations
********************

Choosing an assembly strategy
=============================

SqueezeMeta can be run in four different modes, depending of the type of
multi-metagenome support. These modes are:

-  **Sequential mode**: All samples are treated individually and analysed
   sequentially.

-  **Coassembly mode**: Reads from all samples are pooled and a single
   assembly is performed. Then reads from individual samples are mapped
   to the coassembly to obtain gene abundances in each sample. Binning
   methods allow to obtain genome bins.

-  **Merged mode**: if many big samples are available, co-assembly could
   crash because of memory requirements. This mode achieves a comparable
   resul with a procedure inspired by `the one used by Benjamin Tully for
   analysing TARA Oceans data <https://dx.doi.org/10.17504/protocols.io.hfqb3mw>`_.
   Briefly, samples are assembled individually and the resulting contigs are
   merged in a single co-assembly. Then the analysis proceeds as in the
   co-assembly mode. This is not the recommended procedure (use
   co-assembly if possible) since the possibility of creating chimeric
   contigs is higher. But it is a viable alternative in smaller computers in
   which standard co-assembly is not feasible.

-  **Seqmerge mode**: This is intended to work with more samples than the
   merged mode. Instead of merging all individual assemblies in a single
   step, which can be very computationally demanding, seqmerge works
   sequentially. First, it assembles individually all samples, as in
   merged mode. But then it will merge the two most similar assemblies.
   Similarity is measured as Amino Acid Identity values using the
   wonderful CompareM software by Donovan Parks. After this first
   merging, it again evaluates similarity and merge, and proceeds this
   way until all metagenomes have been merged in one. Therefore, for n
   metagenomes, it will need n-1 merging steps.

Note that the *merged* and *seqmerge* modes work well as a substitute of
coassembly for running small datasets in computers with low memory
(e.g. 16 Gb) but are very slow for analising large datasets (>10
samples) even in workstations with plenty of resources. Still, setting
``-contiglen`` to 1000 or higher can make *seqmerge* a viable strategy
even in those cases. Otherwise, we recommend to use either the
sequential or the co-assembly modes.

Regarding the choice of assembler, MEGAHIT and SPAdes work better with
short Illumina reads, while Canu and Flye support long reads from PacBio
or ONT-Minion. MEGAHIT (the default in SqueezeMeta) is more
resource-efficient than SPAdes, consuming less memory, but SPAdes
supports more analysis modes and produces slightly better assembly
statistics. SqueezeMeta can call SPAdes in three different ways. The
option ``-a spades`` is meant for metagenomic datasets, and will
automatically add the flags ``–meta -k 21,33,55,77,99,127`` to the
*spades.py* call. Conversely, ``-a rnaspades`` will add the flags
``–rna -k 21,33,55,77,99,127``. Finally, the option ``-a spades_base`` will add no
additional flags to the *spades.py* call. This can be used in
conjunction with ``–assembly options`` when one wants to fully customize
the call to SPAdes, e.g. for assembling single cell genomes.

Analizing user-supplied assemblies or bins
==========================================

An user-supplied assembly can be passed to SqueezeMeta with the flag
``-extassembly <your_assembly.fasta>``. The contigs in that fasta file
will be analyzed by the SqueezeMeta pipeline starting from step 2.
With this, you will be able to annotate your assembly, estimate its
abundance in your metagenomes/metatranscriptomes, and perform binning on it.

Additionally, a set of pre-existing genomes and bins can be passed to
SqueezeMeta with the flag ``-extbins <path_to_dir_with_bins>``. This will
work similarly to ``-extassembly``, but SqueezeMeta will treat each fasta
file in the input directory as an individual bin.

Using external databases for functional annotation
=====================================================

In addition to the databases distributed with SqueezeMeta, one or several user-provided databases
can be used for functional annotation. This is invoked using the ``-extdb`` option. Please refer to the 
`documentation <https://squeezemeta.readthedocs.io/en/stable/adv_annotation.html#using-external-databases-for-functional-annotation>`_
for details.

Working with Oxford Nanopore MinION and PacBio reads
====================================================

SqueezeMeta is able to work seamlessly with
single-end reads. In order to obtain better mappings of MinION and
PacBio reads against the assembly, we advise to use minimap2 for read
counting, by including the ``-map minimap2-ont`` (MinION) or ``-map minimap2-pb``
(PacBio) flags when calling SqueezeMeta. We also include
the Canu and Flye assemblers, which are specially tailored to work with
long, noisy reads. They can be selected by including the ``-a canu`` or
``-a flye`` flag when calling SqueezeMeta. As a shortcut, the ``-–minion``
flag will use both Canu and minimap2 for Oxford Nanopore MinION reads.
As an alternative to assembly, we also provide the ``sqm_longreads.pl``
script, which will predict and annotate ORFs within individual long
reads.

Working in a low-memory environment
===================================

In our experience, assembly and DIAMOND alignment against the nr
database are the most memory-hungry parts of the pipeline. By default
SqueezeMeta will set up the right parameters for DIAMOND and the Canu
assembler based on the available memory in the system. DIAMOND memory
usage can be manually controlled via the ``-b`` parameter (DIAMOND will
consume ~5\*\ *b* Gb of memory according to the documentation, but to be
safe we set ``-b`` to *free_ram/8*). Assembly memory usage is trickier, as
memory requirements increase with the number of reads in a sample. We
have managed to run SqueezeMeta with as much as 42M 2x100 Illumina HiSeq
pairs on a virtual machine with only 16Gb of memory. Conceivably, larger
samples could be split an assembled in chunks using the merged mode. We
include the shortcut flag ``-–lowmem``, which will set DIAMOND block size
to 3, and Canu memory usage to 15Gb. This is enough to make SqueezeMeta
run on 16Gb of memory, and allows the *in situ* analysis of Oxford
Nanopore MinION reads. Under such computational limitations, we have
been able to coassemble and analyze 10 MinION metagenomes (taken from
SRA project
`SRP163045 <https://www.ncbi.nlm.nih.gov/sra/?term=SRP163045>`_) in
less than 4 hours.

Tips for working in a computing cluster
=======================================

SqueezeMeta will work fine inside a computing cluster, but there are
some extra things that must be taken into account. Here is a list in
progress based on frequent issues that have been reported.

- Run ``test_install.pl`` to make sure that everything is properly configured

- If using the conda environment, make sure that it is properly activated by your batch script

- If an administrator has set up SqueezeMeta for you (and you have no write privileges in the installation directory), make sure they have run ``make_databases.pl``, ``download_databases.pl`` or ``configure_nodb.pl`` according to the installation instructions. Once again, ``test_install.pl`` should tell you whether things seem to be ok

- Make sure to request enough memory. See the previous section for a rough guide on what is “enough”. If you get a crash during the assembly or during the annotation step, it will be likely because you ran out of memory

- Make sure to manually set the ``-b`` parameter so that it matches the amount of memory that you requested divided by 8. Otherwise, SqueezeMeta will assume that it can use all the free memory in the node in which it is running. This is fine if you got a full node for yourself, but will lead to crashes otherwise

**************************************
Execution, restart and running scripts
**************************************

Scripts location
================

The scripts composing the SqueezeMeta pipeline can be found in the
``/path/to/SqueezeMeta/scripts`` directory. Other utility scripts can be
found in the ``/path/to/SqueezeMeta/utils`` directory.
See `here <https://squeezemeta.readthedocs.io/en/stable/utils.html>`_
for more information on utility scripts.

Execution
=========

The command for running SqueezeMeta has the following syntax:

``SqueezeMeta.pl -m <mode> -p <projectname> -s <equivfile> -f <raw fastq dir> <options>``

Arguments
=========

**Mandatory parameters**

[-m <sequential|coassembly|merged|seqmerge>]
    Mode: See *Choosing an assembly strategy*. (REQUIRED)

[-p <string>]
    Project name (REQUIRED in coassembly and merged modes)

[-s|samples <path>]
    Samples file (REQUIRED)

[-f|-seq <path>]
    Fastq read files’ directory (REQUIRED)

**Restarting**

[-–restart]
    Restarts the given project where it stopped (project must be speciefied with the ``-p`` option) (will NOT overwite previous results, unless ``-–force_overwrite`` is also provided)

[-step <int>]
    In combination with ``–-restart``, restarts the project starting in the given step number (combine with ``force_overwrite`` to regenerate results)

[-–force_overwrite]:
    Do not check for previous results, and overwrite existing ones

**Filtering**

[-–cleaning]
    Filters the input reads with Trimmomatic

[-cleaning_options <string>]
    Options for Trimmomatic (default: ``"LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30"``).
    Please provide all options as a single quoted string

**Assembly**

[-a <megahit|spades|rnaspades|spades-base|canu|flye>]
    assembler (default: ``megahit``)

[-assembly_options <string>]
    Extra options for the assembler (refer to the manual of the specific assembler).
    Please provide all the extra options as a single quoted string
    (e.g. ``-assembly_options "–opt1 foo –opt2 bar"``)

[-c|-contiglen <int>]
    Minimum length of contigs (default: ``200``)

[-extassembly <path>]
    Path to a file containing an external assembly provided by the user. The file must contain contigs
    in the fasta format. This overrides the assembly step of SqueezeMeta

[-extbins <path>]
    Path to a directory containing external genomes/bins provided by the user.
    There must be one file per genome/bin, each containing contigs in the fasta format.
    This overrides the assembly and binning steps

[-–sq|-–singletons]
    Unassembled reads will be treated as contigs and
    included in the contig fasta file resulting from the assembly. This
    will produce 100% mapping percentages, and will increase BY A LOT the
    number of contigs to process. Use with caution

[-contigid <string>]
    Prefix id for contigs (default: *assembler name*)

[–-norename]
    Don't rename contigs (Use at your own risk, characters like ``-`` in contig names may make the pipeline crash)

**Annotation**

[-g <int>]
    Number of targets for DIAMOND global ranking during taxonomic assignment (default: ``100``)

[-db <path>]
    Specifies the location of a new taxonomy database (in DIAMOND format, .dmnd)

[–-nocog]
    Skip COG assignment

[-–nokegg]
    Skip KEGG assignment

[-–nopfam]
    Skip Pfam assignment

[-–fastnr]
    Run DIAMOND in ``-–fast`` mode for taxonomic assignment

[-–euk]
    Drop identity filters for eukaryotic annotation (Default: no). This is recommended for analyses in which the eukaryotic
    population is relevant, as it will yield more annotations (see the
    `documentation <https://squeezemeta.readthedocs.io/en/stable/alg_details.html#taxonomic-annotation-of-eukaryotic-orfs>`_
    for details).
    Note that, regardless of whether this option is selected or not, that result will be available as part of the aggregated
    taxonomy tables generated at the last step of the pipeline and also when loading the project into
    `SQMtools <https://squeezemeta.readthedocs.io/en/stable/SQMtools.html>`_
    so this is only relevant if you are planning to use the intermediate files directly.

[-consensus <float>]
    Minimum percentage of genes assigned to a taxon in order to assign it as the consensus taxonomy
    for that contig (default: ``50``)

[-extdb <path>]
    File with a list of additional user-provided databases for functional annotations. See `Using external databases for functional annotation <https://squeezemeta.readthedocs.io/en/stable/adv_annotation.html#using-external-databases-for-functional-annotation>`_

[–D|–-doublepas]
    Run BlastX ORF prediction in addition to Prodigal (Default: no)

**Mapping**

[-map <bowtie|bwa|minimap2-ont|minimap2-pb|minimap2-sr>]
    Read mapper (default: ``bowtie``)

[-mapping_options <string>]
    Extra options for the mapper (refer to the manual of the specific mapper).
    Please provide all the extra options as a single quoted string
    (e.g. ``-mapping_options "–opt1 foo –opt2 bar"``)

**Binning**

[-binners <string>]
    Comma-separated list with the binning programs to be used (available:
    maxbin, metabat2, concoct) (default: ``concoct,metabat2``)

[–-nobins]
    Skip all binning (Default: no). Overrides ``-binners``

[-–onlybins]
    Run only assembly, binning and bin statistics
    (including GTDB-Tk if requested)

[-extbins <path>]
    Path to a directory containing external genomes/bins provided by the user.
    There must be one file per genome/bin, each containing contigs in the fasta format.
    This overrides the assembly and binning steps

[-–nomarkers]
    Skip retrieval of universal marker genes from bins.
    Note that, while this precludes recalculation of bin
    completeness/contamination in SQMtools for bin refining, you will still
    get completeness/contamination estimates of the original bins obtained
    in SqueezeMeta

[-–gtdbtk]
    Run GTDB-Tk to classify the bins. Requires
    a working GTDB-Tk installation available in your environment

[-gtdbtk_data_path <path>]
    Path to the GTDB database, by default it is assumed to be present in
    ``/path/to/SqueezeMeta/db/gtdb``. Note that the GTDB database is NOT
    included in the SqueezeMeta databases, and must be obtained separately

**Performance**

[-t <integer>]
    Number of threads (default: ``12``)

[-b|-block-size <float>]
    Block size for DIAMOND against the nr database (default: *calculate automatically*)

[-canumem <float>]
    Memory for Canu in Gb (default: ``32``)

[-–lowmem]
    Attempt to run on less than 16 Gb of RAM memory.
    Equivalent to: ``-b 3 -canumem 15``. Note that assembly may still fail due to lack of memory

**Other**

[-–minion]
    Run on MinION reads. Equivalent to
    ``-a canu -map minimap2-ont``. If canu is not working for you consider using
    ``-a flye -map minimap2-ont`` instead

[-test <integer>]
    For testing purposes, stops AFTER the given step number

[-–empty]
    Create an empty directory structure and configuration files WITHOUT
    actually running the pipeline

**Information**

[-v]
    Display version number

[-h]
    Display help

Example SqueezeMeta call
========================

``SqueezeMeta.pl -m coassembly -p test -s test.samples -f mydir --nopfam -miniden 50``

This will create a project “test” for co-assembling the samples
specified in the file “test.samples”, using a minimum identity of 50%
for taxonomic and functional assignment, and skipping Pfam annotation.
The ``-p`` parameter indicates the name under which all results and data
files will be saved. This is not required for sequential mode, where the
name will be taken from the samples file instead. The ``-f`` parameter
indicates the directory where the read files specified in the sample
file are stored.

The samples file
================

The samples file specifies the samples, the names of their corresponding
raw read files and the sequencing pair represented in those files,
separated by tabulators.

It has the format: ``<Sample>   <filename>  <pair1|pair2>``

An example would be

::

   Sample1 readfileA_1.fastq   pair1
   Sample1 readfileA_2.fastq   pair2
   Sample1 readfileB_1.fastq   pair1
   Sample1 readfileB_2.fastq   pair2
   Sample2 readfileC_1.fastq.gz    pair1
   Sample2 readfileC_2.fastq.gz    pair2
   Sample3 readfileD_1.fastq   pair1   noassembly
   Sample3 readfileD_2.fastq   pair2   noassembly

The first column indicates the sample id (this will be the project name
in sequential mode), the second contains the file names of the
sequences, and the third specifies the pair number of the reads. A
fourth optional column can take the ``noassembly`` value, indicating
that these sample must not be assembled with the rest (but will be
mapped against the assembly to get abundances). This is the case for
RNAseq reads that can hamper the assembly but we want them mapped to get
transcript abundance of the genes in the assembly. Similarly, an extra
column with the ``nobinning`` value can be included in order to avoid
using those samples for binning. Notice that a sample can have more than
one set of paired reads. The sequence files can be in fastq or fasta
format, and can be gzipped. If a sample contains paired libraries, it is
the user’s responsability to make sure that the forward and reverse
files are truly paired (i.e. they contain the same number of reads in
the same order). Some quality filtering / trimming tools may produce
unpaired filtered fastq files from paired input files (particularly if
run without the right parameters). This may result in SqueezeMeta
failing or producing incorrect results.

Restart
=======

Any interrupted SqueezeMeta run can be restarted using the program the
flag ``--restart``. It has the syntax:

``SqueezeMeta.pl -p <projectname> --restart``

This command will restart the run of that project by reading the
progress.txt file to find out the point where the run stopped.

Alternatively, the run can be restarted from a specific step by issuing
the command:

``SqueezeMeta.pl -p <projectname> --restart -step <step_to_restart_from>``

By default, already completed steps will not be repeated when
restarting, even if requested with ``-step``. In order to repeat already
completed steps you must also provide the flag ``--force_overwrite``.

e.g. ``SqueezeMeta.pl --restart -p <projectname> -step 6 --force_overwrite``
would restart the pipeline from the taxonomic assignment of genes. The
different steps of the pipeline are listed at the beginning of this documentation.
**NOTE**: When calling SqueezeMeta with ``--restart``, other parameters will be ignored.
If you want to change the configuration of your run, you will need to edit the
``/path/to/project/SqueezeMeta_conf.pl`` and change them there before calling
``SqueezeMeta.pl --restart -p <projectname>``.

Running scripts
===============

Also, any individual script of the pipeline can be run using the same
syntax:

``script <projectname>`` (for instance,
``04.rundiamond.pl <projectname>`` to repeat the DIAMOND run for the
project)

**************************
Alternative analysis modes
**************************

In addition to the main SqueezeMeta pipeline, we provide
`extra scripts <https://squeezemeta.readthedocs.io/en/stable/alt_modes.html>`_
that enable the analysis of individual reads and the annotation of sequences

1) **sqm_reads.pl**: This script performs taxonomic and functional
assignments on individual reads rather than contigs. This can be useful
when the assembly quality is low, or when looking for low abundance
functions that might not have enough coverage to be assembled.

2) **sqm_longreads.pl**: This script performs taxonomic and functional
assignments on individual reads rather than contigs, assuming that more
than one ORF can be found in the same read (e.g. as happens in PacBio or
MinION reads).

3) **sqm_hmm_reads.pl**: This script provides a wrapper to the
`Short-Pair <https://sourceforge.net/projects/short-pair/>`_ software,
which allows to screen the reads for particular functions using an
ultra-sensitive HMM algorithm.

4) **sqm_mapper.pl**: This script maps reads to a given reference using
one of the included sequence aligners (Bowtie2, BWA), and provides
estimation of the abundance of the contigs and ORFs in the reference.
Alternatively, it can be used to filter out the reads mapping to a given
reference.

5) **sqm_annot.pl**: This script performs functional and taxonomic
annotation for a set of genes, for instance these encoded in a genome
(or sets of contigs).

******************************************
Downstream analysis of SqueezeMeta results
******************************************

SqueezeMeta comes with a variety of options to explore the results and
generate different plots. These are fully described in the documentation
and in the `wiki <https://github.com/jtamames/SqueezeMeta/wiki>`_.
Briefly, the three main ways to analyze the output of SqueezeMeta are
the following:

.. image:: https://github.com/jtamames/SqueezeM/blob/images/Figure_1_readmeSQM.svg
   :width: 50%
   :align: right
   :alt: Downstream analysis of SqueezeMeta results

1) **Integration with R:** We provide the
`SQMtools <https://squeezemeta.readthedocs.io/en/stable/SQMtools.html>`_
R package, which allows to easily load a whole SqueezeMeta project and
expose the results into R. The package includes functions to select
particular taxa or functions and generate plots, as well as bindings for
other popular microbiome analysis packages such as
`microeco <https://chiliubio.github.io/microeco/>`_ and
`phyloseq <https://joey711.github.io/phyloseq/>`_. Additionally,
the package exposes all the data generated by SqueezeMeta into R so it can be
used with other third-party R packages or for custom analysis scripts.
See examples
`here <https://github.com/jtamames/SqueezeMeta/wiki/Using-R-to-analyze-your-SQM-results>`_.
**SQMtools can also be used in Mac or Windows**, meaning that you can
run SqueezeMeta in your Linux server and then move the results to your
own computer and analyze them there. See advice for this below.

2) **Integration with the anvi’o analysis pipeline:** We provide a
`compatibility layer <https://squeezemeta.readthedocs.io/en/stable/utils.html#integration-with-anvi-o>`_
for loading SqueezeMeta results into the anvi’o
analysis and visualization platform
(http://merenlab.org/software/anvio/). This includes a built-in query
language for selecting the contigs to be visualized in the anvi’o
interactive interface. See examples
`here <https://github.com/jtamames/SqueezeMeta/wiki/Loading-SQM-results-into-anvi'o>`_.

We also include `utility scripts <https://squeezemeta.readthedocs.io/en/stable/utils.html#integration-with-external-tools>`_
for generating `itol <https://itol.embl.de/>`_ and
`pavian <https://ccb.jhu.edu/software/pavian/>`_ -compatible outputs.

Analyzing SqueezeMeta results in your desktop computer
======================================================

Many users run SqueezeMeta remotely (e.g. in a computing cluster).
However it is easier to explore the results interactively from your own
computer. Since version 1.6.2, we provide an easy way to achieve this.

1) In the system in which you ran SqueezeMeta, run the utility script
``sqm2zip.py /path/to/my_project /output/dir``, where
``/path/to/my_project`` is the path to the output of SqueezeMeta, and
``/output/dir`` an arbitrary output directory

2) This will generate a
file in ``/output/dir`` named ``my_project.zip``, which contains the
essential files needed to load your project into SQMtools. Transfer this
file to your desktop computer.

3) Assuming R is present in your desktop
computer, you can install `SQMtools <https://squeezemeta.readthedocs.io/en/stable/SQMtools.html>`_ with
``if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager")}; BiocManager::install("SQMtools")``.
This will work seamlessly in Windows and Mac computers, for Linux you
may need to previously install the *libcurl* development library.

4) You can load the project directly from the zip file (no need for
decompressing) with
``import(SQMtools); SQM = loadSQM("/path/to/my_project.zip")``.

*****
About
*****

SqueezeMeta is developed by Javier Tamames and Fernando Puente-Sánchez.
Feel free to contact us for support (jtamames@cnb.csic.es,
fernando.puente.sanchez@slu.se).
