***************
Utility scripts
***************

Compressing a SqueezeMeta project into a zip file
=================================================

.. _sqm2zip:
sqm2zip.py
----------
This script generates a compressed zip file with all the essential information needed to load a :doc:`SqueezeMeta <execution>` project into :doc:`SQMtools`. If the directory ``/path/to/project/results/tables`` is not present, it will also run :ref:`sqm2tables.py <sqm2tables>` to generate the required tables (see below).

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
^^^^^
.. code-block:: console

  sqm2zip.py <project_path> <output_dir> [options]

Arguments
^^^^^^^^^

Mandatory parameters (positional)
"""""""""""""""""""""""""""""""""
[project_path <path>]
    Path to the SqueezeMeta run

[output_dir] <path>
    Output directory

Options
"""""""
[--trusted-functions]
    Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables. This will be ignored if the ``/path/to/project/results/tables`` directory already exists

[--ignore-unclassified]
    Ignore reads without assigned functions for TPM calculation (KO, COG, PFAM). This will be ignored if the ``/path/to/project/results/tables`` directory already exists

[--force-overwrite]
    Write results even if the output file already exists

[--doc]
    Print the documentation

Output
^^^^^^
A zip file named ``<project_name>.zip`` in the directory specified by ``output_dir``. This file can be loaded directly into :doc:`SQMtools` using the ``loadSQM`` function.

Generating summary tables
=========================

.. _sqm2tables:
sqm2tables.py
-------------

This script generates tabular outputs from a SqueezeMeta run. It will aggregate the abundances of the ORFs assigned to the same feature (be it a given taxon or a given function) and produce tables with features in rows and samples in columns. Note that if you want to create tables coming from a :ref:`sqm_reads.pl <sqm_reads>` or :ref:`sqm_longreads.pl <sqm_longreads>` run you will need to use the :ref:`sqmreads2tables.py <sqmreads2tables>` script instead.

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

.. note::
  This script is now run automatically with default parameters at the end of a SqueezeMeta run, placing the results in the ``/path/to/project/results/tables`` directory. You may still want to run it on your own if you want to use non-default parameters.

Usage
^^^^^
.. code-block:: console

  sqm2tables.py <project_path> <output_dir> [options]

Arguments
^^^^^^^^^

Mandatory parameters (positional)
"""""""""""""""""""""""""""""""""
[project_path <path>]
    Path to the SqueezeMeta run

[output_dir <path>]
    Output directory

Options
"""""""
[--trusted-functions]
    Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables

[--ignore-unclassified]
    Ignore reads without assigned functions for TPM calculation

[--force-overwrite]
    Write results even if the output file already exists

[--doc]
    Print the documentation

.. _sqm2tables output:
Output
^^^^^^
- ``<project_name>.orfs.sequences.tsv``: ORF sequences

- ``<project_name>.orfs.sequences.tsv``: contig sequences

- ``<project_name>.orf.tax.allfilter.tsv``: taxonomy of each ORF at the different taxonomic levels. Minimum identity cutoffs for taxonomic assignment are applied to all taxa

- ``<project_name>.orf.tax.prokfilter.tsv``: taxonomy of each ORF at the different taxonomic levels. Minimum identity cutoffs for taxonomic assignment are applied to bacteria and archaea, but not to eukaryotes

- ``<project_name>.orf.tax.nofilter.tsv``: taxonomy of each ORF at the different taxonomic levels. No identity cutoffs for taxonomic assignment are applied

- ``<project_name>.orf.marker.genes.tsv``: CheckM1 marker genes present in each ORF

- ``<project_name>.orf.16S.tsv``: RDP taxonomy of the ORFs containing a 16S rRNA gene according to barrnap

- ``<project_name>.contig.tax.allfilter.tsv``: consensus taxonomy of each contig at the different taxonomic levels, based on the taxonomy of their constituent ORFs (applying minimum identity cutoffs to all taxa)

- ``<project_name>.contig.tax.prokfilter.tsv``: consensus taxonomy of each contig at the different taxonomic levels, based on the taxonomy of their constituent ORFs. Minimum identity cutoffs for taxonomic assignment are applied to bacteria and archaea, but not to eukaryotes)

- ``<project_name>.contig.tax.nofilter.tsv``: consensus taxonomy of each contig at the different taxonomic levels, based on the taxonomy of their constituent ORFs. No identity cutoffs for taxonomic assignment are applied

- ``<project_name>.bin.tax.tsv``: consensus taxonomy of each bin at the different taxonomic levels, based on the taxonomy of their constituent contigs

.. note::
   See a deeper discussion on the use of identity cutoffs in taxonomic annotation :ref:`here <euk annot>`.

- ``<project_name>.RecA.tsv``: coverage of RecA (COG0468) in the different samples.

- For each taxonomic rank (superkingdom, phylum, class, order, family, genus, species) the script will produce the following files:
    - ``<project_name>.<rank>.allfilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples, applying the identity cutoffs for taxonomic assignment
    - ``<project_name>.<rank>.prokfilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples. Identity cutoffs for taxonomic assignment are applied to prokaryotic (bacteria + archaea) ORFs but not to Eukaryotes
    - ``<project_name>.<rank>.nofilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples. Identity cutoffs for taxonomic assignment are applied to prokaryotic (bacteria + archaea) ORFs but not to Eukaryotes

- For each functional classification system (KO, COG, PFAM, and any external database provided by the user) the script will produce the following files:
    - ``<project_name>.<classification>.names.tsv``: extended description of the functional categories in that classification system. For KO and COG the file will contain three fields: ID, Name and Path within the functional hierarchy. For external databases, it will contain only ID and Name.
    - ``<project_name>.<classification>.abunds.tsv``: raw read counts of each functional category in the different samples
    - ``<project_name>.<classification>.bases.tsv``: raw base counts of each functional category in the different samples
    - ``<project_name>.<classification>.copyNumber.tsv``: average copy numbers per genome of each functional category in the different samples. Copy numbers are obtained by dividing the aggregate coverage of each function in each sample by the coverage of RecA (COG0468) in each sample.
    - ``<project_name>.<classification>.tpm.tsv``: normalized (TPM) abundances of each functional category in the different samples. This normalization takes into account both sequencing depth and gene length
.. note::
  The ``--ignore_unclassified`` flag can be used to control whether unclassified ORFs are counted towards the total for TPM normalization

.. note::
  There are more advanced ways of calculating copy numbers than normalizing by RecA coverage. These can be accessed through :doc:`SQMtools`

.. _sqmreads2tables:
sqmreads2tables.py
------------------
This script generates tabular outputs from a sqm_reads.pl or sqm_longreads.pl run. It will aggregate the abundances of the ORFs assigned to the same feature (be it a given taxon or a given function) and produce tables with features in rows and samples in columns. It can optionally accept a query argument to generate tables containing only certain taxa and functions.

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
^^^^^

.. code-block:: console

  sqmreads2tables.py <project_path> <output_dir> [options]

Arguments
^^^^^^^^^

Mandatory parameters (positional)
"""""""""""""""""""""""""""""""""
[project_path <path>]
    Path to the SqueezeMeta run

[output_dir <path>]
    Output directory

Options
"""""""

[-q/—query <string>]
    Filter the results based on the provided query in order to create tables containing only certain taxa or functions. See :ref:`query syntax`

[--trusted-functions]
    Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables

[--force-overwrite]
    Write results even if the output file already exists

[--doc]
    Print the documentation

Output
^^^^^^
- For each taxonomic rank (superkingdom, phylum, class, order, family, genus, species) the script will produce the following files:
    - ``<project_name>.<rank>.allfilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples, applying the identity cutoffs for taxonomic assignment
    - ``<project_name>.<rank>.prokfilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples. Identity cutoffs for taxonomic assignment are applied to prokaryotic (bacteria + archaea) ORFs but not to Eukaryotes. See :ref:`euk annot`
    - ``<project_name>.<rank>.nofilter.abund.tsv``: raw abundances of each taxon for that taxonomic rank in the different samples. Identity cutoffs for taxonomic assignment are applied to prokaryotic (bacteria + archaea) ORFs but not to Eukaryotes

- For each functional classification system (KO, COG, PFAM, and any external database provided by the user) the script will produce the following files:
    - ``<project_name>.<classification>.abunds.tsv``: raw abundances of each functional category in the different samples
    - ``<project_name>.<classification>.names.tsv``: extended description of the functions in that classification system. For KO and COG the file will contain three fields: ID, Name and Path within the functional hierarchy. For external databases, it will contain only ID and Name

.. _query syntax:
Query syntax
^^^^^^^^^^^^

.. note::                                                                                                                              This syntax is used by two different scripts:
  - :ref:`sqmreads2tables.py <sqmreads2tables>` script, in order to filter reads annotated with :ref:`sqm_reads_pl <sqm_reads>` or :ref:`sqm_longreads.pl <sqm_longreads>`
  - :ref:`anvi-filter-sqm.py <anvi-filter-sqm>` script, in order to filter an anvi'o database obtained after running :ref:`anvi-load-sqm.py <anvi-load-sqm>` on a SqueezeMeta project

- Please enclose query strings within double brackets.

- Queries are combinations of relational operations in the form of ``<SUBJECT> <OPERATOR> <VALUE>`` (e.g. ``"PHYLUM == Bacteroidota"``) joined by logical operators (``AND``, ``OR``).

- Parentheses can be used to group operations together.

- The ``AND`` and ``OR`` logical operators can't appear together in the same expression. Parentheses must be used to separate them into different  expressions. e.g:
    - ``"GENUS == Escherichia OR GENUS == Prevotella AND FUN CONTAINS iron"`` would not be valid. Parentheses must be used to write either of the following expressions:
        - ``"(GENUS == Escherichia OR GENUS == Prevotella)" AND FUN CONTAINS iron"`` to select features from either *Escherichia* or *Prevotella* which contain ORFs related to iron
        - ``"GENUS == Escherichia OR (GENUS == Prevotella AND FUN CONTAINS iron)"`` to select all features from *Escherichia* and any feature from *Prevotella* which contains ORFs related to iron

- Another example query would be: ``"(PHYLUM == Bacteroidota OR CLASS IN [Alphaproteobacteria, Gammaproteobacteria]) AND FUN CONTAINS iron AND Sample1 > 1"``
    - This would select all the features assigned to either the *Bacteroidota* phylum or the *Alphaproteobacteria* or *Gammaproteobacteria* classes, that also contain the substring ``"iron"`` in the functional annotations of any of their ORFs, and whose abundance in Sample1 is higher than 1

- Possible subjects are:
    - ``FUN``: search within all the combined databases used for functional annotation
    - ``FUNH``: search within the KEGG BRITE and COG functional hierarchies (e.g. ``"FUNH CONTAINS Carbohydrate metabolism"`` will select all the feature containing a gene associated with the broad ``"Carbohydrate metabolism"`` category)
    - ``SUPERKINGDOM``, ``PHYLUM``, ``CLASS``, ``ORDER``, ``FAMILY``, ``GENUS``,  ``SPECIES``: search within the taxonomic annotation at the requested taxonomic rank
    - *<SAMPLE_NAME>* (for :ref:`anvi-filter-sqm.py <anvi-filter-sqm>` only): search within the anvi'o abundances (mean coverage of a split divided by the overall sample mean coverage) in the sample named *<SAMPLE_NAME>* (e.g. if you have two samples named ``Sample1`` and ``Sample2``, the query string ``Sample1 > 0.5 AND Sample2 > 0.5`` would return the splits with an anvi'o abundance higher than 0.5 in both samples)

- Posible relational operators are ``==``,, ``!=``, ``>=``, ``<=``, ``>``, ``<``, ``IN``, ``NOT IN``, ``CONTAINS``, ``DOES NOT CONTAIN``


combine-sqm-tables.py
---------------------
Combine tabular outputs from different projects generated either with SqueezeMeta or *sqm_(long)reads* (but not both at the same time). If the directory ``/path/to/project/results/tables`` is not present, it will also run :ref:`sqm2tables.py <sqm2tables>` or :ref:`sqmreads2tables.py <sqmreads2tables>` to generate the required tables.

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

.. note::
  The recommended way of doing is is now using :doc:`SQMtools`
  - The ``loadSQM`` function accepts an arbitrary number of SqueezeMeta projects, loading them into a single SQM object
  - The ``combineSQMlite`` fucntion can be used to combine previously loaded SqueezeMeta and *sqm_(long)reads* projects into a single object. An advantage of this over ``combine-sqm-tables.py`` is that it can be used to combine projects coming from **both** SqueezeMeta and *sqm_(long)reads* at the same time. 

Usage
^^^^^
.. code-block:: console

  combine-sqm-tables.py <project_paths> [options]

Arguments
^^^^^^^^^

Mandatory parameters (positional)
"""""""""""""""""""""""""""""""""
[project_paths <paths>]
    A space-separated list of paths

Options
"""""""
[-f|--paths-file <path>]
   File containing the paths of the SqueezeMeta projects to combine, one path per line

[-o|--output-dir <path>]
    Output directory (default: ``"combined"``)

[-p|--output-prefix]
    Prefix for the output files (default: ``"combined"``)

[--trusted-functions]
   Include only ORFs with highly trusted KEGG and COG assignments in aggregated functional tables. This will be ignored if the ``/path/to/project/results/tables`` directory already exists

[--ignore-unclassified]
    Ignore reads without assigned functions for TPM calculation. This will be ignored if the ``/path/to/project/results/tables`` directory already exists or if ``--sqmreads`` is passed

[--sqmreads]
    Projects were generated using :ref:`sqm_reads.pl <sqm_reads>` or :ref:`sqm_longreads.pl <sqm_longreads>`

[--force-overwrite]
    Write results even if the output directory already exists

[--doc]
    Print the documentation

Example calls
"""""""""""""
- Combine projects  ``/path/to/proj1`` and ``/path/to/proj2`` and store output in a directory named ``"outputDir"``
    - ``combine-sqm-tables.py /path/to/proj1 /path/to/proj2 -o output_dir``
- Combine a list of projects contained in a file, use default output dir
    - ``combine-sqm-tables.py -f project_list.txt``

Output
^^^^^^
Tables containing aggregated counts and feature names for the different functional hierarchies and taxonomic levels for each sample contained in the different projects that were combined. Tables with the TPM and copy number of functions will also be generated for SqueezeMeta runs, but not for *sqm_(long)reads* runs.

Estimation of the sequencing depth needed for a project
=====================================================

.. _COVER_script:
cover.pl
--------

COVER intends to help in the experimental design of metagenomics by addressing the unavoidable question: How much should I sequence to get good results? Or the other way around: I can spend this much money, would it be worth to use it in sequencing the metagenome?

To answer these questions, COVER allows the estimation of the amount of sequencing needed to achieve a particular objective, being this the coverage attained for the most abundant N members of the microbiome. For instance, how much sequence is needed to reach 5x coverage for the four most abundant members (from now on, OTUs). COVER was first published in 2012 (Tamames *et al.*, 2012, *Environ Microbiol Rep.* **4**:335-41), but we are using a different version of the algorithm described there. Details on this implementation can be found in :ref:`COVER`.

COVER needs information on the composition of the microbiome, and that must be provided as a file containing 16S rRNA sequences obtained by amplicon sequencing of the target microbiome. If you don't have that, you can look for a similar sample already sequenced (for instance, in NCBI's SRA).

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
^^^^^
.. code-block:: console

  cover.pl -i <input_file> [options]

Arguments
^^^^^^^^^

Mandatory parameters
""""""""""""""""""""
[-i <path>]
    FASTA file containing 16S rRNA amplicons

Options
"""""""
[-idcluster <float>]
    Identity threshold for collapsing OTUs (default: ``0.98``)

[-c|-coverage <float>]
    Target coverage (default: ``5``)

[-r|-rank <integer>]
    Rank of target OTU (default: ``4``)

.. note::
   Default values imply looking for 5x coverage for the 4th most abundant 98% OTU

[-cl|-classifier <mothur|rdp>]
    Classifier to use (RDP or Mothur) (default: ``mothur``)

[-d|-dir]
    Output directory (default: ``cover``)

[-t]
    Number of threads (default: ``4``)

  (Default values imply looking for 5 x coverage for the 4 th most abundant OTU)

Output
""""""
The output is a table that first lists the amount of sequencing needed, both uncorrected and corrected by the Good’s estimator:

::
  
  Needed 4775627706 bases, uncorrected
  Correcting by unobserved: 6693800053 bases

And then lists the information and coverages for each OTU, with the following columns:

- OTU: Name of the OTU
- Size: Inferred genomic size of the OTU
- Raw abundance: Number of sequences in the OTU
- Copy number: Inferred 16S rRNA copy number
- Corrected abundance: Abundance n / Σn Abundance
- Pi : Probability of sequencing a base of this OTU
- %Genome sequenced: Percentage of the genome that will be sequenced for that OTU
- Coverage: Coverage that will be obtained for that OTU
- Taxon: Deepest taxonomic annotation for the OTU

Adding new databases to an existing project
===========================================

add_database.pl
---------------
This script adds one or several new databases to the results of an existing project. The list of databases must be provided in an external database file as specified in :ref:`Using external function database`. It must be a tab-delimited file with the following format:

::

   <Database Name>	<Path to database>	<Functional annotation file>

The databases to add must also be formatted in DIAMOND format. See :ref:`Using external function database` for details. If the external database file already exists (because you already used some external databases when running SqueezeMeta), DO NOT create a new one. Instead add the new entries to the existing database file.

The script will run Diamond searches for the new databases, and then will re-run several SqueezeMeta scripts to include the new database(s) to the existing results. The following scripts will be invoked:

- :ref:`fun3 script`
- :ref:`funcover`
- :ref:`ORF table`
- :ref:`stats`
- :ref:`sqm2tables in pipeline`

The outputs of these programs will be regenerated (but all files corresponding to other databases will remain untouched).

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
^^^^^
.. code-block:: console

  add_database.pl <project_path> <database_file>

Integration with external tools
===============================

Integration with itol
---------------------

sqm2itol.pl
^^^^^^^^^^^

This script generates the files for creating a radial plot of abundances using iTOL (https://itol.embl.de/). It can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
"""""
.. code-block:: console

  sqm2itol.pl <project_path> [options]

Arguments
"""""""""

Mandatory parameters (positional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[project_path <path>]
    Path to the SqueezeMeta run

Options
~~~~~~~
[-completion <float>]
    Select only bins with percent completion above that threshold (default: ``30``)

[-contamination <float>]
    Select only bins with percent contamination below that threshold (default: ``100``)

[-classification <metacyc|kegg>]
    Functional classification to use (default: ``metacyc``)

[-functions <path>]
    File containing the name of the functions to be considered (for functional plots). For example:
    ::

     arabinose degradation
     galactose degradation
     glucose degradation

Output
""""""
The script will generate several datafiles that you must upload to https://itol.embl.de/ to produce the figure.

Integration with ipath
----------------------

sqm2ipath.pl
^^^^^^^^^^^^
This script creates data on  the existence of enzymatic reactions that can be plotted in the interactive pathway mapper iPath (http://pathways.embl.de). It can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
"""""
.. code-block:: console

  sqm2ipath.pl <project_path> [options]

Arguments
"""""""""

Mandatory parameters (positional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[project_path <path>]
    Path to the SqueezeMeta run

Options
~~~~~~~
[-taxon <string>]
    Taxon to be plotted  (default: *plot all taxa*)

[-color <string>]
    RGB color to be used in the plot (default: ``red``)

[-c|classification <cog|kegg>]
    Functional classification to use (default: ``kegg``)

[-functions <file>]
    File containing the COG/KEGG identifiers of the functions to be considered. For example:
    ::

     K00036
     K00038
     K00040
     K00052  #ff0000
     K00053

    A second argument following the identifier selects the RGB color to be associated to that ID in the plot

    .. note::
    The plotting colors can be specified by the -color option, or by associating values to each of the IDs in the functions file. In that case, several colors can be used in the same plot. If no color is specified, the default is red.

[-o|out <path>]
    Name of the output file (default: ``ipath.out``)

Output
""""""
A file suitable to be uploaded to http://pathways.embl.de. Several output files can be combined, for instance using different colors for different taxa.

Integration with pavian
-----------------------

sqm2pavian.pl
^^^^^^^^^^^^^
This script produces output files containing abundance of taxa that can be plotted using
the Pavian tool (https://github.com/fbreitwieser/pavian). It works with projects generated with :doc:`SqueezeMeta.pl <execution>`, :ref:`sqm_reads.pl <sqm_reads>` or :ref:`sqm_longreads.pl <sqm_longreads>`. It can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.


Usage
"""""
.. code-block:: console

  sqm2pavian.pl <project_path> [mode]

Arguments
"""""""""

Mandatory parameters (positional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[project_path <path>]
    Path to the SqueezeMeta run

Options (positional)
~~~~~~~~~~~~~~~~~~~~
[mode <reads|bases>]
    Count abundances in reads or bases (default: ``reads``)

Output
""""""
A file named ``<project>.pavian`` that can be uploaded in the pavian app (https://fbreitwieser.shinyapps.io/pavian) or in the pavian R package.

Integration with anvi`o
-----------------------

.. _anvi-load-sqm:
anvi-load-sqm.py
^^^^^^^^^^^^^^^^
This script creates an anvi’o database from a SqueezeMeta project. The database can then be filtered and visually explored using the :ref:`anvi-filter-sqm.py <anvi-filter-sqm>` script. This script can be found in the ``/path/to/SqueezeMeta/utils/anvio_utils`` directory, but if using conda it will be present in your PATH. For this script to work, anvi’o must be installed and present in your PATH.

.. note::
  This has been tested using anvi’o versions 6 and 7. Support is only for released versions, the master/develop branches of anvi’o might (and will likely) not work.

Usage
"""""
.. code-block:: console

  anvi-load-sqm.py -p <project> -o <output> [options]

Arguments
"""""""""

Mandatory parameters
~~~~~~~~~~~~~~~~~~~~
[-p|-project <path>]
    Path to the SqueezeMeta run

[-o|--output <path>]
    Output directory

Options
~~~~~~~
[--num-threads <int>]
    Number of threads (default: ``12``)

[--run-HMMS]
    Run the ``anvi-run-hmms`` command from anvi’o for identifying single-copy core genes

[--run-scg-taxonomy]
    Run the ``anvi-run-scg-taxonomy`` command from anvi’o for assigning taxonomy based on single-copy core genes

[--min-contigs-length <int>]
    Minimum length of contigs (default: ``0``)

[--min-mean-coverage <float>]
    Minimum mean coverage for contigs (default: ``0``)

[--skip-SNV-profiling]
    Skip the profiling of single nucleotide variants

[--profile-SCVs]
    Perform characterization of codon frequencies in genes

[--force-overwrite]
    Force overwrite if the output directory already exists

[--doc]
    Print the documentation

Output
""""""
- ``CONTIGS.db``, ``PROFILE.db``, ``AUXILIARY-DATA.db``: anvi’o databases
- ``<project_name>_anvio_contig_taxonomy.txt``: contig taxonomy to be loaded by :ref:`anvi.filter-sqm.py <anvi-filter-sqm>`

.. _anvi-filter-sqm:
anvi-filter-sqm.py
^^^^^^^^^^^^^^^^^^
This script filters the results of a SqueezeMeta project (previously loaded into to an anvi’o database by the :ref:`anvi-load-sqm.py <anvi-load-sqm>` script) and opens an anvi’o interactive interface to examine them. Filtering criteria can be specified by using
a simple query syntax.  This script can be found in the ``/path/to/SqueezeMeta/utils/anvio_utils`` directory, but if using conda it will be present in your PATH. For this script to work, anvi’o must be installed and present in your PATH.

.. note::                                                                                                                              This has been tested using anvi’o versions 6 and 7. Support is only for released versions, the master/develop branches of anvi’o might (and will likely) not work.

Usage
"""""
.. code-block:: console

  anvi-filter-sqm.py -p <profile_db> -c <contigs_db> -t <contigs_taxonomy_file> -q <query> [options]

Arguments
"""""""""

Mandatory parameters
~~~~~~~~~~~~~~~~~~~~
[-p|--profile-db <path>]
    anvi’o profile db, as generated by :ref:`anvi-load-sqm.py <anvi-load-sqm>`

[-c|--contigs-db <path>]
    anvi’o contigs db, as generated by :ref:`anvi-load-sqm.py <anvi-load-sqm>`

[-t|--taxonomy <path>]
    Contigs taxonomy, as generated by :ref:`anvi-load-sqm.py <anvi-load-sqm>`

[-q/—query <string>]
    Filter the results based on the provided query in order to visualize only certain taxa or functions at certain abundances. See :ref:`query syntax`

Options
~~~~~~~
[-o/--output_dir <path>]
    Output directory for the filtered anvi’o databases (default: ``filteredDB``)

[-m|--max-splits <int>]
    Maximum number of splits to be loaded into anvi'o. If the provided query returns a higher number of splits, the program will stop. By default it is set to ``25,000``, larger values may make the anvi’o interface to respond slowly. Setting ``--max-splits`` to ``0`` will allow an arbitrarily large number of splits to be loaded

[--enforce-clustering]
    Make anvi’o perform an additional clustering based on abundances across samples and sequence composition

[--extra-anvio-args]
    Extra arguments for anvi-interactive, surrounded by quotes (e.g. ``--extra-anvio-args "--taxonomic-level t_phylum --title Parrot"``

[-s <yolo|safe>]
    By default, the script uses an in-house method to subset the anvi’o databases. It's ~5x quicker than using ``anvi-split`` in anvi’o5, and works well for us. However, the night is dark and full of bugs, so if you feel that your anvi’o view is missing some information, you can call the script with ``-s safe`` parameter. This will call ``anvi-split`` which should be much safer than our hacky solution (default: ``yolo``)

[--doc]
    Print the documentation

Output
""""""
The script will produce a subsetted anvi’o database, and call ``anvi-interactive`` to open a browser visualization.

Binning refinement
------------------

.. note::
    Some binning refinement functions are also available in :doc:`SQMtools`

remove_duplicate_markers.pl
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script attempts to reduce the contamination of bins by identifying duplicated
markers (conserved genes for the given taxa that are expected to be single copy but are
found to have more than one) in them. Then, it optimizes the removal of contigs
containing these duplicated markers so that only one copy of the gene is left, and no
other markers are removed.

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
"""""
.. code-block:: console

  remove_duplicate_markers.pl <project_name> [bin_name]

If no bin name is provided, the script will run the analysis for ALL bins in the project.

Output
""""""
The scripts produces a new fasta file for the bin with the name ``refined`` in the binning
directory (usually in ``<project>/results/bins/bins``). It also runs CheckM
again to redo the statistics for the bin(s). The result of that CheckM run is stored in
``<project>/temp/checkm_nodupl.txt``

find_missing_markers.pl
^^^^^^^^^^^^^^^^^^^^^^^
This script intends to improve the completeness of the bin, using the CheckM analysis to
find contigs from the same taxa of the bin that contain missing markers (those that were
not found in any contig of the bin). The user can then decide whether or not including
these contigs in the bin.

This script can be found in the ``/path/to/SqueezeMeta/utils/`` directory, but if using conda it will be present in your PATH.

Usage
"""""
.. code-block:: console

  find_missing_markers.pl <project_name> [bin_name]

If no bin name is provided, the script will run the analysis for ALL bins in the project.

The script also sets the variable $mode that affects the selection of contigs. Mode
``relaxed`` will consider contigs from all taxa not contradicting the taxonomy of the bin,
including these that belong to higher-rank taxa (for instance, if the bin is annotated as
*Escherichia* (genus), the script will consider also contigs classified as
*Enterobacteriaceae* (family), *Gammaproteobacteria* (class), or even Bacteria
(superkingdom), since these assignments are not incompatible with the one of the bin).
Mode ``strict`` will only consider contigs belonging to the same taxa of the bin (in the
example above, only these classified as genus *Escherichia*).

Output
""""""
The script produces a list of contigs containing missing markers for the bin, sorted by the
abundance of markers.

