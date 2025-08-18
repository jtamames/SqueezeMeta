*******************
Advanced annotation
*******************

.. _Using external taxonomy database:
Using a user-supplied database for taxonomic annotation
=======================================================

Currently, we support using a subset of NCBI's nr database (e.g. obtained after clustering it with DIAMOND)

.. _Using external function database:
Using external databases for functional annotation
==================================================

Version 1.0 of SqueezeMeta implements the possibility of using one or several external databases (user-provided) for functional annotation. This is invoked using the ``-extdb`` option. The argument must be a file (external database file) with the following format (tab-separated fields):

::

 <Database Name>    <Path to database>    <Functional annotation file>

For example, we can create the file mydb.list containing information of two databases:

::
 
 DB1	/path/to/my/database1	/path/to/annotations/database1
 DB2	/path/to/my/database2	/path/to/annotations/database2

and give it to SqueezeMeta using ``-extdb mydb.list``.

Each database must be a fasta file of amino acid sequences, in which the sequences must have a header in the format:

::

 >ID|...|Function

Where ID can be any identifier for the entry, and Function is the associated function that will be used for annotation. For example, a KEGG entry could be something like:

::

 >WP_002852319.1|K02835
 MKEFILAKNEIKTMLQIMPKEGVVLQGDLASKTSLVQAWVKFLVLGLDRVDSTPTFSTQKYE...

You can put anything you want between the first and last pipe, because these are the only fields that matter. For instance, the previous entry could also be:

::
 
 >WP_002852319.1|KEGGDB|27/02/2019|K02835
 MKEFILAKNEIKTMLQIMPKEGVVLQGDLASKTSLVQAWVKFLVLGLDRVDSTPTFSTQKYE...

Just remember not to put blank spaces, because they act as field separators in the fasta format.
This database must be formatted for DIAMOND usage. For avoiding compatibility issues between different versions of DIAMOND, it is advisable that you use the DIAMOND that is shipped with SqueezeMeta, and is placed in the bin directory of SqueezeMeta distribution. You can do the formatting with the command:

``/path/to/SqueezeMeta/bin/diamond makedb -d /path/to/ext/db/dbname.dmnd --in /path/to/my/ext/dbname.fasta``

For each database, you can OPTIONALLY provide a file with functional annotations, such as the name of the enzyme or whatever you want. Its location must be specified in the last field of the external database file. It must have only two columns separated by tabulators, the first with the function, the second with the additional information. For instance:

``K02835	peptide chain release factor 1``

The ORF table (see :ref:`ORF table`) will show both the database ID and the associated annotation for each external database you provided.


.. _Extra sensitive ORFs:
Extra-sensitive detection of ORFs
=================================

Version 1.0 implements the -D option (doublepass), that attempts to provide a more sensitive ORF detection by combining the Prodigal prediction with a BlastX search on parts of the contigs where no ORFs were predicted, or where predicted ORFs did not match anything in the taxonomic and functional databases. The procedure starts after the usual taxonomic and functional annotation. It masks the parts of the contigs in which there is a predicted ORF with some (taxonomic and functional) annotation. The remaining sequence corresponds to zones with no ORF prediction or orphan genes (no hits). The first could correspond to missed ORFs, the second to wrongly predicted ORFs. Then a DIAMOND BlastX is run only on these zones, using the same databases. The resulting hits are added as newly predicted ORFs, and the pipeline continues taking into account these new ORFs.

The pros: This procedure is able to detect missing genes or correct errors in gene prediction (for example, these derived from frameshifts). For prokaryotic metagenomes, we estimate a gain of 2-3% in the number of ORFs. This method is especially useful when you suspect that gene prediction can underperform, for instance in cases in which eukaryotes and viruses are present. Prodigal is a prokaryotic gene predictor and its behaviour for other kingdoms is uncertain. In these cases, the gain can be higher than for prokaryotes.

The con: Since it has to do an additional DIAMOND run (and using six frame-Blastx) it slows down the analysis, especially in the case of big and/or many metagenomes.

For more details, see section :ref:`doublepass`.
