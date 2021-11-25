<img align="right" src="https://github.com/jtamames/SqueezeM/blob/images/logo.svg" width="20%">

# SqueezeMeta: a fully automated metagenomics pipeline, from reads to bins

- Find the SqueezeMeta paper at: https://www.frontiersin.org/articles/10.3389/fmicb.2018.03349/full 
- Find a second paper on how to analyse the output of SqueezeMeta at: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03703-2
- Check some papers using SqueezeMeta: https://github.com/jtamames/SqueezeMeta/wiki/Some-papers-using-SqueezeMeta-(non-comprehensive-list)
- Make sure to [check the wiki!](https://github.com/jtamames/SqueezeMeta/wiki)

## 1. What is SqueezeMeta?

SqueezeMeta is a full automatic pipeline for metagenomics/metatranscriptomics, covering all steps of the analysis. SqueezeMeta includes multi-metagenome support allowing the co-assembly of related metagenomes and the retrieval of individual genomes via binning procedures. Thus, SqueezeMeta features several unique characteristics:

1) Co-assembly procedure with read mapping for estimation of the abundances of genes in each metagenome
2) Co-assembly of a large number of metagenomes via merging of individual metagenomes
3) Includes binning and bin checking, for retrieving individual genomes 
4) The results are stored in a database, where they can be easily exported and shared, and can be inspected anywhere using a web interface. 
5) Internal checks for the assembly and binning steps inform about the consistency of contigs and bins, allowing to spot potential chimeras. 
6) Metatranscriptomic support via mapping of cDNA reads against reference metagenomes 

SqueezeMeta can be run in three different modes, depending of the type of multi-metagenome support. These modes are:

- Sequential mode: All samples are treated individually and analysed sequentially.

- Coassembly mode: Reads from all samples are pooled and a single assembly is performed. Then reads from individual samples are mapped to the coassembly to obtain gene abundances in each sample. Binning methods allow to obtain genome bins.

- Merged mode: if many big samples are available, co-assembly could crash because of memory requirements. This mode allows the co-assembly of an unlimited number of samples, using a procedure inspired by the one used by Benjamin Tully for analysing TARA Oceans data (https://dx.doi.org/10.17504/protocols.io.hfqb3mw). Briefly, samples are assembled individually and the resulting contigs are merged in a single co-assembly. Then the analysis proceeds as in the co-assembly mode. This is not the recommended procedure (use co-assembly if possible) since the possibility of creating chimeric contigs is higher. But it is a viable alternative when standard co-assembly is not possible.

- Seqmerge mode: This is intended to work with more samples than the merged mode. Instead of merging all individual assemblies in a single step, which can be very computationally demanding, seqmerge works sequentially. First, it assembles individually all samples, as in merged mode. But then it will merge the two most similar assemblies. Similarity is measured as Amino Acid Identity values using the wonderful CompareM software by Donovan Parks. After this first merging, it again evaluates similarity and merge, and proceeds this way until all metagenomes have been merged in one. Therefore, for n metagenomes, it will need n-1 merging steps.


SqueezeMeta uses a combination of custom scripts and external software packages for the different steps of the analysis:

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
14) Binning with MaxBin
15) Binning with MetaBAT
16) Binning integration with DAS tool
17) Taxonomic assignment of bins, and check for taxonomic disparities
18) Checking of bins with CheckM
19) Merging of previous results to obtain the bin table
20) Merging of previous results to obtain the contig table
21) Prediction of kegg and metacyc patwhays for each bin
22) Final statistics for the run

Detailed information about the different steps of the pipeline can be found in the PDF manual.


## 2. Installation

SqueezeMeta is intended to be run in a x86_64 Linux OS (tested in Ubuntu and CentOS). The easiest way to install it is by using conda.

`conda create -n SqueezeMeta -c conda-forge -c bioconda -c fpusan squeezemeta`

This will create a new conda environment named SqueezeMeta, which must then be activated.

`conda activate SqueezeMeta`

When using conda, all the scripts from the SqueezeMeta distribution will be available on `$PATH`.

Alternatively, just download the latest release from the GitHub repository and uncompress the tarball in a suitable directory. The tarball includes the SqueezeMeta scripts as well as the third-party software redistributed with SqueezeMeta (see section 6). The INSTALL files contain detailed installation instructions, including all the external libraries required to make SqueezeMeta run in a vanilla Ubuntu 14.04 or CentOS7 (DVD iso) installation.

The `test_install.pl` script can be run in order to check whether the required dependencies are available in your environment.

`/path/to/SqueezeMeta/utils/install_utils/test_install.pl`

 
## 3. Downloading or building databases

SqueezeMeta uses several databases. GenBank nr for taxonomic assignment, and eggnog, KEGG and Pfam for functional assignment. 
The script *download_databases.pl* can be run to download a pre-formatted version of all the databases required by SqueezeMeta.

`/path/to/SqueezeMeta/utils/install_utils/preparing_databases/download_databases.pl /download/path`

, where `/download/path` is the destination folder. This is the recommended option, but the files are hosted in our institutional server, which can at times be unreachable.

Alternatively, the script *make_databases.pl* can be run to download from source and format the latest version of the databases.

`/path/to/SqueezeMeta/utils/install_utils/make_databases.pl /download/path`

The databases occupy 200Gb, but we recommend having at least 350Gb free disk space during the building process.

Two directories will be generated after running either `make_databases.pl` or `download_databases.pl`.
- `/download/path/db`, which contains the actual databases.
- `/download/path/test`, which contains data for a test run of SqueezeMeta.

If `download_databases.pl` or `make_databases.pl` can't find our server, you can instead run `make_databases_alt.pl` (same syntax as `make_databases.pl`) which will try to download the data from an alternative site.

If the SqueezeMeta databases are already built in another location in the system, a different copy of SqueezeMeta can be configured to use them with

`/path/to/SqueezeMeta/utils/install_utils/configure_nodb.pl /path/to/db`

, where `/path/to/db` is the route to the `db` folder that was generated by either `make_databases.pl` or `download_databases.pl`.

After configuring the databases, the `test_install.pl` can be run in order to check that SqueezeMeta is ready to work (see previous section).


## 4. Execution, restart and running scripts

### Scripts location
The scripts composing the SqueezeMeta pipeline can be found in the `/path/to/SqueezeMeta/scripts` directory. Other utility scripts can be found in the `/path/to/SqueezeMeta/utils` directory. See the PDF manual for more information on utility scripts.

### Execution

The command for running SqueezeMeta has the following syntax:

`SqueezeMeta.pl -m <mode> -p <projectname> -s <equivfile> -f <raw fastq dir> <options>`

**Arguments** 
*Mandatory parameters* 
* *-m* <sequential, coassembly, merged>: Mode (REQUIRED) 
* *-p* \<string\>: Project name (REQUIRED in coassembly and merged modes) 
* *-s*|*-samples* \<path\>: Samples file (REQUIRED) 
* *-f*|*-seq* \<path\>: Fastq read files' directory (REQUIRED) 
 
*Filtering* 
* *--cleaning*: Filters with Trimmomatic (Default: no) 
* *-cleaning_options* [string]: Options for Trimmomatic (default: "LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30"). Please provide all options as a single quoted string.
 
*Assembly*  
* *-a* [megahit,spades,rnaspades,canu,flye]: assembler (Default:megahit) 
* *-assembly_options* [string]: Extra options for the assembler (refer to the manual of the specific assembler). Please provide all the extra options as a single quoted string (e.g. _-assembly_options “--opt1 foo --opt2 bar”_)

* *-c*|*-contiglen* [number]: Minimum length of contigs (Default:200) 
* *-extassembly* [path]: Path to an external assembly provided by the user. The file must contain contigs in the fasta format. This overrides the assembly step of SqueezeMeta.
* *--singletons*: unassembled reads will be treated as contigs and included in the contig fasta file resulting from the assembly. This will produce 100% mapping percentages, and will increase BY A LOT the number of contigs to process. Use with caution (Default: no)

 
*Annotation* 
* *--nocog*: Skip COG assignment (Default: no) 
* *--nokegg*: Skip KEGG assignment (Default: no) 
* *--nopfam*: Skip Pfam assignment (Default: no) 
* *--euk*: Drop identity filters for eukaryotic annotation (Default: no). This is recommended for analyses in which the eukaryotic population is relevant, as it will yield more annotations. See the manual for details.
* *-extdb* [path]: List of additional user-provided databases for functional annotations. More information can be found in the manual
* *--D*|*--doublepass*: Run BlastX ORF prediction in addition to Prodigal (Default: no) 
 
*Mapping* 
* *-map* [bowtie,bwa,minimap2-ont,minimap2-pb,minimap2-sr]: Read mapper (Default: bowtie) 
* *-mapping_options* [string]: Extra options for the mapper (refer to the manual of the specific mapper). Please provide all the extra options as a single quoted string (e.g. _-mapping_options “--opt1 foo --opt2 bar”_)

*Binning* 
* *--nobins*: Skip binning (Default: no) 
* *--nomaxbin*: Skip MaxBin binning (Default: no) 
* *--nometabat*: Skip MetaBat2 binning (Default: no) 
 
*Performance* 
* *-t* [number]: Number of threads (Default:12) 
* *-b*|*-block-size* [number]: Block size for DIAMOND against the nr database (Default: calculate automatically) 
* *-canumem* [number]: Memory for canu in Gb (Default: 32) 
* *--lowmem*: Run on less than 16 Gb of RAM memory (Default: no). Equivalent to: -b 3 -canumem 15 
 
*Other* 
* *--minion*: Run on MinION reads (Default: no). Equivalent to -a canu -map minimap2-ont 
* *-test* [step]: For testing purposes, stops AFTER the given step number
* *--empty*: Creates an empty directory structure and configuration files. It does not run the pipeline

 
*Information* 
* *-v*: Version number 
* *-h*: Display help 
 
 
**Example SqueezeMeta call:** `SqueezeMeta.pl -m coassembly -p test -s test.samples -f mydir --nopfam -miniden 50`

This will create a project "test" for co-assembling the samples specified in the file "test.samples", using a minimum identity of 50% for taxonomic and functional assignment, and skipping Pfam annotation. The -p parameter indicates the name under which all results and data files will be saved. This is not required for sequential mode, where the name will be taken from the samples file instead. The -f parameter indicates the directory where the read files specified in the sample file are stored.

### The samples file: 

The samples file specifies the samples, the names of their corresponding raw read files and the sequencing pair represented in those files, separated by tabulators.

It has the format: `<Sample>   <filename>  <pair1|pair2>`

An example would be

```
Sample1	readfileA_1.fastq	pair1
Sample1	readfileA_2.fastq	pair2
Sample1	readfileB_1.fastq	pair1
Sample1	readfileB_2.fastq	pair2
Sample2	readfileC_1.fastq.gz	pair1
Sample2	readfileC_2.fastq.gz	pair2
Sample3	readfileD_1.fastq	pair1	noassembly
Sample3	readfileD_2.fastq	pair2	noassembly
```

The first column indicates the sample id (this will be the project name in sequential mode), the second contains the file names of the sequences, and the third specifies the pair number of the reads. A fourth optional column can take the `noassembly` value, indicating that these sample must not be assembled with the rest (but will be mapped against the assembly to get abundances). This is the case for RNAseq reads that can hamper the assembly but we want them mapped to get transcript abundance of the genes in the assembly. Similarly, an extra column with the `nobinning` value can be included in order to avoid using those samples for binning. Notice that a sample can have more than one set of paired reads. The sequence files can be in fastq or fasta format, and can be gzipped.

### Restart

Any interrupted SqueezeMeta run can be restarted using the program restart.pl. It has the syntax:

`restart.pl <projectname>`

This command will restart the run of that project by reading the progress.txt file to find out the point where the run stopped. 
 
Alternatively, the run can be restarted from a specific step by issuing the command:

`restart.pl <projectname> -step <step_to_restart_from>`

e.g. `restart.pl <projectname> -step 6` would restart the pipeline from the taxonomic assignment of genes. The different steps of the pipeline are listed in section 1.

### Running scripts
Also, any individual script of the pipeline can be run using the same syntax: 

`script <projectname>` (for instance, `04.rundiamond.pl <projectname>` to repeat the DIAMOND run for the project)


## 5. Analizing an user-supplied assembly
An user-supplied assembly can be passed to SqueezeMeta with the flag *-extassembly <your_assembly.fasta>*. The contigs in that fasta file will be analyzed by the SqueezeMeta pipeline starting from step 2.


## 6. Using external databases for functional annotation
Version 1.0 implements the possibility of using one or several user-provided databases for functional annotation. This is invoked using the *-extdb* option. Please refer to the manual for details.


## 7. Extra sensitive detection of ORFs
Version 1.0 implements the *--D* option (*doublepass*), that attempts to provide a more sensitive ORF detection by combining the Prodigal prediction with a BlastX search on parts of the contigs where no ORFs were predicted, or where predicted ORFs did not match anything in the taxonomic and functional databases.


## 8. Testing SqueezeMeta
The *download_databases.pl* and *make_databases.pl* scripts also download two datasets for testing that the program is running correctly. Assuming either was run with the directory `/download/path` as its target the test run can be executed with

`cd </download/path/test>`  
`SqueezeMeta.pl -m coassembly -p Hadza -s test.samples -f raw`

Alternatively, `-m sequential` or `-m merged` can be used.


## 9. Working with Oxford Nanopore MinION and PacBio reads
Since version 0.3.0, SqueezeMeta is able to seamlessly work with single-end reads. In order to obtain better mappings of MinION and PacBio reads agains the assembly, we advise to use minimap2 for read counting, by including the *-map minimap2-ont* (MinION) or *-map minimap2-pb* (PacBio) flags when calling SqueezeMeta.
We also include the canu assembler, which is specially tailored to work with long, noisy reads. It can be selected by including the -a *canu* flag when calling SqueezeMeta.
As a shortcut, the *--minion* flag will use both canu and minimap2 for Oxford Nanopore MinION reads. Since version 1.3 we also include flye as an optional assembler for long reads.
As an alternative to assembly, we also provide the sqm_longreads.pl script, which will predict and annotate ORFs within individual long reads.


## 10. Working in a low-memory environment
In our experience, assembly and DIAMOND against the nr database are the most memory-hungry parts of the pipeline. By default SqueezeMeta will set up the right parameters for DIAMOND and the canu assembler based on the available memory in the system. DIAMOND memory usage can be manually controlled via the *-b* parameter (DIAMOND will consume ~5\**b* Gb of memory). Assembly memory usage is trickier, as memory requirements increase with the number of reads in a sample. We have managed to run SqueezeMeta with as much as 42M 2x100 Illumina HiSeq pairs on a virtual machine with only 16Gb of memory. Conceivably, larger samples could be split an assembled in chunks using the merged mode.
We include the shortcut flag *--lowmem*, which will set DIAMOND block size to 3, and canu memory usage to 15Gb. This is enough to make SqueezeMeta run on 16Gb of memory, and allows the *in situ* analysis of Oxford Nanopore MinION reads. Under such computational limitations, we have been able to coassemble and analyze 10 MinION metagenomes (taken from SRA project [SRP163045](https://www.ncbi.nlm.nih.gov/sra/?term=SRP163045)) in less than 4 hours.


## 11. Updating SqueezeMeta
Assuming your databases are not inside the SqueezeMeta directory, just remove it, download the new version and configure it with

`/path/to/SqueezeMeta/utils/install_utils/configure_nodb.pl /path/to/db`


## 12. Downstream analysis of SqueezeMeta results
SqueezeMeta comes with a variety of options to explore the results and generate different plots. These are fully described in the PDF manual and in the [wiki](https://github.com/jtamames/SqueezeMeta/wiki). Briefly, the three main ways to analyze the output of SqueezeMeta are the following:

<img align="right" src="https://github.com/jtamames/SqueezeM/blob/images/Figure_1_readmeSQM.png" width="50%">

**1) Integration with R:** We provide the *SQMtools* R package, which allows to easily load a whole SqueezeMeta project and expose the results into R. The package includes functions to select particular taxa or functions and generate plots. The package also makes the different tables generated by SqueezeMeta easily available for third-party R packages such as *vegan* (for multivariate analysis), *DESeq2* (for differential abundance testing) or for custom analysis pipelines. See examples [here](https://github.com/jtamames/SqueezeMeta/wiki/Using-R-to-analyze-your-SQM-results).

**2) Integration with the anvi'o analysis pipeline:** We provide a compatibility layer for loading SqueezeMeta results into the anvi'o analysis and visualization platform (http://merenlab.org/software/anvio/). This includes a built-in query language for selecting the contigs to be visualized in the anvi'o interactive interface. See examples [here](https://github.com/jtamames/SqueezeMeta/wiki/Loading-SQM-results-into-anvi'o).

We also include utility scripts for generating [itol](https://itol.embl.de/) and [pavian](https://ccb.jhu.edu/software/pavian/) -compatible outputs.


## 13. Alternative analysis modes
In addition to the main SqueezeMeta pipeline, we provide two extra modes that enable the analysis of individual reads.

**1) sqm_reads.pl**: This script performs taxonomic and functional assignments on individual reads rather than contigs. This can be useful when the assembly quality is low, or when looking for low abundance functions that might not have enough coverage to be assembled.

**2) sqm_longreads.pl**: This script performs taxonomic and functional assignments on individual reads rather than contigs, assuming that more than one ORF can be found in the same read (e.g. as happens in PacBio or MinION reads).

**3) sqm_hmm_reads.pl**: This script provides a wrapper to the [Short-Pair](https://sourceforge.net/projects/short-pair/) software, which allows to screen the reads for particular functions using an ultra-sensitive HMM algorithm.

**4) sqm_mapper.pl**: This script maps reads to a given reference using one of the included sequence aligners (Bowtie2, BWA), and provides estimation of the abundance of the contigs and ORFs in the reference.

**5) sqm_annot.pl**: This script performs functional and taxonomic annotation for a set of genes, for instance these encoded in a genome (or sets of contigs).


## 14. License and third-party software
SqueezeMeta is distributed under a GPL-3 license.
Additionally, SqueezeMeta redistributes the following third-party software:
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Megahit](https://github.com/voutcn/megahit)
* [Spades](http://cab.spbu.ru/software/spades)
* [canu](https://github.com/marbl/canu)
* [prinseq](http://prinseq.sourceforge.net)
* [kmer-db](https://github.com/refresh-bio/kmer-db)
* [cd-hit](https://github.com/weizhongli/cdhit)
* [amos](http://www.cs.jhu.edu/~genomics/AMOS)
* [mummer](https://github.com/mummer4/mummer)
* [hmmer](http://hmmer.org/)
* [barrnap](https://github.com/tseemann/barrnap)
* [aragorn](http://130.235.244.92/ARAGORN/)
* [prodigal](https://github.com/hyattpd/Prodigal)
* [DIAMOND](https://github.com/bbuchfink/diamond)
* [bwa](https://github.com/lh3/bwa)
* [minimap2](https://github.com/lh3/minimap2)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [MaxBin](https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
* [MetaBAT](https://bitbucket.org/berkeleylab/metabat)
* [CONCOCT](https://github.com/BinPro/CONCOCT)
* [DAS tool](https://github.com/cmks/DAS_Tool)
* [checkm](http://ecogenomics.github.io/CheckM)
* [comparem](https://github.com/dparks1134/CompareM)
* [MinPath](http://omics.informatics.indiana.edu/MinPath)
* [RDP classifier](https://github.com/rdpstaff/classifier)
* [pullseq](https://github.com/bcthomas/pullseq)
* [Short-Pair](https://sourceforge.net/projects/short-pair/)
* [SAMtools](http://samtools.sourceforge.net/)
* [Mothur](https://mothur.org/)
* [Flye](https://github.com/fenderglass/Flye)


## 15. About
SqueezeMeta is developed by Javier Tamames and Fernando Puente-Sánchez. Feel free to contact us for support (jtamames@cnb.csic.es, fernando.puente.sanchez@slu.se).


