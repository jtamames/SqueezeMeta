# SqueezeM: a fully automated metagenomics pipeline, from reads to bins

## 1. What is squeezeM?

SqueezeM is a full automatic pipeline for metagenomics/metatranscriptomics, covering all steps of the analysis. SqueezeM includes multi-metagenome support allowing the co-assembly of related metagenomes and the retrieval of individual genomes via binning procedures. Thus, squeezeM features several unique characteristics:

1) Co-assembly procedure with read mapping for estimation of the abundances of genes in each metagenome
2) Co-assembly of unlimited number of metagenomes via merging of individual metagenomes
3) Includes binning and bin checking, for retrieving individual genomes 
4) the results are stored in a database, where they can be easily exported and shared, and can be inspected anywhere using a web interface. 
5) Internal checks for the assembly and binning steps inform about the consistency of contigs and bins, allowing to spot potential chimeras. 
6) Metatranscriptomic support via mapping of cDNA reads against reference metagenomes 

SqueezeM can be run in three different modes, depending of the type of multi-metagenome support. These modes are:

-Sequential mode: All samples are treated individually and analysed sequentially. This mode does not include binning.

-Coassembly mode: Reads from all samples are pooled and a single assembly is performed. Then reads from individual samples are mapped to the coassembly to obtain gene abundances in each sample. Binning methods allow to obtain genome bins.

-Merged mode: if many big samples are available, co-assembly could crash because of memory requirements. This mode allows the co-assembly of an unlimited number of samples, using a procedure inspired by the one used by Benjamin Tully for analysing TARA Oceans data (https://dx.doi.org/10.17504/protocols.io.hfqb3mw ). Briefly, samples are assembled individually and the resulting contigs are merged in a single co-assembly. Then the analysis proceeds as in the co-assembly mode. This is not the recommended procedure (use co-assembly if possible) since the possibility of creating chimeric contigs is higher. But it is a viable alternative when standard co-assembly is not possible.

SqueezeM uses a combination of custom scripts and external software packages for the different steps of the analysis:

1) Assembly 
2) RNA prediction and classification
3) ORF (CDS) prediction
4) Homology searching against taxonomic and functional databases
5) Hmmer searching against Pfam database
6) Taxonomic assignment of genes 
7) Functional assignment of genes
8) Taxonomic assignment of contigs, and chimera checking
9) Coverage and abundance estimation for genes and contigs 
10) Estimation of taxa abundances
11) Estimation of function abundances
12) Merging of previous results to obtain the ORF table
13) Binning with Maxbin
14) Binning with metabat2
15) Taxonomic assignment of bins, and chimera checking
16) Checking of bins
17) Merging of previous results to obtain the bin table
18) Merging of previous results to obtain the contig table
19) Final statistics for the run
20) Prediction of kegg and metacyc patwhays in each bin


## 2. Installation

For installing squeezeM, download the latest release from the GitHub repository and uncompress the tarball in a suitable directory. The tarball includes the squeezeM scripts as well as the third-party software redistributed with squeezeM (see section 6). The INSTALL file contains detailed installation instructions, including all the external libraries required to make squeezeM run in a vanilla Ubuntu 14.04 installation.
 
 
## 3. Building databases

SqueezeM uses several databases. GenBank nr for taxonomic assignment, and eggnog, KEGG and Pfam for functional assignment. The script make_databases.pl must be run to download and format all these databases.

`<installpath>/squeezeM/scripts/preparing_databases/make_databases.pl <datapath>`

, where `<datapath>` is the destination folder. The process will take about a day. The databases occupy 130Gb, but we recommend having at least 300Gb free disk space during the building process.


## 4. Execution, restart and running scripts
### Scripts location
The scripts composing the squeezeM pipeline can be found in the `.../squeezeM/scripts` directory. We recommend adding it to your $PATH environment variable.
### Execution

The command for running squeezeM has the following syntax:

`squeezeM.pl -m <mode> -p <projectname> -s <equivfile> -f <raw fastq dir> <options>`

**Arguments**
* -m: Mode (sequential, coassembly, merged) (REQUIRED) 
* -p: Project name (REQUIRED in coassembly and merged modes) 
* -s|-samples: Samples file (REQUIRED) 
* -f|-seq: Fastq read files' directory (REQUIRED) 
* -t: Number of threads (Default:12) 
* -a: assembler [megahit,spades] (Default:megahit) 
* -c|-contiglen: Minimum length of contigs (Default:1200)
* -map: Read mapper [bowtie,bwa,minimap2-ont,minimap2-pb,minimap2-sr] (Default: bowtie)
* --nocog: Skip COG assignment (Default: no) 
* --nokegg: Skip KEGG assignment (Default: no) 
* --nopfam: Skip Pfam assignment (Default: no) 
* --nobins: Skip binning (Default: no) 
* -e|-evalue: max evalue for diamond run (Default: 1e-03) 
* -miniden: minimum identity perc for diamond run (Default: 50) 
* --megahit_options: Options for megahit assembler 
* --spades_options: Options for spades assembler 

**Example squeezeM call:** `squeezeM.pl -m coassembly -p test -s test.samples -f mydir --nopfam -miniden 60`

This will create a project "test" for co-assembling the samples specified in the file "test.samples", using a minimum identity of 60% for taxonomic and functional assignment, and skipping Pfam annotation. The -p parameter indicates the name under which all results and data files will be saved. This is not required for sequential mode, where the name will be taken from the samples file instead. The -f parameter indicates the directory where the read files specified in the sample file are stored.

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

The first column indicates the sample id (this will be the project name in sequential mode), the second contains the file names of the sequences, and the third specifies the pair number of the reads. A fourth optional column can take the "noassembly" value, indicating that these sample must not be assembled with the rest (but will be mapped against the assembly to get abundances). This is the case for RNAseq reads that can hamper the assembly but we want them mapped to get transcript abundance of the genes in the assembly. Notice also that paired reads are expected, and that a sample can have more than one set of paired reads. The sequence files can be in fastq or fasta format, and can be gzipped.

### Restart

Any interrupted squeezeM run can be restarted using the program restart.pl. It has the syntax:

`restart.pl <projectname>`

This command must be issued in the upper directory to the project <projectname>, and will restart the run of that project by reading the progress.txt file to find out the point where the run stopped.

### Running scripts
Also, any individual script of the pipeline can be run in the upper directory to the project using the same syntax: 

`script <projectname>` (for instance, `04.rundiamond.pl <projectname>` to repeat the Diamond run for the project)

## 5. Testing squeezeM
The make_databases.pl script also downloads two datasets for testing that the program is running correctly. Assuming make_databases.pl was run with the directory `<datapath>` as its target the test run can be executed with

`cd <datapath>`
`squeezeM.pl -m coassembly -p Hadza -s test.samples -f raw`

Alternative `-m sequential`, `-m merged` can be used.

## 6. Working with Oxford MinION and PacBio reads.
Since version 0.3.0, squeezeM is able to seamlessly work with single-end reads. In order to obtain better mappings of MinION and PacBio reads agains the assembly, we advise to use minimap2 for read counting, by including the -map *minimap2-ont* (MinION) or *-map minimap2-pb* (PacBio) flags when calling squeezeM.

## 7. Setting up the MySQL database.
SqueezeM includes a built in MySQL database that can be queried via a web-based interface, in order to facilitate the exploration of metagenomic results. Code and instruction installations can be found at https://github.com/jtamames/SqueezeMdb.

## 8. License and third-party software
SqueezeM is distributed with a GPL-3 license.
Additionally, squeezeM redistributes the following third-party software:
* [Megahit](https://github.com/voutcn/megahit)
* [Spades](http://cab.spbu.ru/software/spades/)
* [prinseq](http://prinseq.sourceforge.net/)
* [prodigal](https://github.com/hyattpd/Prodigal)
* [cd-hit](https://github.com/weizhongli/cdhit)
* [amos](http://www.cs.jhu.edu/~genomics/AMOS/)
* [mummer](https://github.com/mummer4/mummer/)
* [hmmer](http://hmmer.org/)
* [diamond](https://github.com/bbuchfink/diamond)
* [bedtools](https://github.com/arq5x/bedtools2)
* [bwa](https://github.com/lh3/bwa)
* [minimap2](https://github.com/lh3/minimap2)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [barrnap](https://github.com/tseemann/barrnap)
* [maxbin](https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
* [metabat](https://bitbucket.org/berkeleylab/metabat)
* [checkm](http://ecogenomics.github.io/CheckM/)
* [MinPath](http://omics.informatics.indiana.edu/MinPath)
* [RDP classifier](https://github.com/rdpstaff/classifier)

## 9. About
SqueezeM is developed by Javier Tamames with collaboration from Fernando Puente-Sánchez. Feel free to contact us for support (jtamames@cnb.csic.es, fpuente@cnb.csic.es).


