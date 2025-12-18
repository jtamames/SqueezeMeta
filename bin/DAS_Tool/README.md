# DAS Tool for genome resolved metagenomics

![DAS Tool](img/logo.png)

DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly.

# Reference

Christian M. K. Sieber, Alexander J. Probst, Allison Sharrar, Brian C. Thomas, Matthias Hess, Susannah G. Tringe & Jillian F. Banfield (2018). [Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy.](https://www.nature.com/articles/s41564-018-0171-1) Nature Microbiology. [https://doi.org/10.1038/s41564-018-0171-1.](https://doi.org/10.1038/s41564-018-0171-1)

# Usage

```
DAS_Tool [options] -i <contig2bin> -c <contigs_fasta> -o <outputbasename>

Options:
   -i --bins=<contig2bin>                   Comma separated list of tab separated contigs to bin tables.
   -c --contigs=<contigs>                   Contigs in fasta format.
   -o --outputbasename=<outputbasename>     Basename of output files.
   -l --labels=<labels>                     Comma separated list of binning prediction names.
   --search_engine=<search_engine>          Engine used for single copy gene identification (diamond/blastp/usearch) [default: diamond].
   -p --proteins=<proteins>                 Predicted proteins (optional) in prodigal fasta format (>contigID_geneNo).
                                            Gene prediction step will be skipped.
   --write_bin_evals                        Write evaluation of input bin sets.
   --write_bins                             Export bins as fasta files.
   --write_unbinned                         Write unbinned contigs.
   -t --threads=<threads>                   Number of threads to use [default: 1].
   --score_threshold=<score_threshold>      Score threshold until selection algorithm will keep selecting bins (0..1) [default: 0.5].
   --duplicate_penalty=<duplicate_penalty>  Penalty for duplicate single copy genes per bin (weight b).
                                            Only change if you know what you are doing (0..3) [default: 0.6].
   --megabin_penalty=<megabin_penalty>      Penalty for megabins (weight c). Only change if you know what you are doing (0..3) [default: 0.5].
   --dbDirectory=<dbDirectory>              Directory of single copy gene database [default: db].
   --resume                                 Use existing predicted single copy gene files from a previous run.
   --debug                                  Write debug information to log file.
   -v --version                             Print version number and exit.
   -h --help                                Show this.

```


### Input file format
- Bins [\--bins, -i]: Tab separated files of contig-IDs and bin-IDs.
Contigs to bin file example:
```
Contig_1	bin.01
Contig_8	bin.01
Contig_42	bin.02
Contig_49	bin.03
```
- Contigs [\--contigs, -c]: Assembled contigs in fasta format:
```
>Contig_1
ATCATCGTCCGCATCGACGAATTCGGCGAACGAGTACCCCTGACCATCTCCGATTA...
>Contig_2
GATCGTCACGCAGGCTATCGGAGCCTCGACCCGCAAGCTCTGCGCCTTGGAGCAGG...
```

- Proteins (optional) [\--proteins]: Predicted proteins in prodigal fasta format. Header contains contig-ID and gene number:
```
>Contig_1_1
MPRKNKKLPRHLLVIRTSAMGDVAMLPHALRALKEAYPEVKVTVATKSLFHPFFEG...
>Contig_1_2
MANKIPRVPVREQDPKVRATNFEEVCYGYNVEEATLEASRCLNCKNPRCVAACPVN...
```

### Output files
- Summary of output bins including quality and completeness estimates (*_DASTool_summary.tsv).
- Contigs to bin file of output bins (*_DASTool_contigs2bin.tsv).
- Quality and completeness estimates of input bin sets, if ```--write_bin_evals```  is set (*_allBins.eval).
- Bins in fasta format if ```--write_bins``` is set (DASTool_bins).



### Examples: Running DAS Tool on sample data.

**Example 1:**  Run DAS Tool on binning predictions of MetaBAT, MaxBin, CONCOCT and tetraESOMs. Output files will start with the prefix *DASToolRun1*:
```
DAS_Tool  -i sample_data/sample.human.gut_concoct_contigs2bin.tsv,\
sample_data/sample.human.gut_maxbin2_contigs2bin.tsv,\
sample_data/sample.human.gut_metabat_contigs2bin.tsv,\
sample_data/sample.human.gut_tetraESOM_contigs2bin.tsv \
-l concoct,maxbin,metabat,tetraESOM \
-c sample_data/sample.human.gut_contigs.fa \
-o sample_output/DASToolRun1
```

**Example 2:** Run DAS Tool again with different parameters. Use the proteins predicted in Example 1 to skip the gene prediction step, output evaluations of input bins, set the number of threads to 2 and score threshold to 0.6. Output files will start with the prefix *DASToolRun2*:
```
DAS_Tool -i sample_data/sample.human.gut_concoct_contigs2bin.tsv,\
sample_data/sample.human.gut_maxbin2_contigs2bin.tsv,\
sample_data/sample.human.gut_metabat_contigs2bin.tsv,\
sample_data/sample.human.gut_tetraESOM_contigs2bin.tsv \
-l concoct,maxbin,metabat,tetraESOM \
-c sample_data/sample.human.gut_contigs.fa \
-o sample_output/DASToolRun2 \
--proteins sample_output/DASToolRun1_proteins.faa \
--write_bin_evals \
--threads 2 \
--score_threshold 0.6
```


# Dependencies

- R (>= 3.2.3): https://www.r-project.org
- R-packages: data.table (>= 1.9.6), magrittr (>= 2.0.1), docopt (>= 0.7.1)
- ruby (>= v2.3.1): https://www.ruby-lang.org
- Pullseq (>= 1.0.2): https://github.com/bcthomas/pullseq
- Prodigal (>= 2.6.3): https://github.com/hyattpd/Prodigal
- coreutils (only macOS/ OS X): https://www.gnu.org/software/coreutils
- One of the following search engines:
  - DIAMOND (>= 0.9.14): https://ab.inf.uni-tuebingen.de/software/diamond
  - BLAST+ (>= 2.5.0): https://blast.ncbi.nlm.nih.gov/Blast.cgi
  - USEARCH* (>= 8.1): http://www.drive5.com/usearch/download.html

\*) The free version of USEARCH only can use up to 4Gb RAM. Therefore, the use of DIAMOND or BLAST+ is recommended for big datasets.


# Installation

```
# Download and extract DASTool.zip archive:
unzip DAS_Tool-1.x.x.zip
cd ./DAS_Tool-1.x.x

# Unzip SCG database:
unzip ./db.zip -d db

# Run DAS Tool:
./DAS_Tool -h
```


Installation of dependent R-packages:
```
$ R
> repo='http://cran.us.r-project.org' #select a repository
> install.packages('data.table', repos=repo, dependencies = T)
> install.packages('magrittr', repos=repo, dependencies = T)
> install.packages('docopt', repos=repo, dependencies = T)
> q() #quit R-session
```

# Installation using conda or homebrew
DAS Tool now can also be installed via bioconda and homebrew.

## Bioconda

 Bioconda repository: https://bioconda.github.io/recipes/das_tool/README.html. Thanks @[keuv-grvl]("https://github.com/keuv-grvl") and @[silask]("https://github.com/SilasK")!.

Add bioconda channel:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```
Install DAS Tool using conda:
```
conda install -c bioconda das_tool
```

## Homebrew

Homebrew-bio repository: https://github.com/brewsci/homebrew-bio. Thanks @[gaberoo]("https://github.com/gaberoo")!

Install DAS Tool using Homebrew:
```
brew install brewsci/bio/das_tool
```

# Docker
It is also possible to run DAS Tool using Docker. A Docker image can be built using the Dockerfile included in the repository:
```
cd ./DAS_Tool-1.x.x
docker build -t cmks/das_tool .
```
To test the build run:
```
docker run --rm -it -v $(pwd)/sample_data:/sample_data -v $(pwd)/sample_output:/sample_output cmks/das_tool \
DAS_Tool -i /sample_data/sample.human.gut_concoct_contigs2bin.tsv,\
/sample_data/sample.human.gut_maxbin2_contigs2bin.tsv,\
/sample_data/sample.human.gut_metabat_contigs2bin.tsv,\
/sample_data/sample.human.gut_tetraESOM_contigs2bin.tsv \
-l concoct,maxbin,metabat,tetraESOM \
-c /sample_data/sample.human.gut_contigs.fa \
-o /sample_output/dockerTest
```

# Preparation of input files

Not all binning tools provide results in a tab separated file of contig-IDs and bin-IDs. A helper script can be used to convert a set of bins in fasta format to tabular contigs2bin file, which can be used as input for DAS Tool: `src/Fasta_to_Contigs2Bin.sh -h`.

### Usage:
```
Fasta_to_Contigs2Bin: Converts genome bins in fasta format to contigs-to-bin table.

Usage: Fasta_to_Contigs2Bin.sh -e fasta > my_contigs2bin.tsv

   -e, --extension            Extension of fasta files. (default: fasta)
   -i, --input_folder         Folder with bins in fasta format. (default: ./)
   -h, --help                 Show this message.
```

### Example: Converting MaxBin fasta output into tab separated contigs2bin file:
```
$ ls /maxbin/output/folder
maxbin.001.fasta   maxbin.002.fasta   maxbin.003.fasta...

$ src/Fasta_to_Contigs2Bin.sh -i /maxbin/output/folder -e fasta > maxbin.contigs2bin.tsv

$ head gut_maxbin2_contigs2bin.tsv
NODE_10_length_127450_cov_375.783524	maxbin.001
NODE_27_length_95143_cov_427.155298	maxbin.001
NODE_51_length_78315_cov_504.322425	maxbin.001
NODE_84_length_66931_cov_376.684775	maxbin.001
NODE_87_length_65653_cov_460.202156	maxbin.001
```

Some binning tools (such as CONCOCT) provide a comma separated tabular output. To convert a comma separated file into a tab separated file a one liner can be used: `perl -pe "s/,/\t/g;" contigs2bin.csv > contigs2bin.tsv`.

### Example: Converting CONCOCT csv output into tab separated contigs2bin file:
```
$ head concoct_clustering_gt1000.csv
NODE_2_length_147519_cov_33.166976,42
NODE_3_length_141012_cov_38.678171,42
NODE_4_length_139685_cov_35.741896,42

$ perl -pe "s/,/\tconcoct./g;" concoct_clustering_gt1000.csv > concoct.contigs2bin.tsv

$ head concoct.contigs2bin.tsv
NODE_2_length_147519_cov_33.166976	concoct.42
NODE_3_length_141012_cov_38.678171	concoct.42
NODE_4_length_139685_cov_35.741896	concoct.42
```

# Troubleshooting/ known issues

### Docopt issue

**Problem:** When executing DAS Tool a truncated version of the help message is displayed (see below). This is a known bug of the current version of the `docopt` R package, which occurs if the command-line syntax is violated.
```
Error: DAS Tool

Usage:
  DAS_Tool [options] -i <contig2bin> -c <contigs_fasta> -o <outputbasename>
  DAS_Tool [--help]

Options:
   -i --bins=<contig2bin>                   Comma separated list of tab separated contigs to bin tables.
   -c --contigs=<contigs>                   Contigs in fasta format.
   -o --outputbasename=<outputbasename>     Basename of output files.
   -l --labels=<labels>                     Comma separated list of binning prediction names.
   --search_engine=<search_engine>          Engine used for single copy gene identification (di
Execution halted
```

**Solution:** Check command line for any typos.



### Dependencies not found

**Problem:** All dependencies are installed and the environmental variables are set but DAS Tool still claims that specific depencendies are missing.

**Solution:** Make sure that the dependency executable names are correct. For example USEARCH has to be executable with the command
If your USEARCH binary is called differently (e.g. `usearch9.0.2132_i86linux32`) you can either rename it or add a symbolic link called usearch:

```
$ ln -s usearch9.0.2132_i86linux32 usearch
```

### Memory limit of 32-bit usearch version exceeded

**Problem:** Running DAS Tool with the free version of USEARCH on a large metagenomic dataset results in the following error:
```
---Fatal error---
Memory limit of 32-bit process exceeded, 64-bit build required
makeblastdb did not work for my_proteins.faa, please check your input file
```

**Solution:** Use DIAMOND or BLAST as alignment tool (`--search_engine diamond` or `--search_engine blast`):
