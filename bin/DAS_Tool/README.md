# DAS Tool for genome resolved metagenomics

![DAS Tool](img/logo.png)

DAS Tool is an automated method that integrates the results of a flexible number of binning algorithms to calculate an optimized, non-redundant set of bins from a single assembly.

# Reference

Christian M. K. Sieber, Alexander J. Probst, Allison Sharrar, Brian C. Thomas, Matthias Hess, Susannah G. Tringe & Jillian F. Banfield (2018). [Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy.](https://www.nature.com/articles/s41564-018-0171-1) Nature Microbiology. [https://doi.org/10.1038/s41564-018-0171-1.](https://doi.org/10.1038/s41564-018-0171-1)

# Usage

```
DAS_Tool -i methodA.scaffolds2bin,...,methodN.scaffolds2bin
         -l methodA,...,methodN -c contigs.fa -o myOutput

   -i, --bins                 Comma separated list of tab separated scaffolds to bin tables.
   -c, --contigs              Contigs in fasta format.
   -o, --outputbasename       Basename of output files.
   -l, --labels               Comma separated list of binning prediction names. (optional)
   --search_engine            Engine used for single copy gene identification [blast/diamond/usearch].
                              (default: usearch)
   --write_bin_evals          Write evaluation for each input bin set [0/1]. (default: 1)
   --create_plots             Create binning performance plots [0/1]. (default: 1)
   --write_bins               Export bins as fasta files  [0/1]. (default: 0)
   --proteins                 Predicted proteins in prodigal fasta format (>scaffoldID_geneNo).
                              Gene prediction step will be skipped if given. (optional)
   --score_threshold          Score threshold until selection algorithm will keep selecting bins [0..1].
                              (default: 0.5)
   --duplicate_penalty        Penalty for duplicate single copy genes per bin (weight b).
                              Only change if you know what you're doing. [0..3]
                              (default: 0.6)
   --megabin_penalty          Penalty for megabins (weight c). Only change if you know what you're doing. [0..3]
                              (default: 0.5)
   --db_directory             Directory of single copy gene database. (default: install_dir/db)
   --resume                   Use existing predicted single copy gene files from a previous run [0/1]. (default: 0)
   --debug                    Write debug information to log file.
   -t, --threads              Number of threads to use. (default: 1)
   -v, --version              Print version number and exit.
   -h, --help                 Show this message.

```


### Input file format
- Bins [\--bins, -i]: Tab separated files of scaffold-IDs and bin-IDs.
Scaffolds to bin file example:
```
Scaffold_1	bin.01
Scaffold_8	bin.01
Scaffold_42	bin.02
Scaffold_49	bin.03
```
- Contigs [\--contigs, -c]: Assembled contigs in fasta format:
```
>Scaffold_1
ATCATCGTCCGCATCGACGAATTCGGCGAACGAGTACCCCTGACCATCTCCGATTA...
>Scaffold_2
GATCGTCACGCAGGCTATCGGAGCCTCGACCCGCAAGCTCTGCGCCTTGGAGCAGG...
```

- Proteins (optional) [\--proteins]: Predicted proteins in prodigal fasta format. Header contains scaffold-ID and gene number:
```
>Scaffold_1_1
MPRKNKKLPRHLLVIRTSAMGDVAMLPHALRALKEAYPEVKVTVATKSLFHPFFEG...
>Scaffold_1_2
MANKIPRVPVREQDPKVRATNFEEVCYGYNVEEATLEASRCLNCKNPRCVAACPVN...
```

### Output files
- Summary of output bins including quality and completeness estimates (DASTool_summary.txt).
- Scaffolds to bin file of output bins (DASTool_scaffolds2bin.txt).
- Quality and completeness estimates of input bin sets, if ```--write_bin_evals 1```  is set ([method].eval).
- Plots showing the amount of high quality bins and score distribution of bins per method, if ```--create_plots 1``` is set (DASTool_hqBins.pdf, DASTool_scores.pdf).
- Bins in fasta format if ```--write_bins 1``` is set (DASTool_bins).



### Examples: Running DAS Tool on sample data.

**Example 1:**  Run DAS Tool on binning predictions of MetaBAT, MaxBin, CONCOCT and tetraESOMs. Output files will start with the prefix *DASToolRun1*:
```
$ ./DAS_Tool  -i sample_data/sample.human.gut_concoct_scaffolds2bin.tsv,
                 sample_data/sample.human.gut_maxbin2_scaffolds2bin.tsv,
                 sample_data/sample.human.gut_metabat_scaffolds2bin.tsv,
                 sample_data/sample.human.gut_tetraESOM_scaffolds2bin.tsv
              -l concoct,maxbin,metabat,tetraESOM
              -c sample_data/sample.human.gut_contigs.fa
              -o sample_output/DASToolRun1
```

**Example 2:** Run DAS Tool again with different parameters. Use the proteins predicted in Example 1 to skip the gene prediction step, disable writing of bin evaluations, set the number of threads to 2 and score threshold to 0.6. Output files will start with the prefix *DASToolRun2*:
```
$ ./DAS_Tool -i sample_data/sample.human.gut_concoct_scaffolds2bin.tsv,
                sample_data/sample.human.gut_maxbin2_scaffolds2bin.tsv,
                sample_data/sample.human.gut_metabat_scaffolds2bin.tsv,
                sample_data/sample.human.gut_tetraESOM_scaffolds2bin.tsv
             -l concoct,maxbin,metabat,tetraESOM
             -c sample_data/sample.human.gut_contigs.fa
             -o sample_output/DASToolRun2
             --proteins sample_output/DASToolRun1_proteins.faa
             --write_bin_evals 0
             --threads 2
             --score_threshold 0.6
```


# Dependencies

- R (>= 3.2.3): https://www.r-project.org
- R-packages: data.table (>= 1.9.6), doMC (>= 1.3.4), ggplot2 (>= 2.1.0)
- ruby (>= v2.3.1): https://www.ruby-lang.org
- Pullseq (>= 1.0.2): https://github.com/bcthomas/pullseq
- Prodigal (>= 2.6.3): https://github.com/hyattpd/Prodigal
- coreutils (only macOS/ OS X): https://www.gnu.org/software/coreutils
- One of the following search engines:
	- USEARCH (>= 8.1): http://www.drive5.com/usearch/download.html
    - DIAMOND (>= 0.9.14): https://ab.inf.uni-tuebingen.de/software/diamond
	- BLAST+ (>= 2.5.0): https://blast.ncbi.nlm.nih.gov/Blast.cgi


# Quick installation

```
# Download and extract DASTool.zip archive:
unzip DAS_Tool-1.x.x.zip
cd ./DAS_Tool-1.x.x

# Install R-packages:
R CMD INSTALL ./package/DASTool_1.x.x.tar.gz

# Unzip SCG database:
unzip ./db.zip -d db

# Run DAS Tool:
./DAS_Tool -h
```


# Installation of dependent R-packages

```
$ R
> repo='http://cran.us.r-project.org' #select a repository
> install.packages('doMC', repos=repo, dependencies = T)
> install.packages('data.table', repos=repo, dependencies = T) > install.packages('ggplot2', repos=repo, dependencies = T)
> q() #quit R-session
```

After installing all dependent R-packages, the DAS Tool R-functions can be installed in a bash terminal:
```
$ R CMD INSTALL ./package/DASTool_1.x.x.tar.gz
```
...or in an R-session:
```
$ R
> install.packages('package/DASTool_1.x.x.tar.gz')
> q() #quit R-session
```

# Preparation of input files

Not all binning tools provide results in a tab separated file of scaffold-IDs and bin-IDs. A helper script can be used to convert a set of bins in fasta format to tabular scaffolds2bin file, which can be used as input for DAS Tool: `src/Fasta_to_Scaffolds2Bin.sh -h`.

### Usage:
```
Fasta_to_Scaffolds2Bin: Converts genome bins in fasta format to scaffolds-to-bin table.

Usage: Fasta_to_Scaffolds2Bin.sh -e fasta > my_scaffolds2bin.tsv

   -e, --extension            Extension of fasta files. (default: fasta)
   -i, --input_folder         Folder with bins in fasta format. (default: ./)
   -h, --help                 Show this message.
```

### Example: Converting MaxBin fasta output into tab separated scaffolds2bin file:
```
$ ls /maxbin/output/folder
maxbin.001.fasta   maxbin.002.fasta   maxbin.003.fasta...

$ src/Fasta_to_Scaffolds2Bin.sh -i /maxbin/output/folder -e fasta > maxbin.scaffolds2bin.tsv

$ head gut_maxbin2_scaffolds2bin.tsv
NODE_10_length_127450_cov_375.783524	maxbin.001
NODE_27_length_95143_cov_427.155298	maxbin.001
NODE_51_length_78315_cov_504.322425	maxbin.001
NODE_84_length_66931_cov_376.684775	maxbin.001
NODE_87_length_65653_cov_460.202156	maxbin.001
```

Some binning tools (such as CONCOCT) provide a comma separated tabular output. To convert a comma separated file into a tab separated file a one liner can be used: `perl -pe "s/,/\t/g;" scaffolds2bin.csv > scaffolds2bin.tsv`.

### Example: Converting CONCOCT csv output into tab separated scaffolds2bin file:
```
$ head concoct_clustering_gt1000.csv
NODE_2_length_147519_cov_33.166976,42
NODE_3_length_141012_cov_38.678171,42
NODE_4_length_139685_cov_35.741896,42

$ perl -pe "s/,/\tconcoct./g;" concoct_clustering_gt1000.csv > concoct.scaffolds2bin.tsv

$ head concoct.scaffolds2bin.tsv
NODE_2_length_147519_cov_33.166976	concoct.42
NODE_3_length_141012_cov_38.678171	concoct.42
NODE_4_length_139685_cov_35.741896	concoct.42
```

# Trouble shooting and FAQs

### Dependencies not found

**Problem:** All dependencies are installed and the environmental variables are set but DAS Tool still claims that specific depencendies are missing.
**Solution:** Make sure that the dependency executable names are correct. For example USEARCH has to be executable with the command
If your USEARCH binary is called differently (e.g. `usearch9.0.2132_i86linux32`) you can either rename it or add a symbolic link called usearch:

```$ ln -s usearch9.0.2132_i86linux32 usearch```
