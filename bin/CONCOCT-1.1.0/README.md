## CONCOCT 1.1.0 [![Build Status](https://travis-ci.org/BinPro/CONCOCT.png?branch=master)](https://travis-ci.org/BinPro/CONCOCT)

A program for unsupervised binning of metagenomic contigs by using nucleotide composition,
coverage data in multiple samples and linkage data from paired end reads.

## Please Cite ##
If you use CONCOCT in your publication, please cite:

Johannes Alneberg, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie Schirmer, Joshua Quick, Umer Z Ijaz, Leo Lahti, Nicholas J Loman, Anders F Andersson & Christopher Quince. 2014. Binning metagenomic contigs by coverage and composition. *Nature Methods*, doi: 10.1038/nmeth.3103

## Documentation ##
A comprehensive documentation for concoct is hosted on [readthedocs](https://concoct.readthedocs.org).

## Basic Usage ##
Cut contigs into smaller parts
```bash
cut_up_fasta.py original_contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
```

Generate table with coverage depth information per sample and subcontig.
This step assumes the directory 'mapping' contains sorted and indexed bam files where each sample has been mapped against the original contigs.
```bash
concoct_coverage_table.py contigs_10K.bed mapping/Sample*.sorted.bam > coverage_table.tsv
```

Run concoct
```bash
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
```

Merge subcontig clustering into original contig clustering
```bash
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
```

Extract bins as individual FASTA
```bash
mkdir concoct_output/fasta_bins
extract_fasta_bins.py original_contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
```

## Support ##
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg?style=flat-square)](https://gitter.im/BinPro/CONCOCT)
If you are having trouble running CONCOCT or interpretting any results, please don't hesitate to write a question in our gitter channel.

## Contribute ##

 - Issue Tracker: [github](https://github.com/binpro/CONCOCT/issues)
 - Source Code: [github](https://github.com/binpro/CONCOCT)
