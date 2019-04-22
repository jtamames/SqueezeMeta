Version 1.0
===================
This software was last updated 01/24/2017. For any questions, comments, suggestions, etc, please send emails
to yannisun@msu.edu.

Installation
===================
1. Untar the source file called "ShortPair.tar.gz".
2. g++ compiler is required in your Unix system. To install Short-Pair, run the
Makeme file using the following command:
make
3. Dependencies:
  1) HMMER3 (http://hmmer.janelia.org/). The binary file
hmmsearch should be in the path.
  2) Python version 2.6 or later.

Run Short-Pair
===================
To run Short-Pair, use the following command:

./Short-Pair.py -m <HMMER3 HMM file> -s <seed file corresponding to HMM file> -x <fasta file1> -y <fasta file2> 
  Options:
    -h:  show this message
    -t:  probability threshold, default: 0.4;
    -o:  output file name

1. The hmm file can contain multiple hmm models and should be in HMMER3.0's hmm file format. All the hmm files of Pfam database can be downloaded from Pfam's website.
2. If you build the hmm file yourself using hmm-build in HMMER3, please make sure you have a accession number (the line that begins with ACC) as its unique identifier. Otherwise, please manually add it. The ACC should be 7-character long. 
3. The seed file contains seed sequences in stockholm format. All the seed files of Pfam database can be downloaded from Pfam's website.
4. The nucleotide sequence file should be in fasta format. 
5. Fasta file1 contains first end of reads. Fasta file name1 should end with ".1.fasta".
6. Fasta file2 contains the other end of reads. Fasta file name2 should end with ".2.fasta".
7. The format of paired-end reads should be in ".1" and ".2" notation. An example of a paired-end read will be gnl|SRA|SRR360147.1.1 and gnl|SRA|SRR360147.1.2.
 
Output
===================
The output includes putative pairs of reads that probability are above the threshold. 
For example,
gnl|SRA|SRR360147.8116738.1 PF00069 24.6 8.1e-10 178 201 1 24 -, gnl|SRA|SRR360147.8116738.2 PF00069 25.5 4e-10 104 125 4 25 +

For each line in output file, the information of each end of paired-end reads is separated by comma(,).
The first column is read name:
gnl|SRA|SRR360147.8116738.1 is read_name.1 (First end).
gnl|SRA|SRR360147.8116738.2 is read_name.2 (Second end).
The second column is the pfam domain family.
PF00069 is pfam domain family that the read align to.
The third column is the alignment score. For example, 24.6.
The fourth column is the e-value. For example, 8.1e-10.
The fifth column is the starting of alignment on pfam domain model. For example, 178.
The sixth column is the ending of alignment on pfam domain model. For example, 201.
The seventh column is the starting of alignment on read. For example, 1.
The eight column is the ending of alignment on read. For example, 24.
The nineth column is the direction of the read on the strand. We have 2 directions: - or +.


License
===================
Copyright (C) 2017 Prapaporn Techa-Angkoon, Yanni Sun, and Jikai Lei.


