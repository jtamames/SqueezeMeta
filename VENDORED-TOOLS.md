This is a list of all the tools redistributed with SqueezeMeta, and a brief description of the custom modifications (if any) that were applied to each tool.

We vendor third-party software since
- The pipeline is complex and we want to minimize the burden on our users. We aim for SqueezeMeta to depend only on libraries that can be installed via standard packaging tools (apt, yum, etc)
- Some tools require modifications (e.g. parametrized rather than hardcoded database locations) to work well within our pipeline.

A given tool _should_ be replaceable by its original version if
- It has no custom patch listed
- It has ONLY the "Work within the SQM directory structure" patch listed

In order to control which software is called by SqueezeMeta, modify the "External software" section of the `SqueezeMeta/scripts/SqueezeMeta_conf.pl`

E.g. changing `$spades_soft = "$installpath/bin/SPAdes/spades.py";` to `$spades_soft = "spades.py";` will make SqueezeMeta use the SPAdes version in `$PATH` rather than the vendored one.

Note that some of these tools require additional software and libraries to be available via `PATH` and `LD_LIBRARY_PATH`. This is also indicated in the `SqueezeMeta_conf.pl` file.

SqueezeMeta redistributes the following third-party software:
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* [Megahit](https://github.com/voutcn/megahit)
* [Spades](http://cab.spbu.ru/software/spades)
  - Work within the SQM directory structure
* [canu](https://github.com/marbl/canu)
  - Work within the SQM directory structure
* [prinseq](http://prinseq.sourceforge.net)
* [kmer-db](https://github.com/refresh-bio/kmer-db)
* [CD-HIT](https://github.com/weizhongli/cdhit)
  - Recompile with MAX_SEQ=20000000
* [amos](http://www.cs.jhu.edu/~genomics/AMOS)
  - Work within the SQM directory structure
  - Add multithreading in nucmer calls (minimus2)
  - Add a custom minimus2 script for the SQM-seqmerge mode
* [mummer](https://github.com/mummer4/mummer)
* [hmmer](http://hmmer.org/)
* [barrnap](https://github.com/tseemann/barrnap)
  - Work within the SQM directory structure
  - Add `-dbdir` as an additional command line argument
* [aragorn](http://130.235.244.92/ARAGORN/)
* [prodigal](https://github.com/hyattpd/Prodigal)
* [DIAMOND](https://github.com/bbuchfink/diamond)
* [bwa](https://github.com/lh3/bwa)
* [minimap2](https://github.com/lh3/minimap2)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [MaxBin](https://downloads.jbei.org/data/microbial_communities/MaxBin/MaxBin.html)
  - Work within the SQM directory structure
  - Add `-markerpath` as an additional command line argument
* [MetaBAT](https://bitbucket.org/berkeleylab/metabat)
* [DAS tool](https://github.com/cmks/DAS_Tool)
  - Add extra logging, remove some superfluous error messages
  - Explicitly load `library(methods)` in DAS_Tool.R since Rscript does not load it on startup (even if R console does)
* [checkm](http://ecogenomics.github.io/CheckM)
  - Work within the SQM directory structure
  - Port to python3
* [comparem](https://github.com/dparks1134/CompareM)
  - Work within the SQM directory structure
  - Port to python3
* [MinPath](http://omics.informatics.indiana.edu/MinPath)
  - Work within the SQM directory structure
  - Port to python3
* [RDP classifier](https://github.com/rdpstaff/classifier)
* [pullseq](https://github.com/bcthomas/pullseq)
* [Short-Pair](https://sourceforge.net/projects/short-pair/)
  - Work within the SQM directory structure
  - Port to python3
* [SAMtools](http://samtools.sourceforge.net/)
* [Mothur](https://mothur.org/)
* [Flye](https://github.com/fenderglass/Flye)
