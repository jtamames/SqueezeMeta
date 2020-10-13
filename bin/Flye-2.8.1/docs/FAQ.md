Flye FAQ
========

Which datasets Flye was tested on?
----------------------------------

Flye was extensively tested on various whole genome PacBio and ONT datasets.
In particular, we used Flye to assemble PacBio's P5C3, P6C4, Sequel and Sequel II;
ONT's R7-R10 basecalled with Albacore, Guppy and Flipflop (fast mode).
We typically use raw reads, however error-corrected input is also supported. 
Flye is designed for all kinds of genomes,
for various bacteria to mammalians. We are also testing Flye on 
metagenomic datasets (mock and real). You can also check the table with the Flye 
benchmarks in the [Usage file](USAGE.md).

We have NOT extensively tested Flye on targeted sequencing (not whole meta/genome).
Also, support of very short sequences, such as viruses, phages, mitochondria
is limited. On the other hand, assembly of short plasmids IS supported
through the `--plasmids` option.

Are diploid genomes supported?
------------------------------

Currently Flye does not explicitly support diploid assemblies. If heterzygosity
is low, it should not be a problem for contiguity; however similar alternative
haplotypes could be collapsed. Likewise, SNPs and short indels between
the alternative haplotypes will not be captured. If heterozygosity is high,
Flye will likely recover alternative haplotypes, but will not phase them.
Because we do not attempt to reconstruct pseudo-haplotypes (like FALCON),
this will also reduce the overall contiguity.

Are metagenomes supported?
--------------------------

Yes, use the `--meta` option. This option also should be applied if you
expect highly non-uniform read coverage in your dataset.
In this mode, some graph procedures will be less aggressive,
which can decrease the contiguity of a single genome assembly.

Are very large genomes supported?
---------------------------------

In theory - yes, but RAM usage is the limit. We have tested Flye on
many human assemblies, which typically require ~450Gb of RAM for ONT and ~140Gb for HiFi.
Memory consumption grows linearly with genome size and reads coverage.
Thus, genomes beyond ~10Gb is size might be problemmatic to assmeble.

Are PacBio CCS/HiFi reads supported?
-------------------------------

Yes, use the `--pacbio-hifi` option. 

Is Flye suitable for assembly of very short sequences, such as viruses, phages, mitochndria etc.?
-------------------------------------------------------------------------------------------------

No guarantees. Currently, the support of very short sequences, such as viruses, 
phages, mitochondria is limited. We expect this to improve in the future versions.
Assembly of short plasmids *is* supported through the `--plasmids` option.

How much resources (CPU / RAM) do I need for my genome?
-------------------------------------------------------

For a typical bacterial assembly with ~100x read coverage, 
Flye needs <10 Gb of RAM and finishes within an hour 
using ~30 threads. This will scale linearly with the increase in
read coverage. If you coverage is above 100x, consider use
`--asm-coverage 100` to use the longest 100x reads for disjointig
assembly - this should speed things up.

Mid-size eukaryotes (like C. elegans or D. melanogaster) 
with coverage around 50x might require 2-3 hours and 30-50Gb RAM
to assemble (30 threads).

Mammalian assemblies with 40x coverage need ~450Gb of RAM (for ONT)
and typically finishes within 3-4 days using 30 threads.

Various benchmarks are also given in the [Usage file](USAGE.md).
Typically, the time and memory usage usually scale linearly with
genome size and read coverage. However, highly repetitive
genomes might require more memory and be slower to assemble.
Most of the Flye stages run in parallel, so the more threads
you use, the faster it will be. We typically use 30-50 threads
on our hardware.

Note that you can also use `--asm-coverage` option to
reduce the memory usage by sampling the longest reads
for the initial disjointig assembly.

What is minimum read coverage required for Flye?
------------------------------------------------

One can typically get satisfying assembly contiguity
with 40x+ PacBio / ONT reads, if the read length is
sufficient. You might need higher coverage to improve
the consensus quality.

Depending on the technology and read length distribution, 
you might have success with 20-30x long reads. Assembly
of datasets with coverage below 10x is not recommended.

How do I select genome size if I don't know it in advance?
----------------------------------------------------------

(Note: enome size parameter is not required since the version 2.8)

The genome size estimate is used for solid k-mer selection in the
initial disjointig assembly stage. Flye is not very sensitive to
this parameter, and the estimate could be rough.
It is ok if the parameter is within 0.5x-2x of the actual genome size.
If the final assembly size is very different from the initial guess,
consider re-running the pipeline with an updated estimate for
better results.

An alternative option is to run Flye in `--meta` mode, which uses
a different approach for solid k-mer selection. This
mode is almost independent from the genome size parameter
(you still need to provide an estimate for the selection of
some other parameters). When assembly is completed,
you can re-run in the normal mode with the inferred genome size.

I have a lot of reads, but (almost) no sequence was assembled. Why is that?
---------------------------------------------------------------------------

First, make sure that your dataset type is supported (see above), 
and the parameters are set correctly. In particular, make sure you 
are correctly using either raw or corrected reads mode, 
and the genome size parameter is reasonable.

Secondly, make sure that coverage and read length is sufficient.
Flye generally expects coverage to be more than 10x, and reads
N90 over 1kb (5kb+ recommended).

Another problem with disjointig assembly could be that solid
k-mers are not properly selected because of some bias / contamination.
For example, a lot of nonsense (artificial) reads may confuse k-mer
selection. In this case, you can try the following. First, if you have 
sufficient read coverage, try to use the `--asm-coverage` option
to subsample the longest reads, which typically have better quality.
Secondly, try the `--meta` option that represents an alternative
k-mer selection strategy.

My assembly size / contiguity is not what I expected. What parameters can I tweak to fix it?
--------------------------------------------------------------------------------------------

We designed Flye to work on a wide range of datasets
using the default parameters. We thus do not expose most of the
technical parameters to the user. This also ensures the reproducibility
of Flye assemblies in different environments.

If the quality of your assembly worse than expected, first
make sure that all *required* parameters are set correctly 
(e.g. check the FAQ questions above). Make sure that input reads
have sufficient quality, coverage and length.

A notable exception is the `--min-overlap` parameter. Intuitively,
we want keep it as high as possible (e.g. 5kb) to reduce the complexity
of a repeat graph. However, if the read length is not sufficient, 
this might lead to gaps in assembly. Flye automatically
selects this parameter based on the read length distribution,
and for the most of datasets the selected value works well.
In some rare cases, this parameter needs to be adjusted manually,
for example if the read length distribution is skewd.

Currently, the minimum overlap is automatically selected with 
the 1kb-5kb range. For some datasets (such as ultra-long ONT reads)
it makes sense to manually increase minimum overlap to 10k.
We will likely include this as an automatic rule in the 
future releases.

Can I use both PacBio and ONT reads for assembly?
-------------------------------------------------

You can do this as follows: first, run the pipeline with all your reads
in the `--pacbio-raw` mode (you can specify multiple files, no need to 
merge all you reads into one). Also add `--iterations 0` to stop the pipeline
before polishing.

Once the assembly finishes, run polishing using either PacBio or ONT reads only.
Use the same assembly options, but add `--resume-from polishing`. Here is an
example of a script that should do the job (thanks to @jvhaarst):

```
flye --pacbio-raw $PBREADS $ONTREADS --iterations 0 --out-dir $OUTPUTDIR --genome-size $SIZE --threads $THREADS
flye --pacbio-raw $PBREADS --resume-from polishing --out-dir $OUTPUTDIR  --genome-size $SIZE --threads $THREADS
```

Do I still need Illumina polishing or long-read polishing is good enough?
-------------------------------------------------------------------------

It is a somewhat difficult question to answer. Flye does include
polishing step, and it producing high quality consensus on bacterial
PacBio datasets with high coverage. For example, see this recent 
[evaluation by Ryan Wick](https://github.com/rrwick/Long-read-assembler-comparison).
On the other hand, PacBio has specialized Quiver/Arrow tools that
are more advanced, since they use PacBio-specific signal information. 
It is possible, that you can get a bit of improvement after applying these tools.

For ONT data, Flye still has 0.5-1% errors in consensus due to the systematic error 
pattern of ONT reads. One might need an external polisher 
(such as nanopolish/Medaka) to get higher consensus quality.

Illumina correction can fix many of the remaining errors and improve
the assembly quality for both PacBio and ONT, for example, using Pilon or Racon.
But it should be applied with caution to prevent over-correction of repetitive regions. 
Also see [Watson and Warr paper](https://www.nature.com/articles/s41587-018-0004-z) 
for a discussion on the assembly quality.

Should I use raw or error-corrected reads?
------------------------------------------

Flye was primarily designed and tested  using raw reads, so it is always the recommended option.
Should you decide to use error-corrected reads, it might be a good idea to perform another assembly
using raw read input and compare the results. The only exception is PacBio HiFi reads, 
which should be assembled using `--pacbio-hifi`.


Do I need to preprocess my reads in any way?
--------------------------------------------

No, usually it is not necessary. Flye automatically filters out
chimeric reads or reads with bad ends. If the read coverage is very high,
you can use the built-in `--asm-coverage` option for subsampling the longest ones.

Note that in PacBio mode, Flye assumes that the input files represent PacBio subreads,
e.g. adaptors and scraps are removed and multiple passes of the same insertion
sequence are separated. This is typically handled by PacBio instruments/toolchains,
however we saw examples of problemmatic raw -> fastq conversions, 
which resulted into incorrectl subreads. In this case, 
consider using [pbclip](https://github.com/fenderglass/pbclip) to fix your Fasta/q reads.

Are cluster environments (SGE / Slurm etc.) supported?
------------------------------------------------------

Currently, cluster environments are not supported. Flye was designed to run on a single
high-memory node, and it will be difficult to make it run in a distributed environment.
Note that Flye pipeline has multiple consecutive stages, that are could be resumed and
run on different machines, if desired.

Can I run the Flye polisher on an existing assembly?
----------------------------------------------------

Yes, you can use the `--polish-target` option. Here is an example of 
polishing using PacBio reads:

```
flye --polish-target SEQ_TO_POLISH --pacbio-raw READS --iterations NUM_ITER --out-dir OUTPUTDIR --threads THREADS
```

My question is not listed, how do I get help?
---------------------------------------------

Please post your question to the [issue tracker](https://github.com/fenderglass/Flye/issues). 
In case you prefer personal communcation, you can contact Mikhail at fenderglass@gmail.com.
If you reporting a problem, please include the `flye.log` file and provide some 
details about your dataset (if possible).
