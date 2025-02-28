*********************************
Adding new binners and assemblers
*********************************

Version 1.5 of SqueezeMeta allows connecting other binning tools than the shipped ones (*MaxBin*, *MetaBAT2* and *CONCOCT*). To do so, you must create a script for putting the new binning results in a new directory under project/intermediate/binners. For instance, suppose that you have created a script ``amazingbinner.py`` that runs the new binning tool *amazingbinner*, and your project will be named ``myrun``. Your script must create a directory ``myrun/intermediate/binners/amazingbinner`` and put the fasta files resulting from the binning in that directory (these fasta files MUST have extension ``.fasta`` or ``.fa``).

Next, edit the ``SqueezeMeta_conf.pl`` file in the ``scripts`` directory of the SqueezeMeta installation. You will see a line like this:

::
  
  %binscripts=('maxbin',"$installpath/lib/SqueezeMeta/bin_maxbin.pl",'metabat2',"$installpath/lib/SqueezeMeta/bin_metabat2.pl",'concoct',"$installpath/lib/SqueezeMeta/bin_concoct.pl");

This tells SqueezeMeta the available scripts for running binners. Add your new script:

::

  %binscripts=('maxbin',"$installpath/lib/SqueezeMeta/bin_maxbin.pl",'metabat2',"$installpath/lib/SqueezeMeta/bin_metabat2.pl",'concoct',"$installpath/lib/SqueezeMeta/bin_concoct.pl",’amazingbinner’,"mylocation/amazingbinner.py");

Where ``"mylocation"`` is the directory where you put your script (You may want to move it to the ``lib/SqueezeMeta`` directory in the SqueezeMeta installation, to have all binning scripts in one place).

Next, give execution permissions to the script (for instance, ``chmod a+x mylocation/amazingbinner.py``).Note that the script can be written in any language (python in this example), as long as it is executable and provides the results in the expected place.

And you are done! Now, for running SqueezeMeta with your new binner, just mention it using the ``-binners`` option, for instance:

``SqueezeMeta.pl -p myrun -s mysamples -f /path/to/sequences -m coassembly -binners maxbin,metabat,concoct,amazingbinner``

From version 1.6, SqueezeMeta also allows the connection of other assemblers than the ones shipped with the distro (MEGAHIT, SPAdes, Canu and Flye). Here I will teach you a practical example of how to do it, showing the plugging of the IDBA-UD assembler (https://github.com/loneknightpy/idba) (Peng et al, Bioinformatics 2012, 28:111420–1428; https://doi.org/10.1093/bioinformatics/bts174).

I will assume that you already installed the IDBA-UD software and put it somewhere in your system, for instance the ``/software/idba`` directory (my choice, but you can put it wherever you want).

What you need to do is to create a script to run the assembler. Your script will be called by the SqueezeMeta pipeline. I named my script ``assembly_idba.pl``, and it looks like this:

::

  #!/usr/bin/perl
  use strict;
  print "Running IDBA assembly\n";
  #-- By default, SqueezeMeta will pass the following arguments to your script:`

  my $projectdir=$ARGV[0];	# First argument: Directory of the project
  my $sample=$ARGV[1];		# Second argument: Name of the sample (in sequential mode) or project (in the rest)
  my $par1name=$ARGV[2];		# Third argument: Name of the pair1 file
  my $par2name=$ARGV[3];		# Fourth argument: Name of the pair2 file
  my $numthreads=12;		#-- In addition, we can define other parameters, for instance number of threads
  #-- IDBA wants data as an interlaced fasta file`

  #-- Fortunately, they provide a fq2fa script converting our fastq files to that format`

  #-- But, our fastq files are gzipped. Therefore first thing is to gunzip them`

  #-- We define $g1 and $g2 as variables containing the name of the gunzipped files`
  #-- Simply remove the ".gz" extension to get the gunzipped name`

  my $g1=$par1name; $g1=~s/\.gz$//; my $g2=$par2name; $g2=~s/\.gz$//;
  my $fastafile="temp.fasta";	#-- And we define $fastafile as the resulting interlaced fasta file
  #-- Now, we gunzip the files and run the fq2fa script`

  my $merge_command="gunzip $par1name; gunzip $par2name; /software/idba/bin/fq2fa --merge $g1 $g2 $fastafile"; system($merge_command);
  #-- And then we can run the IDBA assembler, just providing the input filename ($fastafile)`

  #-- We could the desired assembler options to this command line.`

  #-- For instance, we added the number of threads`

  #-- The results will be stored in a directory we named "tempidba"`

  my $assembly_command="/software/idba/bin/idba -r $fastafile --num_threads $numthreads -o tempidba"; system($assembly_command);
  #-- Finally, we have to move the resulting fasta file to the "results" directory of the SqueezeMeta project`

  #-- IDBA names the file "scaffold.fa"`

  #-- Keep in mind that the file must be named "01.project.fasta"`

  my $mv_command="mv tempidba/scaffold.fa $projectdir/results/01.$sample.fasta"; system($mv_command);
  #-- To finish, we clean up things we don't need anymore (the temporal directory and fasta files)
  my $rm_command="rm -r tempidba; rm $g1; rm $g2"; system($rm_command);
  print "All done here! Have fun!\n";

Take into account that when SqueezeMeta calls your script, it will pass four arguments to it: The project directory, the project name, and the read files (two paired-end, gunzipped fastq or fasta files). This is probably all you need to know to call the assembler.

In the script I run a formatting script ``fq2fa`` provided by IDBA-UD, to put the runs in the format it wants them. Then I run the assembler, and finally I move the resulting contig file to the results directory in the SqueezeMeta project. This is very important because the rest of the pipeline will look for the contig file there. Also, take into account that the name of the contig file must be ``01.<project>.fasta`` (where "project" is your project name).

To plug this into SqueezeMeta, the first thing to do is to move your script to the place where all other assembly scripts are, which is the ``<installpath>/lib/SqueezeMeta`` directory (where "installpath" is the installation directory of SqueezeMeta. You will see there other scripts for running assemblers, like ``assembly_megahit.pl``, ``assembly_spades.pl``, etc). Then, edit the ``SqueezeMeta_conf.pl`` file in the ``scripts`` directory of the SqueezeMeta installation. You will see a line like this:

::

 %assemblers = ("megahit","assembly_megahit.pl","spades", "assembly_spades.pl","canu","assembly_canu.pl","flye", "assembly_flye.pl");

This line is a hash (equivalent to a dict in python), telling SqueezeMeta the names of the available assemblers and the associated scripts for running them. Just add yours. Remember that the name you specify will be the one to run the assembler:

::

  %assemblers = ("megahit","assembly_megahit.pl","spades", "assembly_spades.pl","canu","assembly_canu.pl","flye", "assembly_flye.pl",”idba”,”assembly_idba.pl”);

Save it, and you are done. Now you can run a SqueezeMeta project using your new “idba” assembler:

``SqueezeMeta.pl -m coassembly -f mydir -s mysamples.samples -p idba_test -a idba``
