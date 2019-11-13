#!/usr/bin/env perl -w

### Modified from MaxBin version 2.2.6 
### by Fernando Puente-Sanchez on 2019-MAR-28, for compatibility with the SqueezeMeta pipeline.
### SqueezeMeta is released under the GPL-3 license.

use strict;
use Cwd;
use LWP::Simple;
use FindBin qw($Bin);

my $currdir = getcwd;

my $VERSION_NUM = "2.2.6";

require("$Bin\/_getmarker.pl");
require("$Bin\/_getabund.pl");
require("$Bin\/_sepReads.pl");

my $MAX_MARKER = 0;
my $LOGOUT;

#Path info modified by Fernando Puente-Sánchez, 27-IV-2018
#my $SETTING_FILE = "setting";
my $BOWTIE2BUILD = "$Bin\/..\/bowtie2\/bowtie2-build";
my $BOWTIE2 = "$Bin\/..\/bowtie2\/bowtie2";
my $HMMSEARCH = "$Bin\/..\/hmmer\/hmmsearch";
my $RUNFRAG = "$Bin\/auxiliary\/run_FragGeneScan.pl";
#my $VELVETH = "velveth";
#my $VELVETG = "velvetg";
my $IDBA_UD = "$Bin\/auxiliary\/idba_ud";

checkProgram();

my $RSCRIPT = "Rscript";
my $HEATMAP_R = "$Bin\/heatmap.r";
my $MARKERHMM = "$Bin\/..\/..\/db\/marker.hmm"; #This is overriden by the -markerpath parameter.";
my $MARKERNUM = 107;
my $MAXBIN = "$Bin\/MaxBin"; # FPS
my $ABUND_OUTPUT = "";
my $THREADNUM = 1;
my $REASSEMBLY = 0;
my $RPLOT = 0;
my $KMERLEN = 55;
my $RECURSION_MAX = 5;
my $MIN_SEQ_LENGTH = 1000;
my $MIN_BIN_SIZE = 100000;
my $FINISH = "FINISH";

my $USAGE = qq(MaxBin - a metagenomics binning software.
Usage:
  run_MaxBin.pl
    -contig (contig file)
    -out (output file)

   (Input reads and abundance information)
    [-reads (reads file) -reads2 (readsfile) -reads3 (readsfile) -reads4 ... ]
    [-abund (abundance file) -abund2 (abundfile) -abund3 (abundfile) -abund4 ... ]

   (You can also input lists consisting of reads and abundance files)
    [-reads_list (list of reads files)]
    [-abund_list (list of abundance files)]

   (Other parameters)
    [-min_contig_length (minimum contig length. Default 1000)]
    [-max_iteration (maximum Expectation-Maximization algorithm iteration number. Default 50)]
    [-thread (thread num; default 1)]
    [-prob_threshold (probability threshold for EM final classification. Default 0.9)]
    [-plotmarker]
    [-markerset (marker gene sets, 107 (default) or 40.  See README for more information.)]

  (for debug purpose)
    [-version] [-v] (print version number)
    [-verbose]
    [-preserve_intermediate]

  Please specify either -reads or -abund information.
  You can input multiple reads and/or abundance files at the same time.
  Please read README file for more details.\n);


main();

sub main
{
	my $contig_f = "";
	my $contig_name = "";
	my @reads_f;
	my $readscount = 0;
	my @abund_f;
	my $abundcount = 0;
	my $reads_list = "";
	my $abund_list = "";
	my $out_f = "";
	my $outname = "";
	my $outdir = "";
	my $verbose = 0;
	my $preserve = 0;
	my $maxem = -1;
	my $remodel_contig = 0;
	my $prob_threshold = -1;
	my $old_contig_name;
	my $new_contig;
	my $isfastq;
	my @finisharr;

	my $starttime = time();
	my $endtime;
	
	my $ARGC = @ARGV;
	my $i;
	my $j;
	my $k;
	my @arr;
	my @tmparr;
	my $cmd;
	my $line;
	my $tmp;
	my $param;
	my $param2;

	# Create a temporary log file
	$k = 0;
	while ($k == 0)
	{
		$k = int(rand(100000000));
		if (-e "$k.log")
		{
			$k = 0;
		}
	}
	open(TMPLOG, ">$k.log");


	# Check if MaxBin program exist. If not then build everything, including 3rd-party software
	if (!(-e $MAXBIN))
	{
		print "Program not built. Please run \"make\" under src directory to build MaxBin program.\n";
		close(TMPLOG);
		unlink("$k.log");
		exit(-1);
	}
	
	print "MaxBin $VERSION_NUM\n";
	print TMPLOG "MaxBin $VERSION_NUM\n";
	for ($i = 0; $i < $ARGC; $i++)
	{
		if ($ARGV[$i] eq "-contig")
		{
			$i++;
			$contig_f = $ARGV[$i];
			$contig_name = $ARGV[$i];
			print "Input contig: $contig_f\n";
			print TMPLOG "Input contig: $contig_f\n";
		}
		elsif ($ARGV[$i] eq "-reads_list")
		{
			$i++;
			$reads_list = $ARGV[$i];
		}
		elsif ($ARGV[$i] eq "-abund_list")
		{
			$i++;
			$abund_list = $ARGV[$i];
		}
		elsif ($ARGV[$i] =~ /^\-reads/)
		{
			$i++;
			if (-e $ARGV[$i])
			{
				$reads_f[$readscount] = $ARGV[$i];
				$readscount++;
				print "Located reads file [$ARGV[$i]]\n";
				print TMPLOG "Located reads file [$ARGV[$i]]\n";
			}
			else
			{
				print "Cannot find reads file [$ARGV[$i]]\n";
				close(TMPLOG);
				unlink("$k.log");
				exit;
			}
		}
		elsif ($ARGV[$i] =~ /^\-abund/)
		{
			$i++;
			if (-e $ARGV[$i])
			{
				$abund_f[$abundcount] = $ARGV[$i];
				$abundcount++;
				print "Located abundance file [$ARGV[$i]]\n";
				print TMPLOG "Located abundance file [$ARGV[$i]]\n";
			}
			else
			{
				print "Cannot find abundance file [$ARGV[$i]]\n";
				close(TMPLOG);
				unlink("$k.log");
				exit;
			}
		}
		elsif ($ARGV[$i] eq "-out")
		{
			$i++;
			$out_f = $ARGV[$i];
			$j = rindex($out_f, "/");
			if ($j == -1)
			{
				$outname = $out_f;
			}
			else
			{
				$outname = substr($out_f, $j + 1);
				$outdir = substr($out_f, 0, $j);
			}
			print "out header: $ARGV[$i]\n";
			print TMPLOG "out header: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-thread")
		{
			$i++;
			$THREADNUM = $ARGV[$i];
			print "Thread: $ARGV[$i]\n";
			print TMPLOG "Thread: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-max_iteration")
		{
			$i++;
			$maxem = $ARGV[$i];
			print "Max iteration: $ARGV[$i]\n";
			print TMPLOG "Max iteration: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-prob_threshold")
		{
			$i++;
			if ($ARGV[$i] >= 0)
			{
				$prob_threshold = $ARGV[$i];
				print "Probability threshold: $ARGV[$i]\n";
				print TMPLOG "Probability threshold: $ARGV[$i]\n";
			}
			else
			{
				print "Invalid value $ARGV[$i]. Probability threshold not changed.\n";
			}
		}
		elsif ($ARGV[$i] eq "-min_contig_length")
		{
			$i++;
			$MIN_SEQ_LENGTH = $ARGV[$i];
			print "Min contig length: $ARGV[$i]\n";
			print TMPLOG "Min contig length: $ARGV[$i]\n";
		}
		elsif ($ARGV[$i] eq "-reassembly") # Reassembly currently still does not support FASTQ reads
		{
			$REASSEMBLY = 1;
			print "Reassembly: 1\n";
			print TMPLOG "Reassembly: 1\n";
		}
		elsif ($ARGV[$i] eq "-markerset")
		{
			$i++;
			if ($ARGV[$i] == 40)
			{
				print "Switch to 40 marker genes universal for bacteria and archaea.\n";
				print TMPLOG "Switch to 40 marker genes universal for bacteria and archaea.\n";
				$MARKERHMM = "$Bin/bacar_marker.hmm";
				$MARKERNUM = 40;
			}
		}
		elsif ($ARGV[$i] eq "-verbose")
		{
			$verbose = 1;
		}
		elsif ($ARGV[$i] eq "-plotmarker")
		{
			$RPLOT = 1;
		}
		elsif ($ARGV[$i] eq "-preserve_intermediate")
		{
			$preserve = 1;
		}
                elsif ($ARGV[$i] eq "-markerpath") #FPS, PARAMETRIZE DATABASE PATH!
                {
                        $i++;
                        $MARKERHMM = $ARGV[$i];
                }
		elsif ($ARGV[$i] eq "-version" || $ARGV[$i] eq "-v")
		{
			print "MaxBin $VERSION_NUM\n";
			exit;
		}
		else
		{
			print "Unrecognized token \[$ARGV[$i]\]\n";
			print $USAGE;
			close(TMPLOG);
			unlink("$k.log");
			exit;
		}
	}
	
	if ($contig_f eq "")
	{
		print "No Contig file. Please specify contig file by -contig\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if ($out_f eq "")
	{
		print "Please specify output file by -out.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}
	if (-e "$out_f.log")
	{
		unlink "$out_f.log";
	}

	# Processing reads_list and abund_list if any
	if ($reads_list ne "")
	{
		open(FILE, "<$reads_list") || die "Cannot open reads list file $reads_list\n";
		while(defined($line = <FILE>))
		{
			chomp($line);
			if (-e $line)
			{
				print "Located reads file [$line]\n";
				print TMPLOG "Located reads file [$line]\n";
				$reads_f[$readscount] = $line;
				$readscount++;
			}
			elsif ($line ne "")
			{
				print "Cannot find reads file [$line]. Stop.\n";
				exit;
			}
		}
		close(FILE);
	}
	if ($abund_list ne "")
	{
		open(FILE, "<$abund_list") || die "Cannot open abundance list file $abund_list\n";
		while(defined($line = <FILE>))
		{
			chomp($line);
			if (-e $line)
			{
				print "Located abundance file [$line]\n";
				print TMPLOG "Located abundance file [$line]\n";
				$abund_f[$abundcount] = $line;
				$abundcount++;
			}
			elsif ($line ne "")
			{
				print "Cannot find abundance file [$line]. Stop.\n";
				exit;
			}
		}
		close(FILE);
	}

	if ($readscount == 0 && $abundcount == 0)
	{
		print "Please input at least one abundance file or reads file. You may also specify a reads or abundance file list.\n";
		print $USAGE;
		close(TMPLOG);
		unlink("$k.log");
		exit;
	}

	$new_contig = checkContig($contig_f, $MIN_SEQ_LENGTH, "$out_f.contig.tmp", "$out_f.tooshort");
	if ($new_contig ne "")
	{
		$old_contig_name = $contig_f;
		$contig_f = $new_contig;
		$remodel_contig = 1;
	}

	$param = "";
	$j = 1;
	for ($i = 0; $i < $abundcount; $i++)
	{
		$line = "$contig_f.abund" . $j;
		checkAbundFile($abund_f[$i], $line);
		if ($j == 1)
		{
			$param = $param . " -abund " . $line;
		}
		else
		{
			$param = $param . " -abund$j " . $line;
		}
		$j++;
	}

	openLOG("$out_f.log");
	close(TMPLOG);
	open(TMPLOG, "<$k.log");
	while(defined($line = <TMPLOG>))
	{
		print $LOGOUT $line;
	}
	close(TMPLOG);
	unlink("$k.log");
	
	$k = 0;
	for ($i = 0; $i < $readscount; $i++)
	{
		writeLOG("Running Bowtie2 on reads file \[$reads_f[$i]\]...this may take a while...\n");
		# Run bowtie2 to find abundance information
		if ($k == 0)
		{
			if (!(-e "$out_f.idx.$FINISH"))
			{
				formatBowtie2($old_contig_name, $out_f, "$out_f.sam$i");
			}
			touch("$out_f.idx.$FINISH");
			$k = 1;
		}
		if (!(-e "$contig_f.reads.abund$j.$FINISH"))
		{
			runBowtie2($old_contig_name, $reads_f[$i], $out_f, "$out_f.sam$i");
			$j = $i + 1;
			$line = "$contig_f.reads.abund$j";
			getsam("$out_f.sam$i", $line);
			touch("$contig_f.reads.abund$j.$FINISH");
		}
		else
		{
			$j = $i + 1;
			$line = "$contig_f.reads.abund$j";
		}
		if ($j == 1)
		{
			$param = $param . " -abund " . $line;
		}
		else
		{
			$param = $param . " -abund$j " . $line;
		}
		$j++;
	}

	my @binarr = ();
	my $currbin;
	my $currout;
	my $currnum;
	my $maxnum;
	my @summarr = ();
	my @allsumarr = ();
	my @reclassifyarr = ();
	my @noclassarr = ();
	push(@binarr, $contig_f);
	# Push some dummy into result array for the original dataset
	push(@summarr, "dummy");
	push(@allsumarr, "dummy");
	push(@reclassifyarr, 1);
	$currnum = 0;
	$maxnum = 1;

	# Run HMMER3 to identify seed contigs
	writeLOG("Searching against $MARKERNUM marker genes to find starting seed contigs for [$contig_name]...\n");
	if (!(-e "$contig_f.hmmout.$FINISH"))
	{
		if (-e "$contig_f.frag.faa")
		{
			unlink("$contig_f.frag.faa");
		}
		getHMM($contig_f, "$contig_f.hmmout");
		touch("$contig_f.hmmout.$FINISH");
	}

	while ($currnum < $maxnum)
	{
		$currbin = $binarr[$currnum];

		if ($currbin eq $contig_f)
		{
			$currout = $out_f;
		}
		else
		{
			$currout = substr($currbin, 0, length($currbin) - 6) . ".out";
		}
		if ($MAX_MARKER == 1)
		{
			$i = gethmmmarker("$contig_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", 1);
		}
		else
		{
			$i = gethmmmarker("$contig_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed");
		}
		# Check if current file exceeds the current file processing recursion limitation
		# Check for how many "####.out", in which # represent digits
		@tmparr = $currbin  =~ /[0-9]{4}.out/g;
		$j = scalar @tmparr;

		if ($currbin eq $contig_f && $i == -1)
		{
			writeLOG("Try harder to dig out marker genes from contigs.\n");
			$i = gethmmmarker("$contig_f.hmmout", $currbin, $MIN_SEQ_LENGTH, "$currout.seed", 1);
			if ($i == -1)
			{
				writeLOG("Marker gene search reveals that the dataset cannot be binned (the medium of marker gene number <= 1). Program stop.\n");
				closeLOG("$out_f.log");
				exit(-1);
			}
		}
		elsif ($currbin ne $contig_f && $j >= $RECURSION_MAX)
		{
			$currnum++;
			next;
		}
		elsif ($currbin ne $contig_f && $i == -1)
		{
			$currnum++;
			next;
		}
		#elsif ($j > $RECURSION_MAX)
		#{
		#	$currnum++;
		#	next;
		#}
		elsif ($i != -1)
		{
			$reclassifyarr[$currnum] = 1;
		}

		# Running MaxBin
		writeLOG("Done data collection. Running MaxBin...\n");
		$param2 = "";
		if ($maxem != -1)
		{
			$param2 = "-max_run $maxem ";
		}
		if ($prob_threshold != -1)
		{
			$param2 = $param2 . "-prob_threshold $prob_threshold ";
		}
		if ($verbose == 1)
		{
			$param2 = $param2 . "-verbose ";
		}
		if ($THREADNUM > 1)
		{
			$param2 = $param2 . "-thread $THREADNUM";
		}
		$cmd = "$MAXBIN -fasta $currbin $param -seed $currout.seed -out $currout -min_contig_length $MIN_SEQ_LENGTH $param2";
		writeLOG("Command: $cmd\n");
		if (!(-e "$currout.$FINISH"))
		{
			system($cmd);
		}

		if (checkResult("$currout.summary") == -1)
		{
			if ($currbin eq $contig_f)
			{
				writeLOG("Error encountered while running core MaxBin program. Error recorded in $currout.log.\nProgram Stop.\n");
				exit(-1);
			}
		}
		push(@noclassarr, $currout);

		# Read in summary file and putting bins into stack
		open(FILE, "<$currout.summary");
		while(defined($line = <FILE>))
		{
			chomp($line);
			if ($line =~ /^Bin \[([A-Za-z0-9._\\\/\@\!\|\#\$\%\^\?\<\>\[\]\{\}\(\)\+\-]+)\]\t([0-9.\t]+)/)
			{
				push(@binarr, $1);
				@arr = split(/\t/, $2);
				push(@summarr, $arr[0]);
				push(@allsumarr, $2);
				push(@reclassifyarr, 0);
				$maxnum++;
			}
		}
		close(FILE);

		$currnum++;
		touch("$currout.$FINISH");
		push(@finisharr, "$currout.$FINISH");
	}

	# Remove progress recording files that ends with $FINISH
	unlink("$out_f.idx.$FINISH");
	for ($j = 1; $j <= $readscount; $j++)
	{
		unlink("$contig_f.reads.abund$j.$FINISH");
	}
	unlink("$contig_f.hmmout.$FINISH");
	foreach $i (@finisharr)
	{
		unlink($i);
	}

	# Put everything in order

	# Remove re-classified bins
	my @tmpgenomesize = ();
	my @tmpgc = ();
	my @genomesize = ();
	my @gc = ();
	my @bin_to_noclass = ();
	for ($k = 0; $k < $maxnum; $k++)
	{
		($tmpgenomesize[$k], $tmpgc[$k]) = getBinInfo($binarr[$k]);
	}
	for ($k = 0; $k < $maxnum; $k++)
	{
		if ($reclassifyarr[$k] == 1 || $tmpgenomesize[$k] < $MIN_BIN_SIZE)
		{
			if ($reclassifyarr[$k] == 1 && $binarr[$k] ne $contig_f)
			{
				unlink($binarr[$k]);
			}
			elsif ($reclassifyarr[$k] == 0 && $tmpgenomesize[$k] < $MIN_BIN_SIZE)
			{
				push(@bin_to_noclass, $binarr[$k]);
			}
			splice(@reclassifyarr, $k, 1);
			splice(@binarr, $k, 1);
			splice(@allsumarr, $k, 1);
			splice(@summarr, $k, 1);
			splice(@tmpgenomesize, $k, 1);
			splice(@tmpgc, $k, 1);
			$k--;
			$maxnum--;
		}
	}

	# First sort @summarr and store the index in another array
	my @sortarr = sort {$b <=> $a} @summarr;

	#my %indexhash;
	#@indexhash{@sortarr} = (0..$#summarr);

	my @sorted;
	my @indexarr;

	for ($i = 0; $i < $maxnum; $i++)
	{
		for ($j = 0; $j < $maxnum; $j++)
		{
			if (!(defined($sorted[$j])) && $summarr[$i] == $sortarr[$j])
			{
				$indexarr[$i] = $j;
				$sorted[$j] = 1;
				last;
			}
		}
	}
	@sortarr = ();

	my @filearr;
	for ($k = 0; $k < $maxnum; $k++)
	{
		#$i = $indexhash{$summarr[$k]};
		$i = $indexarr[$k];
		$j = getBinNum($i + 1, $maxnum);
		# rename original fasta file
		rename ($binarr[$k], "$out_f.$j.fasta.tmp");
		$filearr[$i] = "$out_f.$j.fasta.tmp";
		$genomesize[$i] = $tmpgenomesize[$k];
		$gc[$i] = $tmpgc[$k];
		$sortarr[$i] = $allsumarr[$k];
	}
	# Write summary
	open(TMPSUM, ">$out_f.tmp.summary");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$j = getBinNum($i + 1, $maxnum);
		print TMPSUM "Bin \[$out_f.$j.fasta\]\t$sortarr[$i]\n";
	}
	close(TMPSUM);

	# Collect all no-classes and log file
	open(OUT, ">$out_f.tmp.noclass");
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		open(FILE, "<$noclassarr[$i].noclass");
		while(defined($line = <FILE>))
		{
			print OUT $line;
		}
		close(FILE);

		open(FILE, "<$noclassarr[$i].log");
		while(defined($line = <FILE>))
		{
			writeLOG($line, 0);
		}
		writeLOG("\n");
		close(FILE);
	}
	foreach $currbin (@bin_to_noclass)
	{
		if (-e $currbin)
		{
			open(FILE, "<$currbin");
			while(defined($line = <FILE>))
			{
				print OUT $line;
			}
			close(FILE);
			unlink($currbin);
		}
		else
		{
			writeLOG("File $currbin not found.\n");
		}
	}
	close(OUT);

	# delete all normal files
	#for ($i = 0; $i < $maxnum; $i++)
	#{
	#	if ($binarr[$i] ne $contig_f && (-e $binarr[$i]))
	#	{
	#		unlink($binarr[$i]);
	#	}
	#}
	$j = @noclassarr;
	for ($i = 0; $i < $j; $i++)
	{
		unlink("$noclassarr[$i].log");
		unlink("$noclassarr[$i].summary");
		unlink("$noclassarr[$i].noclass");
		#unlink("$noclassarr[$i].prob_dist");
		#unlink("$noclassarr[$i].dist");
		unlink("$noclassarr[$i].seed");
	}

	# rename all tmp files to normal files
	rename("$out_f.tmp.summary", "$out_f.summary");
	rename("$out_f.tmp.noclass", "$out_f.noclass");
	for ($i = 0; $i < $maxnum; $i++)
	{
		$currbin = substr($filearr[$i], 0, length($filearr[$i]) - 4);
		rename($filearr[$i], $currbin);
		$filearr[$i] = $currbin;
	}

	# Separate reads for reassembly
	my %lenhash;
	my $n50_before;
	my $n50_after;
	if ($REASSEMBLY == 1)
	{
		if ($readscount <= 0)
		{
			writeLOG("Cannot perform reassembly since reads file is not provided.\n");
			$REASSEMBLY = 0;
			last;
		}
		writeLOG("Performing reassembly. Reads file found.\n");
		writeLOG("Separating reads according to the bins...\n");
		for ($i = 0; $i < $readscount; $i++)
		{
			$isfastq = checkFastq($reads_f[$i]);
			sepReads($reads_f[$i], $i + 1, $out_f, "$out_f.sam$i", $isfastq);
		}
		mkdir("$out_f.tmp");
		mkdir("$out_f.reassem");
		open(N50OUT, ">$out_f.reassem/N50.txt");
		open(FILE, "<$out_f.summary");
		print N50OUT "Bin\tN50 before reassem\tN50 after reassem\tN50 improvement\n";
		while(defined($line = <FILE>))
		{
			if ($line =~ /Bin \[$out_f.([0-9]+).fasta\]\t([0-9.]+)/)
			{
				$j = $1;
				writeLOG("Reassembling bin $j\n");
				$cmd = "$IDBA_UD -r $out_f.reads.$j -o $out_f.tmp --pre_correction --num_threads $THREADNUM 1>$out_f.idba.out 2>$out_f.idba.err";
				writeLOG("Command: $cmd\n");
				system($cmd);
				if (checkResult("$out_f.tmp\/scaffold.fa") == -1)
				{
					writeLOG("Error occurs while reassembling Bin $out_f.$j.fasta.\nError logged in $out_f.idba.out and $out_f.idba.err.\nProgram stop.\n");
					closeLOG("$out_f.log");
					exit;
				}

				writeLOG("Placing reassembled contigs into $out_f.reassem\n");
				$cmd = "cp $out_f.tmp\/scaffold.fa $out_f.reassem/$out_f.$j.fasta";
				writeLOG("$cmd\n");
				system($cmd);

				$n50_before = getN50("$out_f.$j.fasta");
				$n50_after = getN50("$out_f.reassem/$out_f.$j.fasta");
				printf N50OUT "$out_f.$j.fasta\t$n50_before\t$n50_after\t%0.2f%%\n", (($n50_after - $n50_before) / $n50_before) * 100;

				#unlink("$out_f.reads.$1");
				$cmd = "mv $out_f.reads.$j $out_f.reassem/";
				system($cmd);
			}
		}
		close(FILE);
		close(N50OUT);
		$cmd = "mv $out_f.reads.noclass $out_f.reassem/";
		system($cmd);
	}

	# Count marker genes for bins
	my $completearr;
=cut
	if ($REASSEMBLY == 1)
	{
		$cmd = "cat $out_f.*.fasta > $out_f.all.fasta";
		print "$cmd\n";
		system($cmd);
		getHMM("$out_f.all.fasta", "$contig_f.hmmout", "cuttc");
		unlink("$out_f.all.fasta");
	}
	else
	{
		#getHMM($contig_f, "$contig_f.hmmout", "cuttc");
	}
=cut
	$completearr = countmarker("$contig_f.hmmout", $out_f, $outname, $MARKERHMM, "$out_f.marker");
	$i = 0;

	# Re-calculate the genome sizes and gc content after reassembly
	if ($REASSEMBLY == 1)
	{
		foreach $currbin (@filearr)
		{
			($genomesize[$i], $gc[$i]) = getBinInfo($currbin);
			$i++;
		}
	}

	# Re-read summary file and write other info: completeness, Genome size, GC
	open(FILE, "<$out_f.summary");
	open(SUMOUT, ">$out_f.summary.tmp");
	$j = $readscount + $abundcount;
	if ($j == 1)
	{
		print SUMOUT "Bin name\tAbundance\tCompleteness\tGenome size\tGC content\n";
	}
	else
	{
		print SUMOUT "Bin name\tCompleteness\tGenome size\tGC content\n";
		open(SUMABUND, ">$out_f.abundance");
		print SUMABUND "Bin name";
		for ($i = 0; $i < $abundcount; $i++)
		{
			print SUMABUND "\t$abund_f[$i]";
		}
		for ($i = 0; $i < $readscount; $i++)
		{
			print SUMABUND "\t$reads_f[$i]";
		}
		print SUMABUND "\n";
	}
	$i = 0;
	while(defined($line = <FILE>))
	{
		if ($line =~ /^Bin \[$out_f.([0-9]+).fasta\]\t([0-9.\t]+)/)
		{
			if ($j == 1)
			{
				printf SUMOUT "$outname.$1.fasta\t%0.2f\t%0.1f%%\t$genomesize[$i]\t%0.1f\n", $2, $$completearr[$i] * 100, $gc[$i];
			}
			else
			{
				print SUMOUT "$outname.$1.fasta";
				printf SUMOUT "\t%0.1f%%\t$genomesize[$i]\t%0.1f\n", $$completearr[$i] * 100, $gc[$i];
				print SUMABUND "$outname.$1.fasta";
				@arr = split(/\t/, $2);
				foreach $tmp (@arr)
				{
					printf SUMABUND "\t%0.2f", $tmp;
				}
				print SUMABUND "\n";
			}
			$i++;
		}
	}
	close(FILE);
	close(SUMOUT);
	if ($j > 1)
	{
		close(SUMABUND);
	}
	rename("$out_f.summary.tmp", "$out_f.summary");

	if ($RPLOT == 1)
	{
		$cmd = "$RSCRIPT $HEATMAP_R $out_f.marker $out_f.marker.pdf";
		print "$cmd\n";
		system($cmd);
	}

	# Get markers for each bin
	listmarker_bybin("$contig_f.hmmout", $out_f, "$contig_f.frag.faa", $MARKERHMM, $out_f);
	if ($outdir ne "")
	{
		chdir($outdir);
	}
	$cmd = "tar zcvf $outname.marker_of_each_bin.tar.gz $outname.*.marker.fasta";
	system($cmd);
	$cmd = "rm $outname.*.marker.fasta";
	system($cmd);
	if ($outdir ne "")
	{
		chdir($currdir);
	}

	if ($preserve == 0)
	{
		writeLOG("Deleting intermediate files.\n");
		for ($i = 0; $i < $readscount; $i++)
		{
			unlink("$out_f.sam$i");
		}
		$cmd = "rm -rf $out_f.tmp";
		system($cmd);
		#unlink("$contig_f.prob_dist");
		#unlink("$contig_f.seed.m8");
		unlink("$contig_f.seed");
		for ($i = 0; $i < $abundcount; $i++)
		{
			$j = $i + 1;
			$line = "$contig_f.abund" . $j;
			unlink($line);
		}
		for ($i = 0; $i < $readscount; $i++)
		{
			$j = $i + 1;
			$line = "$contig_f.reads.abund$j";
			$tmp = $out_f . ".abund" . $j;
			rename($line, $tmp);
		}
		unlink("$contig_f.hmmout");
		if (-e "$contig_f.frag.faa")
		{
			unlink("$contig_f.frag.faa");
		}
	}

	unlink("$out_f.idx.1.bt2");
	unlink("$out_f.idx.2.bt2");
	unlink("$out_f.idx.3.bt2");
	unlink("$out_f.idx.4.bt2");
	unlink("$out_f.idx.rev.1.bt2");
	unlink("$out_f.idx.rev.2.bt2");
	if ($remodel_contig == 1)
	{
		unlink($contig_f);
	}

	# Write post instruction to users
	writeLOG("\n\n========== Job finished ==========\nYielded $maxnum bins for contig (scaffold) file $contig_name\n\n");
	writeLOG("Here are the output files for this run.\nPlease refer to the README file for further details.\n\n");
	$i = $maxnum;
	$j = getBinNum($i, $maxnum);
	$tmp = getBinNum(1, $maxnum);
	writeLOG("Summary file: $out_f.summary\n");
	$k = $readscount + $abundcount;
	if ($k > 1)
	{
		writeLOG("Genome abundance info file: $out_f.abundance\n");
	}
	writeLOG("Marker counts: $out_f.marker\nMarker genes for each bin: $out_f.marker_of_each_gene.tar.gz\nBin files: $out_f.$tmp.fasta - $out_f.$j.fasta\nUnbinned sequences: $out_f.noclass\n");

	if ($RPLOT == 1 && -e "$out_f.marker.pdf")
	{
		writeLOG("Marker plot: $out_f.marker.pdf\n");
	}
	if ($readscount > 0)
	{
		writeLOG("\n");
		for ($i = 0; $i < $readscount; $i++)
		{
			$j = $i + 1;
			$tmp = $out_f . ".abund" . $j;
			writeLOG("Store abundance information of reads file \[$reads_f[$i]\] in \[$tmp\].\n");
		}
	}
	writeLOG("\n\n========== Elapsed Time ==========\n");
	$endtime = time();
	$line = getElapsedTime($endtime - $starttime);
	writeLOG("$line\n");

	closeLOG("$out_f.log");
}

sub formatBowtie2
{
	my $contig_f = $_[0];
	my $outheader = $_[1];
	my $out_f = $_[2];
	my $cmd = "$BOWTIE2BUILD $contig_f $outheader.idx 1>$out_f.bowtie2build.out 2>$out_f.bowtie2build.err";
	system($cmd);
}

sub runBowtie2
{
	my $contig_f = $_[0];
	my $reads_f = $_[1];
	my $outheader = $_[2];
	my $out_f = $_[3];
	my $cmd;
	my $isfastq;
	my $isgz;

	# Run bowtie2 to find abundance information
	#$cmd = "$BOWTIE2BUILD $contig_f $outheader.idx 1>$out_f.bowtie2build.out 2>$out_f.bowtie2build.err";
	#system($cmd);
	$isfastq = checkFastq($reads_f);
	if ($isfastq eq "fq" || $isfastq eq "fqgz")
	{
		$cmd = "$BOWTIE2 -p $THREADNUM -x $outheader.idx -U $reads_f -S $out_f 1>$out_f.bowtie2.out 2>$out_f.bowtie2.err";
	}
	elsif ($isfastq eq "fa" || $isfastq eq "fagz")
	{
		$cmd = "$BOWTIE2 -f -p $THREADNUM -x $outheader.idx -U $reads_f -S $out_f 1>$out_f.bowtie2.out 2>$out_f.bowtie2.err";
	}
	system($cmd);
	if (checkResult("$out_f") == -1)
	{
		print "Error running Bowtie2. Bowtie2 log recorded in $out_f.bowtie2build.out, $out_f.bowtie2build.err, $out_f.bowtie2.out, and $out_f.bowtie2.err\nProgram stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$out_f.bowtie2build.out");
		unlink("$out_f.bowtie2build.err");
		unlink("$out_f.bowtie2.out");
		unlink("$out_f.bowtie2.err");
	}
}

sub getHMM
{
	my $contig_f = $_[0];
	my $out_f = $_[1];
	my $cutmethod = $_[2];
	my $cmd;

	if (!(-e "$contig_f.frag.faa"))
	{
		print "Running FragGeneScan....\n";
		$cmd = "$RUNFRAG -genome=$contig_f -out=$contig_f.frag -complete=0 -train=complete -thread=$THREADNUM 1>$contig_f.frag.out 2>$contig_f.frag.err";
		system($cmd);
	}
	if (checkResult("$contig_f.frag.faa") == -1)
	{
		print "Error running FragGeneScan. Output recorded in $contig_f.frag.out and $contig_f.frag.err.\n";
		print "=== Please make sure that you are running FragGeneScan 1.18 or above. ===\n";
		print "=== There are known bugs in v1.17 and before that will crash FragGeneScan program. ===\n";
		exit(-1);
	}
	else
	{
		unlink("$contig_f.frag.out");
		unlink("$contig_f.frag.err");
	}

=cut # Prodigal gene prediction
	print "Running Prodigal....\n";
	$cmd = "$PRODIGAL -a $contig_f.faa -i $contig_f -m -o $contig_f.prodigal.tmp -p meta -q 1>$contig_f.prodigal.out 2>$contig_f.prodigal.err";
	system($cmd);
	if (checkResult("$contig_f.faa") == -1)
	{
		print "Error running Prodigal. Output recorded in $contig_f.prodigal.out and $contig_f.prodigal.err.\nProgram Stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$contig_f.prodigal.out");
		unlink("$contig_f.prodigal.err");
		unlink("$contig_f.prodigal.tmp");
	}
=cut

	print "Running HMMER hmmsearch....\n";
	if ($MARKERNUM == 107)
	{
		#$cmd = "$HMMSEARCH --tblout $out_f --cut_tc --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
		#if (defined($cutmethod) && $cutmethod eq "cuttc")
		#{
			$cmd = "$HMMSEARCH --domtblout $out_f --cut_tc --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
		#}
		#else
		#{
		#	$cmd = "$HMMSEARCH --domtblout $out_f -E 1e-5 --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
		#}
	}
	else
	{
		$cmd = "$HMMSEARCH --domtblout $out_f -E 1e-3 --cpu $THREADNUM $MARKERHMM $contig_f.frag.faa 1>$out_f.out 2>$out_f.err";
	}
	system($cmd);
	if (checkResult("$out_f") == -1)
	{
		print "Error running Hmmer3. Output recorded in $out_f.out and $out_f.err.\nProgram Stop.\n";
		exit(-1);
	}
	else
	{
		unlink("$out_f.out");
		unlink("$out_f.err");
	}
	#unlink("$contig_f.frag.faa");
	unlink("$contig_f.frag.ffn");
	unlink("$contig_f.frag.out");
	unlink("$contig_f.frag.gff");
}

sub getBinNum
{
	my $snum = $_[0];
	my $smax = $_[1];
	my $tmp = "";
	my $ret = "";
	if ($smax > 1000)
	{
		if ($snum < 10)
		{
			$tmp = "000" . $snum;
		}
		elsif ($snum < 100)
		{
			$tmp = "00" . $snum;
		}
		elsif ($snum < 1000)
		{
			$tmp = "0" . $snum;
		}
		else
		{
			$tmp = $snum;
		}
	}
	else
	{
		if ($snum < 10)
		{
			$tmp = "00" . $snum;
		}
		elsif ($snum < 100)
		{
			$tmp = "0" . $snum;
		}
		else
		{
			$tmp = $snum;
		}
	}

	#if ($smax > 1000 && $snum < 1000)
	#{
	#	$ret = "0" . $tmp;
	#}
	#else
	#{
		$ret = $tmp;
	#}
	return $ret;
}

sub checkResult
{
	my $s = -s $_[0];
	if (-e $_[0] && $s > 0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

sub checkContig
{
	my $contig_f = $_[0];
	my $min_length = $_[1];
	my $faout = $_[2];
	my $below_f = $_[3];
	my $i;
	my $len;
	my $faline;
	my $header;
	my $seq;
	my $gz = "";
	my $FA_LIMIT = 128;
	my $FASTA_LEN = 70;

	if ($contig_f =~ /.gz$/)
	{
		$gz = $contig_f . ".fa";
		$i = "gunzip -c $contig_f > $gz";
		system($i);
		$contig_f = $gz;
	}

	open(FILE, "<$contig_f") || die "Cannot open file $contig_f\n";
	$i = 0;

	# Remodel input contig file
	seek(FILE, 0, 0);
	$header = "";
	#$faout = $contig_f . ".tmp";
	open(FAOUT, ">$faout") || die "Cannot write into specified output file/directory. Please check your settings or disk space.\n";
	open(FABELOW, ">$below_f");
	while(defined($faline = <FILE>))
	{
		chomp($faline);
		if ($faline =~ /^>/)
		{
			if ($header ne "" && ($len = length($seq)) >= $min_length)
			{
				print FAOUT "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FAOUT substr($seq, $i);
						print FAOUT "\n";
						$i = $len;
					}
					else
					{
						print FAOUT substr($seq, $i, $FASTA_LEN);
						print FAOUT "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			elsif ($header ne "" && ($len = length($seq)) < $min_length)
			{
				print FABELOW "$header\n";
				$i = 0;
				while ($i < $len)
				{
					if ($i + $FASTA_LEN > $len)
					{
						print FABELOW substr($seq, $i);
						print FABELOW "\n";
						$i = $len;
					}
					else
					{
						print FABELOW substr($seq, $i, $FASTA_LEN);
						print FABELOW "\n";
						$i += $FASTA_LEN;
					}
				}
			}
			$header = $faline;
			$seq = "";
		}
		else
		{
			$seq = $seq . $faline;
		}
	}
	$len = length($seq);
	if ($len >= $min_length)
	{
		print FAOUT "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FAOUT substr($seq, $i);
				print FAOUT "\n";
				$i = $len;
			}
			else
			{
				print FAOUT substr($seq, $i, $FASTA_LEN);
				print FAOUT "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	else
	{
		print FABELOW "$header\n";
		$i = 0;
		while ($i < $len)
		{
			if ($i + $FASTA_LEN > $len)
			{
				print FABELOW substr($seq, $i);
				print FABELOW "\n";
				$i = $len;
			}
			else
			{
				print FABELOW substr($seq, $i, $FASTA_LEN);
				print FABELOW "\n";
				$i += $FASTA_LEN;
			}
		}
	}
	close(FILE);
	close(FAOUT);
	close(FABELOW);

	if ($gz ne "")
	{
		unlink($gz);
	}

	return $faout;
}

# (Genome size, GC content) = getBinInfo(fasta file)
sub getBinInfo
{
	my $all = 0;
	my $gc = 0;
	my $i;
	my $line;
	open(FAFILE, "<$_[0]");
	while(defined($line = <FAFILE>))
	{
		if ($line !~ /^>/)
		{
			chomp($line);
			$i = ($line =~ tr/[cgCG]//);
			$gc += $i;
			$i = ($line =~ tr/[atcgATCG]//);
			$all += $i;
		}
	}
	close(FAFILE);

	if ($all == 0)
	{
		$gc = 0;
	}
	else
	{
		$gc = $gc / $all * 100;
	}
	return ($all, $gc);
}

# (Read and Check program directory)
sub checkProgram #Modified by Fernando Puente-Sánchez, 27-IV-2018. Settings file parsing is overriden, instead, program locations are provided above. Checking will still occur.
{
	my $line;
	my $tmpstr;
	my $tmpname =  "tmp_" . time();

#	open(FILE, "<$Bin\/$SETTING_FILE");
#	while(defined($line = <FILE>))
#	{
#		chomp($line);
#		if ($line =~ /\[([A-Za-z0-9_]+)\] ([A-Za-z0-9._\(\)\[\]\{\}\|\$\!\=\-\+\\\/]+)/)
#		{
#			if ($1 eq "FragGeneScan")
#			{
#				if (-d $2 && -e "$2\/$RUNFRAG")
#				{
#					$RUNFRAG = $2 . "\/" . $RUNFRAG;
#				}
#			}
#			elsif ($1 eq "Bowtie2")
#			{
#				if (-d $2 && -e "$2\/$BOWTIE2" && -e "$2\/$BOWTIE2BUILD")
#				{
#					$BOWTIE2 = $2 . "\/" . $BOWTIE2;
#					$BOWTIE2BUILD = $2 . "\/" . $BOWTIE2BUILD;
#				}
#			}
#			elsif ($1 eq "HMMER3")
#			{
#				if (-d $2 && -e "$2\/$HMMSEARCH")
#				{
#					$HMMSEARCH = $2 . "\/" . $HMMSEARCH;
#				}
#			}
#			elsif ($1 eq "IDBA_UD")
#			{
#				if (-d $2 && -e "$2\/$IDBA_UD")
#				{
#					$IDBA_UD = $2 . "\/" . $IDBA_UD;
#				}
#			}
#		}
#	}
#	close(FILE);

	# Check program
	# FragGeneScan
	$line = "$RUNFRAG 1>$tmpname 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /FragGeneScan/)
	{
		print "Cannot run FragGeneScan. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	# Bowtie2
	$line = "$BOWTIE2 1>/dev/null 2>$tmpname";
	system($line);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /bowtie/)
	{
		print "Cannot run Bowtie2. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	# HMMER3
	$line = "$HMMSEARCH 1>$tmpname 2>/dev/null";
	system($line);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /hmmsearch/)
	{
		print "Cannot run HMMER3. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	# IDBA_UD
	$line = "$IDBA_UD 1>/dev/null 2>$tmpname";
	system($line);
	$tmpstr = "";
	open(FILE, "<$tmpname");
	while(<FILE>)
	{
		$tmpstr .= $_;
	}
	if ($tmpstr !~ /IDBA\-UD/)
	{
		print "Cannot run IDBA-UD. Please indicate the file directory in \'setting\' file.\n";
		exit;
	}

	unlink("$tmpname");
}

sub checkFastq
{
	my $input_f = $_[0];
	my $fastq_line;
	my $fastq = "";
	if ($input_f =~ /q.gz$/)
	{
		return "fqgz";
	}
	if ($input_f =~ /q.bz2$/)
	{
		return "fqgz";
	}
	elsif ($input_f =~ /a.gz$/)
	{
		return "fagz";
	}
	elsif ($input_f =~ /a.bz2$/)
	{
		return "fagz";
	}
	open(CHECK, "<$input_f") || die "Reads file [$input_f] not found. Check path again.\n";
	while(defined($fastq_line = <CHECK>))
	{
		chomp($fastq_line);
		if ($fastq_line =~ /^>/)
		{
			$fastq = "fa";
		}
		elsif ($fastq_line =~ /^\@/)
		{
			$fastq = "fq";
		}
		elsif ($fastq_line ne "")
		{
			print "Error header found in reads file.\n======\n$fastq_line\n======\nMaxBin stop.\n";
			exit;
		}
		if ($fastq ne "")
		{
			last;
		}
	}
	close(CHECK);
	return($fastq);
}

sub checkAbundFile
{
	my $input_f = $_[0];
	my $out_f = $_[1];
	my $line;

	if ($input_f eq $out_f)
	{
		return;
	}

	#Check header
	open(ABUNDIN, "<$input_f") || die "Cannot open abund file $input_f\n";
	open(ABUNDOUT, ">$out_f") || die "Cannot open tmp abund file $out_f\n";
	while(defined($line = <ABUNDIN>))
	{
		if ($line =~ /^>/)
		{
			print ABUNDOUT substr($line, 1);
		}
		else
		{
			print ABUNDOUT $line;
		}
	}
	close(ABUNDOUT);
	close(ABUNDIN);
}

sub getElapsedTime
{
	my $input = $_[0];
	my $hour;
	my $min;
	my $sec;
	my $str;

	$sec = $input % 60;
	$input = int($input / 60);
	$min = $input % 60;
	$input = int($input / 60);
	$hour = $input;

	$str = $hour . " hours " . $min . " minutes and " . $sec . " seconds.\n";
	return $str;
}

sub openLOG
{
	my $log_f = $_[0];
	open($LOGOUT, ">$log_f.tmp");
}

# Currently the writeLOG procedure will NOT write to STDOUT only if there is a second parameter, say 1. Will be revised later.
sub writeLOG
{
	my $msg = $_[0];
	my $is_print = $_[1];
	if (!(defined($is_print)))
	{
		print $msg;
	}
	print $LOGOUT $msg;
}

sub closeLOG
{
	my $log_f = $_[0];
	close($LOGOUT);
	rename("$log_f.tmp", $log_f);
}

sub getN50
{
	my $f = $_[0];
	my $line;
	my $N50_len;
	my $count;
	my @len;
	my $i;
	my $sum = 0;

	# Read contig length
	$count = -1;
	open(N50FILE, "<$f") || return -1;
	while(defined($line = <N50FILE>))
	{
		if ($line =~ /^>/)
		{
			$count++;
		}
		else
		{
			chomp($line);
			$i = length($line);
			$len[$count] += $i;
			$sum += $i;
		}
	}
	close(N50FILE);
	$count++;

	my @sort = sort{$b <=> $a} @len;
	my $halfcount = 0;
	for ($i = 0; $i < $count; $i++)
	{
		$halfcount = $halfcount + $sort[$i];
		if ($halfcount >= $sum / 2)
		{
			last;
		}
	}
	$N50_len = $sort[$i];
	return($N50_len);
}

sub touch
{
	my $cmd = "touch $_[0]";
	system($cmd);
}

