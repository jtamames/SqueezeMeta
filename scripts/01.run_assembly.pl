#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades). Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$megahit_soft,$assembler_options,$extassembly,$contigid,$numthreads,$spades_soft,$flye_soft,$prinseq_soft,$mappingfile,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$scriptdir,$singletons,$methodsfile,$syslogfile,$norename);

my($seqformat,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name,%extassemblies);

if(-e "$datapath/raw_fastq/par1.fastq.gz") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fasta.gz") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta.gz"; $par2name="$datapath/raw_fastq/par2.fasta.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fastq") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
elsif(-e "$datapath/raw_fastq/par1.fasta") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta"; $par2name="$datapath/raw_fastq/par2.fasta"; }
else { die "Can't find read files in $datapath/raw_fastq\n"; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

open(infile1,$mappingfile);  #-- To check for extassemblies in sequential mode
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	$_=~s/\r//g;
	my($sample,$file,$iden,$mapreq)=split(/\t/,$_);
	if($mapreq) {
		my @k=split(/\,\s?/,$mapreq);
		foreach my $d(@k) {
			$d=~s/\"//g;
			if($d=~/extassembly\=(.*)/) { $extassemblies{$sample}=$1; }
			}
		}
	}
close infile1;


if($extassembly) {
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
	$outassembly=$extassembly; 
	}
elsif($extassemblies{$projectname}) {
	$extassembly=$extassemblies{$projectname};
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
	$outassembly=$extassembly; 
	}
	
else {

	#-- Checks the assembler

	if($assembler=~/megahit/i) { 
		system("rm -r $datapath/megahit > /dev/null 2>&1"); 
		$outassembly="$datapath/megahit/final.contigs.fa";
		if(-e $par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }
		else {  $command="$megahit_soft $assembler_options -r $par1name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }  #-- Support for single reads
		print outmet "Assembly was done using Megahit (Li et al 2015, Bioinformatics 31(10):1674-6)\n";
		}
	elsif($assembler=~/spades/i) { # works for spades and rnaspades
		system("rm -r $datapath/spades > /dev/null 2>&1");
	        if($assembler eq "spades") { $outassembly="$datapath/spades/contigs.fasta"; }
		elsif($assembler eq "rnaspades") { $outassembly="$datapath/spades/transcripts.fasta"; }
		my $command_start;
		if($assembler eq 'spades') { $command_start="$spades_soft --meta"; }
		elsif($assembler eq 'rnaspades') { $command_start="$spades_soft --rna" }
		if(-e $par2name) { $command="$command_start --pe1-1 $par1name --pe1-2 $par2name -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile 2>&1"; }
		else { $command="$command_start --s1 $par1name  -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile"; } #-- Support for single reads
		print outmet "Assembly was done using SPAdes (Bankevich et al 2012, J Comp Biol 19(5):455-77)\n";
		}
	elsif($assembler=~/canu/i) {
                system("rm -r $datapath/canu > /dev/null 2>&1");
		$outassembly="$datapath/canu/contigs.fasta";
		if($canumem eq "NF") {
			print "  Setting available memory for Canu\n";
			my %mem=get_mem_info;
			my $ram=$mem{"MemAvailable"}/(1024*1024);
			my $ramstr=sprintf('%.2f',$ram);
			$canumem=sprintf('%.0f',int($ram));
			$canumem*=0.8;
			print "  AVAILABLE (free) RAM memory: $ramstr Gb. We will set canu to use $canumem Gb.\n  You can override this setting using the -canumem option when calling SqueezeMeta.pl\n";
			print outsyslog "canumem set to $canumem (Free Mem $ramstr Gb)\n";
			}
     	   	$command="$canu_soft  -p $project -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=$canumem oeaMemory=$canumem batMemory=$canumem mhapThreads=$numthreads mmapThreads=$numthreads ovlThreads=$numthreads ovbThreads=$numthreads ovsThreads=$numthreads corThreads=$numthreads oeaThreads=$numthreads redThreads=$numthreads batThreads=$numthreads gfaThreads=$numthreads merylThreads=$numthreads $assembler_options -nanopore-raw  $par1name > $syslogfile 2>&1; "; 
	   	$command.="mv $datapath/canu/$project.contigs.fasta $outassembly"; 
 		print outmet "Assembly was done using Canu (Koren et al 2017, Genome Res 27(5):722-36)\n";
     	  }
	elsif($assembler=~/flye/i) {
                system("rm -r $datapath/flye > /dev/null 2>&1");
                $outassembly="$datapath/flye/contigs.fasta";
                $command="$flye_soft $assembler_options -o $datapath/flye --plasmids --meta --genome-size 2.6g --threads $numthreads --nano-raw $par1name > $syslogfile 2>&1; "; 
                $command.="mv $datapath/flye/assembly.fasta $outassembly";
                print outmet "Assembly was done using Flye (Kolmogorov et al 2019, Nature Biotech 37, 540-546)\n";
          }


	else { die "Unrecognized assembler\n"; }
	
	#-- Run assembly

	print "  Running assembly with $assembler\n";
	print outsyslog "Running assembly with $assembler: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	}

#-- Run prinseq_lite for removing short contigs

if($extassembly) { system("cp $outassembly $contigsfna"); }
else {
	$command="$prinseq_soft -fasta $outassembly -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna > /dev/null 2>&1";
	print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen \n";
	print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n  ";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	}

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
print outsyslog "Running prinseq for contig statistics: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";

#-- Standardization of contig names

if(!$norename) {
	print "  Renaming contigs\n";
	open(infile1,$contigsfna) || die "Can't open $contigsfna\n";
	my $provcontigs="$tempdir/contigs.prov";
	open(outfile1,">$provcontigs") || die "Can't open $provcontigs for writing\n";
	my $cocount;
	if(!$contigid) { $contigid="$assembler"; }
	while(<infile1>) {
		chomp;
		next if !$_;
		if($_=~/^\>/) {
			$cocount++;
			my $newcontigname="$contigid\_$cocount";
			print outfile1 ">$newcontigname\n";
			}
		else { print outfile1 "$_\n"; }
		}
	close infile1;
	close outfile1;
	system("mv $provcontigs $contigsfna");
	}
	
#-- Counts length of the contigs (we will need it later)

my $numc;
print "  Counting length of contigs\n";
open(outfile2,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile2,$contigsfna) || die "Can't open $contigsfna\n";
while(<infile2>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
		$numc++;
		$thisname=$1;
		if($contigname) {
			$len=length $seq;
			print outfile2 "$contigname\t$len\n"; 
			}
		$seq="";
		$contigname=$thisname;
		}
	else { $seq.=$_;}
	}
close infile2;
if($contigname) { $len=length $seq; print outfile2 "$contigname\t$len\n"; }
close outfile2;
close outmet;

print "  Contigs stored in $contigsfna\n  Number of contigs: $numc\n";
#system("rm $datapath/raw_fastq/par1.$format.gz; rm $datapath/raw_fastq/par2.$format.gz");

close outsyslog;
