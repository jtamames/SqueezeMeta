#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades). Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;
use lib "."; 

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$megahit_soft,$assembler_options,$numthreads,$spades_soft,$prinseq_soft,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$contigsfna,$contigslen,$cleaning,$cleaningoptions);

my($seqformat,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name);

if(-e "$datapath/raw_fastq/par1.fastq.gz") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fasta.gz") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta.gz"; $par2name="$datapath/raw_fastq/par2.fasta.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fastq") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
elsif(-e "$datapath/raw_fastq/par1.fasta") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta"; $par2name="$datapath/raw_fastq/par2.fasta"; }
else { die "Cannot find read files in $datapath/raw_fastq\n"; }

#-- trimmomatic commands

if($cleaning) {
	my $orig1=$par1name;
	my $orig2=$par2name;
	$orig1=~s/\.fastq/\.original.fastq/;
	$orig1=~s/\.fasta/\.original.fasta/;
	$orig2=~s/\.fastq/\.original.fastq/;
	$orig2=~s/\.fasta/\.original.fasta/;
	my $tcommand="mv $par1name $orig1; mv $par2name $orig2";
	system $tcommand; 
	if(-e $orig2) { $trimmomatic_command="$trimmomatic_soft PE -threads $numthreads -phred33 $orig1 $orig2 $par1name $par1name.removed $par2name $par2name.removed $cleaningoptions"; }
	else { $trimmomatic_command="$trimmomatic_soft SE -threads $numthreads -phred33 $orig1 $par1name $cleaningoptions"; }

	if($cleaning) {
		print "Running trimmomatic: $trimmomatic_command\n";
		my $ecode = system $trimmomatic_command;
		if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
		}
	}

#-- Checks the assembler

if($assembler=~/megahit/i) { 
	$outassembly="$datapath/megahit/final.contigs.fa";
	if(-e $par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name -t $numthreads -o $datapath/megahit"; }
	else {  $command="$megahit_soft $assembler_options -r $par1name -t $numthreads -o $datapath/megahit"; }  #-- Support for single reads
	}
elsif($assembler=~/spades/i) { 
	$outassembly="$datapath/spades/contigs.fasta";
	if(-e $par2name) { $command="$spades_soft $assembler_options --meta --pe1-1 $par1name --pe1-2 $par2name -m 400 -k 21,33,55,77,99,127 -t $numthreads -o $datapath/spades"; }
	else { $command="$spades_soft $assembler_options --meta --s1 $par1name  -m 400 -k 21,33,55,77,99,127 -t $numthreads -o $datapath/spades"; } #-- Support for single reads
	}
elsif($assembler=~/canu/i) {
        $outassembly="$datapath/canu/contigs.fasta";
        $command="rm -r $datapath/canu; $canu_soft $assembler_options -p $project -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=$canumem oeaMemory=$canumem batMemory=$canumem mhapThreads=$numthreads mmapThreads=$numthreads ovlThreads=$numthreads ovbThreads=$numthreads ovsThreads=$numthreads corThreads=$numthreads oeaThreads=$numthreads redThreads=$numthreads batThreads=$numthreads gfaThreads=$numthreads merylThreads=$numthreads -nanopore-raw  $par1name;"; 
	$command.="mv $datapath/canu/$project.contigs.fasta $outassembly"; 
        }


else { die "Unrecognized assembler\n"; }

	
#-- Run assembly

print "Running assembly with $assembler: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }


#-- Run prinseq_lite for removing short contigs

$command="$prinseq_soft -fasta $outassembly -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna";
print "Running prinseq: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }


#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $resultpath/01.$project.stats";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }


#-- Counts length of the contigs (we will need it later)

print "Counting length of contigs\n";
open(outfile1,">$contigslen") || die;
open(infile1,$contigsfna) || die;
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
		$thisname=$1;
		if($contigname) {
			$len=length $seq;
			print outfile1 "$contigname\t$len\n"; 
			}
		$seq="";
		$contigname=$thisname;
		}
	else { $seq.=$_;}
	}
close infile1;
if($contigname) { $len=length $seq; print outfile1 "$contigname\t$len\n"; }
close outfile1;

print "Contigs stored in $contigsfna\n";
#system("rm $datapath/raw_fastq/par1.$format.gz; rm $datapath/raw_fastq/par2.$format.gz");
