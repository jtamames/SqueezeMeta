#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades) for several metagenomes that will be merged in the next step (merged mode).
#-- Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; } 
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$mappingfile,$tempdir,$megahit_soft,$assembler_options,$numthreads,$spades_soft,$canu_soft,$canumem,$prinseq_soft,$trimmomatic_soft,$mincontiglen,$resultpath,$contigsfna,$contigslen,$cleaning,$cleaningoptions);

#-- Read all the samples and store file names

my %ident;
my %samplefiles;

open(infile1,$mappingfile) || die "Cannot open samples file $mappingfile\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my($sample,$file,$iden,$flag)=split(/\t/,$_);
	next if($flag eq "noassembly");
	$ident{$file}=$iden;
	$samplefiles{$sample}{$file}=$iden;
	}
close infile1;

#-- Start working with samples, one at the time

	#-- Prepare files for the assembly

my($command,$trimmomatic_command);
foreach my $thissample(sort keys %samplefiles) {
	my($par1name,$par2name);
	print "Working for sample $thissample\n";
	system("rm $tempdir/par*fast*");
	my($seqformat,$gzformat,$numfiles,$cat1,$cat2);
	foreach my $thisfile(sort keys %{ $samplefiles{$thissample} }) {
		if($thisfile=~/gz$/) { $gzformat=".gz";  }	
		if($thisfile=~/fasta/) { $seqformat="fasta";  }
		elsif($thisfile=~/fastq/) { $seqformat="fastq";  }	
		if($ident{$thisfile} eq "pair1") { $par1name="$tempdir/par1.$seqformat$gzformat"; $cat1.="$datapath/raw_fastq/$thisfile "; $numfiles++; }
		elsif($ident{$thisfile} eq "pair2") { $par2name="$tempdir/par2.$seqformat$gzformat"; $cat2.="$datapath/raw_fastq/$thisfile "; }
		}
	print "Now merging files\n";
	$command="cat $cat1 > $par1name";
	print "$command\n";
	system $command;		
	if($cat2) {		#-- Support for single reads
		$command="cat $cat2 > $par2name";
		print "$command\n";
		system $command;		
		}
		

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


	#-- Run the assembly
	#-- For megahit

        my $assemblyname;
        my $ecode;
	if($assembler=~/megahit/i) { 
		system("rm -r $datapath/megahit"); 
		$assemblyname="$datapath/megahit/$thissample.final.contigs.fa";
		if(-e $par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name --k-list 29,39,59,79,99,119,141 -t $numthreads -o $datapath/megahit"; }
		else { $command="$megahit_soft $assembler_options -r $par1name --k-list 29,39,59,79,99,119,141 -t $numthreads -o $datapath/megahit"; }	#-- Support for single reads
		print "Running Megahit for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		system("mv $datapath/megahit/final.contigs.fa $assemblyname");
	}

	#-- For spades

	if($assembler=~/spades/i) { 
		system("rm -r $datapath/spades"); 
		$assemblyname="$datapath/spades/$thissample.contigs.fasta";
		if(-e $par2name) { $command="$spades_soft $assembler_options --meta --pe1-1 $par1name --pe1-2 $par2name -m 400 -t $numthreads -o $datapath/spades"; }
		else { $command="$spades_soft $assembler_options --meta --s1 $par1name -m 400 -t $numthreads -o $datapath/spades"; } #-- Support for single reads
		print "Running Spades for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		system("mv $datapath/spades/contigs.fasta $assemblyname");
	}
 
       #-- For canu

        if($assembler=~/canu/i) {
                system("rm -r $datapath/canu");
                $assemblyname="$datapath/canu/$thissample.contigs.fasta";
		$command="$canu_soft $assembler_options -p $project -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=$canumem oeaMemory=$canumem batMemory=$canumem mhapThreads=$numthreads mmapThreads=$numthreads ovlThreads=$numthreads ovbThreads=$numthreads ovsThreads=$numthreads corThreads=$numthreads oeaThreads=$numthreads redThreads=$numthreads batThreads=$numthreads gfaThreads=$numthreads merylThreads=$numthreads -nanopore-raw $par1name";
                print "Running canu for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
                system("mv $datapath/canu/$project.contigs.fasta $assemblyname");
        }




	#-- Run prinseq_lite for removing short contigs

	$contigsfna="$resultpath/01.$project.$thissample.fasta";	#-- Contig file from assembly
	$contigslen="$resultpath/01.$project.$thissample.lon";
	$command="$prinseq_soft -fasta $assemblyname -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna.prov";
	print "Running prinseq: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }

	#-- Now we need to rename the contigs for minimus2, otherwise there will be contigs with same names in different assemblies

	print "Renaming contigs\n"; 
	open(outfile1,">$contigsfna") || die;
	open(infile2,"$contigsfna.prov") || die;
	while(<infile2>) {
		chomp;
		if($_=~/^\>([^ ]+)/) { 
			my $tname=$1; 
			$_=~s/$tname/$tname\_$thissample/; 
		}
	print outfile1 "$_\n";
	}
	close infile2;
	close outfile1;
	system("rm $contigsfna.prov");

	#-- Run prinseq_lite for statistics

	$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $resultpath/01.$project.$thissample.stats";
        my $ecode = system $command;
        if($ecode!=0) { die "Error running command:    $command"; }
	

	#-- Counts length of the contigs (we will need it later)

	print "Counting lengths\n";
	my($seq,$thisname,$contigname);
	open(outfile2,">$contigslen") || die;
	print outfile2 "#-- Created by $0, ",scalar localtime,"\n";
	open(infile3,$contigsfna) || die;
	while(<infile3>) {
		chomp;
		next if !$_;
		if($_=~/^\>([^ ]+)/) {
			$thisname=$1;
			if($contigname) {
				my $len=length $seq;
				print outfile2 "$contigname\t$len\n"; 
			}
			$seq="";
			$contigname=$thisname;
		}
		else { $seq.=$_; }
	}
close infile3;
if($contigname) { my $len=length $seq; print outfile2 "$contigname\t$len\n"; }
close outfile2;

print "Contigs for sample $thissample stored in $contigsfna\n";

}                              #-- End of current sample

# system("rm $tempdir/par1.fastq.gz; rm $tempdir/par2.fastq.gz");
