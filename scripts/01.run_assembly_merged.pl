#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades) for several metagenomes that will be merged in the next step (merged mode).
#-- Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 
use Term::ANSIColor qw(:constants);

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";


#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$mappingfile,$extassembly,$tempdir,$interdir,$megahit_soft,$assembler_options,$numthreads,$spades_soft,$canu_soft,$canumem,$prinseq_soft,$trimmomatic_soft,$mincontiglen,$resultpath,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$methodsfile,$syslogfile);

#-- Read all the samples and store file names

exit if $extassembly;

my %ident;
my %samplefiles;

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

open(infile1,$mappingfile) || die "Can't open samples file $mappingfile\n";
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
	print BOLD "\n  Sample $thissample\n"; print RESET;
	system("rm $tempdir/par*fast* > /dev/null 2>&1");
	my($seqformat,$gzformat,$numfiles,$cat1,$cat2);
	foreach my $thisfile(sort keys %{ $samplefiles{$thissample} }) {
		if($thisfile=~/gz$/) { $gzformat=".gz";  }	
		if($thisfile=~/fasta/) { $seqformat="fasta";  }
		elsif($thisfile=~/fastq/) { $seqformat="fastq";  }	
		if($ident{$thisfile} eq "pair1") { $par1name="$tempdir/par1.$seqformat$gzformat"; $cat1.="$datapath/raw_fastq/$thisfile "; $numfiles++; }
		elsif($ident{$thisfile} eq "pair2") { $par2name="$tempdir/par2.$seqformat$gzformat"; $cat2.="$datapath/raw_fastq/$thisfile "; }
		}
	# print "Now merging files\n";
	$command="cat $cat1 > $par1name";
	# print "$command\n";
	print outsyslog "Merging read files: $command\n";
	system $command;		
	if($cat2) {		#-- Support for single reads
		$command="cat $cat2 > $par2name";
		# print "$command\n";
		print outsyslog "Merging read files: $command\n";
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
		if(-e $orig2) { $trimmomatic_command="$trimmomatic_soft PE -threads $numthreads -phred33 $orig1 $orig2 $par1name $par1name.removed $par2name $par2name.removed $cleaningoptions > /dev/null 2>&1"; }
		else { $trimmomatic_command="$trimmomatic_soft SE -threads $numthreads -phred33 $orig1 $par1name $cleaningoptions > /dev/null 2>&1"; }

		if($cleaning) {
			print "  Running Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20) for filtering reads\n";
			print outsyslog "Running Trimmomatic: $trimmomatic_command\n";
			my $ecode = system $trimmomatic_command;
			if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
			print outmet "Quality filtering was done using Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20)\n";
			}
		}


	#-- Run the assembly
	#-- For megahit

        my $assemblyname;
        my $ecode;
	if($assembler=~/megahit/i) { 
		system("rm -r $datapath/megahit > /dev/null 2>&1"); 
		$assemblyname="$datapath/megahit/$thissample.final.contigs.fa";
		if(-e $par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name --k-list 29,39,59,79,99,119,141 -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }
		else { $command="$megahit_soft $assembler_options -r $par1name --k-list 29,39,59,79,99,119,141 -t $numthreads -o $datapath/megahit > /dev/null 2>&1"; }	#-- Support for single reads
		print "  Running Megahit (Li et al 2015, Bioinformatics 31(10):1674-6) for $thissample\n";
		print outsyslog "Running Megahit for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		system("mv $datapath/megahit/final.contigs.fa $assemblyname");
	}

	#-- For spades

	if($assembler=~/spades/i) { 
		system("rm -r $datapath/spades > /dev/null 2>&1"); 
		$assemblyname="$datapath/spades/$thissample.contigs.fasta";
		if(-e $par2name) { $command="$spades_soft --meta --pe1-1 $par1name --pe1-2 $par2name -m 400 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile  2>&1"; }
		else { $command="$spades_soft --meta --s1 $par1name -m 400 $assembler_options -t $numthreads -o $datapath/spades > /dev/null 2>&1"; } #-- Support for single reads
		print "  Running Spades (Li et al 2015, Bioinformatics 31(10):1674-6) for $thissample\n";
		print outsyslog "Running Spades for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		system("mv $datapath/spades/contigs.fasta $assemblyname");
	}
 
       #-- For canu

        if($assembler=~/canu/i) {
                system("rm -r $datapath/canu > /dev/null 2>&1");
                $assemblyname="$datapath/canu/$thissample.contigs.fasta";
		if($canumem eq "NF") {
			print "  Setting available memory for Canu\n";
			my %mem=get_mem_info;
			my $ram=$mem{"MemAvailable"};
			$canumem=sprintf('%.1f',$ram/1000000);
			$canumem*=0.8;
			print "AVAILABLE (free) RAM memory: $ram\nWe will set canu to $canumem. You can override this setting using the -canumem option\n";
			print outsyslog "canumem set to $canumem (Free Mem $ram bytes)\n";
			}
 		$command="$canu_soft -p $project -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=$canumem oeaMemory=$canumem batMemory=$canumem mhapThreads=$numthreads mmapThreads=$numthreads ovlThreads=$numthreads ovbThreads=$numthreads ovsThreads=$numthreads corThreads=$numthreads oeaThreads=$numthreads redThreads=$numthreads batThreads=$numthreads gfaThreads=$numthreads merylThreads=$numthreads -nanopore-raw $assembler_options $par1name >> $syslogfile 2>&1";
                print "  Running canu (Koren et al 2017, Genome Res 27(5):722-36) for $thissample\n";
		print outsyslog "Running canu for $thissample: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
                system("mv $datapath/canu/$project.contigs.fasta $assemblyname");
        }




	#-- Run prinseq_lite for removing short contigs

	$contigsfna="$interdir/01.$project.$thissample.fasta";	#-- Contig file from assembly
	$contigslen="$interdir/01.$project.$thissample.lon";
	$command="$prinseq_soft -fasta $assemblyname -min_len $mincontiglen -out_good $tempdir/prinseq; mv $tempdir/prinseq.fasta $contigsfna.prov";
	print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen\n";
	print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	if($mincontiglen>200) { print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n"; }

	#-- Now we need to rename the contigs for minimus2, otherwise there will be contigs with same names in different assemblies

	print "  Renaming contigs\n"; 
	open(outfile1,">$contigsfna") || die "Can't open $contigsfna for writing\n";
	open(infile2,"$contigsfna.prov") || die "Can't open $contigsfna.prov\n";
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

	$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$project.$thissample.stats";
	print outsyslog "Run prinseq for statistics: $command\n";
        my $ecode = system $command;
        if($ecode!=0) { die "Error running command:    $command"; }
	

	#-- Counts length of the contigs (we will need it later)

	print "  Counting contig lengths\n";
	my($seq,$thisname,$contigname);
	open(outfile2,">$contigslen") || die "Can't open $contigslen for writing\n";
	print outfile2 "#-- Created by $0, ",scalar localtime,"\n";
	open(infile3,$contigsfna) || die "Can't open $contigsfna\n";
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

print "  Contigs for sample $thissample stored in $contigsfna\n";

}                              #-- End of current sample

# system("rm $tempdir/par1.fastq.gz; rm $tempdir/par2.fastq.gz");
if($assembler=~/megahit/i) { print outmet "Assembly was done using Megahit (Li et al 2015, Bioinformatics 31(10):1674-6)\n"; }
elsif($assembler=~/spades/i) { print outmet "Assembly was done using SPAdes (Bankevich et al 2012, J Comp Biol 19(5):455-77)\n"; }
elsif($assembler=~/canu/i) { print outmet "Assembly was done using Canu (Koren et al 2017, Genome Res 27(5):722-36)\n"; }
if($mincontiglen>200) { print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
close outmet; 
close outsyslog;
