#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Merges individual assemblies using minimus2, for merged mode. It also uses cd-hit for excluding identical contigs

use strict;
use warnings;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

my $project=$ARGV[0];
$project=~s/\/$//;
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; } 
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($scriptdir, $resultpath,$interdir,$tempdir,$cdhit_soft,$extassembly,$minimus2_soft,$toamos_soft,$prinseq_soft,$numthreads);

open(out_tr,">$tempdir/merge.order");

#-- Merges the assemblies in a single dataset

my $finalcontigs="$resultpath/01.$project.fasta";
my($ecode,$command,$mergestep,$merged);
opendir(indir0,$interdir);
my @indassemblies=grep(/01.*fasta$/,readdir indir0);
my $numassem=$#indassemblies+1;
closedir indir0;

if($extassembly) { 
	print "External assembly provided: $extassembly. Overriding assembly\n";
	system("cp $extassembly $finalcontigs");
	}
else {
	if(-e "$tempdir/mergelog") { system("rm $tempdir/mergelog"); }
	while($numassem>1) {
		$mergestep++;	
		if(-e $finalcontigs) { system("rm $finalcontigs"); }
		$merged="$interdir/merged_$mergestep.$project.fasta";
		my $tomerge=system("$scriptdir/kmerdist.pl $project");
		open(infile0,"$tempdir/$project.2merge") || die;
		$_=<infile0>;
		chomp;
		my($sample1,$sample2,$mdist)=split(/\t/,$_);
		$sample1.=".fasta";
		$sample2.=".fasta";
		close infile0;
		print "MERGE $mergestep, $sample1 and $sample2\n";
		print out_tr "MERGE $mergestep, $sample1 and $sample2 ($mdist)\n";
		$command="cat $interdir/$sample1  $interdir/$sample2 > $merged";
		system $command;
		if(-z $merged) { die "$merged is empty\n"; }
		

			#-- Uses cd-hit to identify and remove contigs contained in others

			my $merged_clustered="$tempdir/mergedassemblies.$project.99.fasta";
			$command="$cdhit_soft -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9 > /dev/null 2>&1";
			print "Running cd-hit-est: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $merged_clustered) { die "$merged_clustered is empty\n"; }

			#-- Uses Amos to chage format to afg (for minimus2)

			my $afg_format="$tempdir/mergedassemblies.$project.afg";
			$command="$toamos_soft -s $merged_clustered -o $afg_format";
			print "Transforming to afg format: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Now separate reference and query
		
			parseafg($afg_format);
		
			#-- Uses minimus2 to assemble overlapping contigs

			$command="$minimus2_soft\_mod $tempdir/mergedassemblies.$project -D OVERLAP=100 -D MINID=95 -D THREADS=$numthreads > /dev/null 2>&1";
			print "Merging with minimus2: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Create the final result (overlapping contigs plus singletons)

			system("cat $tempdir/mergedassemblies.$project.fasta $tempdir/mergedassemblies.$project.singletons.seq > $merged.prov");

			open(outfile0,">$merged") || die;
			open(infile0,"$merged.prov") || die;
			while(<infile0>) {
				chomp;
				if($_=~/^\>(\d+)/) {
					my $newname="Merged\_$1\_$mergestep";
					print outfile0 ">$newname\n"; 
					}
				else { print outfile0 "$_\n"; }
				}
			close infile0;
			close outfile0;
			
			$numassem--;
			open(mlog,">>$tempdir/mergelog") || die;
			print mlog "$sample1\n$sample2\n";
			close mlog;
			#system("mv $interdir/$sample1 $interdir/$sample1.orig");
			#system("mv $interdir/$sample2 $interdir/$sample2.orig");
			#opendir(indir0,$interdir);
			#@indassemblies=grep(/fasta$/,readdir indir0);
			#$numassem=$#indassemblies+1;
			}
				

		#-- Remove files from temp

		 system("rm -r $tempdir/mergedassemblies*");
	}

system("mv $merged $finalcontigs");

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $finalcontigs -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
	

#-- Count length of contigs (needed later)

my $contigslen="$interdir/01.$project.lon";
print "Counting lengths\n";
open(outfile1,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile1,$finalcontigs) || die "Can't open $finalcontigs\n";
my($thisname,$contigname,$seq);
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
	$thisname=$1;
	if($contigname) {
		my $len=length $seq;
		print outfile1 "$contigname\t$len\n"; 
		}
	$seq="";
	$contigname=$thisname;
	}
	else { $seq.=$_; }
}
close infile1;
close outfile1;
close out_tr;

sub parseafg {
	my $inafg=shift;
	my($inred,$inpos,$accum);
	my(%order,%samples);
	open(inp1,$inafg) || die;
	while(<inp1>) {
	chomp;
	next if !$_;
	if($_=~/\{RED/) { $inred=1; }
	elsif($_=~/\}/) { $inred=0; }
	next if(!$inred);
	if($_=~/eid\:(.*)/) {
		$inpos++;
		my @m=split(/\_/,$1);
		my $ts=$m[$#m];
		$order{$inpos}=$ts;
		$samples{$ts}=1;
		}
	}
	close inp1;
	foreach my $p(sort keys %samples) { $accum++; $samples{$p}=$accum; }

	open(out1,">$inafg.1") || die;
	open(out2,">$inafg.2") || die;
	my($headsw,$intofrg,$numfrg,$intored,$numred,$samp,$tofile);
	open(inp2,$inafg) || die;
	while(<inp2>) {
		chomp;
		next if !$_;
		if($_=~/\{FRG/) { $headsw=1; $intofrg=1; $numfrg++; }
		if($_=~/\{RED/) { $headsw=1; $intofrg=0; $intored=1; $numred++; }
		if(!$headsw) {
			print out1 "$_\n";
			print out2 "$_\n";
			next;
			}
		if($intofrg) {	
			$samp=$order{$numfrg};
			$tofile=$samples{$samp}; 
			if($tofile==1) { print out1 "$_\n"; } else { print out2 "$_\n"; } 
			}
		elsif($intored) {
			$samp=$order{$numred};
			$tofile=$samples{$samp}; 
			if($tofile==1) { print out1 "$_\n"; } else { print out2 "$_\n"; } 
			}
		
		}
	close inp2;
	}
		


