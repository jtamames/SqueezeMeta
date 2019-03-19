#!/usr/bin/perl

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

our($resultpath,$tempdir,$cdhit_soft,$minimus2_soft,$toamos_soft,$prinseq_soft,$numthreads);

#-- Merges the assemblies in a single dataset

my $finalcontigs="$resultpath/01.$project.fasta";
if(-e $finalcontigs) { system("rm $finalcontigs"); }
my $merged="$tempdir/mergedassemblies.$project.fasta";
my $command="cat $resultpath/01*fasta > $merged";
system $command;
if(-z $merged) { die "$merged is empty\n"; }

#-- Uses cd-hit to identify and remove contigs contained in others

my $merged_clustered="$tempdir/mergedassemblies.$project.99.fasta";
$command="$cdhit_soft -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9";
print "Running cd-hit-est: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
if(-z $merged_clustered) { die "$merged_clustered is empty\n"; }

#-- Uses Amos to chage format to afg (for minimus2)

my $afg_format="$tempdir/mergedassemblies.$project.99.afg";
$command="$toamos_soft -s $merged_clustered -o $afg_format";
print "Transforming to afg format: $command\n";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
if(-z $afg_format) { die "$afg_format is empty\n"; }

#-- Uses minimus2 to assemble overlapping contigs

$command="$minimus2_soft $tempdir/mergedassemblies.$project.99 -D OVERLAP=100 -D MINID=95 -D THREADS=$numthreads";
print "Merging with minimus2: $command\n";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
if(-z $afg_format) { die "$afg_format is empty\n"; }

#-- Create the final result (overlapping contigs plus singletons)

system("cat $tempdir/mergedassemblies.$project.99.fasta $tempdir/mergedassemblies.$project.99.singletons.seq > $finalcontigs");

#-- Remove files from temp

system("rm -r $tempdir/mergedassemblies*");

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $finalcontigs -stats_len -stats_info -stats_assembly > $resultpath/01.$project.stats";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
	

#-- Count length of contigs (needed later)

my $contigslen="$resultpath/01.$project.lon";
print "Counting lengths\n";
open(outfile1,">$contigslen") || die;
open(infile1,$finalcontigs) || die;
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
