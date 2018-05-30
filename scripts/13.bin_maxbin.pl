#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs MaxBin for binning

use strict;
use warnings;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($contigsfna,%bindirs,$contigcov,$maxbin_soft,$numthreads);

my %allcontigs;

	#-- Reading contigs

open(infile1,$contigsfna) || die;
while(<infile1>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { $allcontigs{$1}=1; }
	}
close infile1;

	#-- Creating binning directory

my $dirbin=$bindirs{maxbin};
if(-d $dirbin) {} else { system "mkdir $dirbin"; }

	#-- Creating abundance file

my $abundlist="$dirbin/abund.list";
open(outfile1,">$abundlist") || die;	#-- Stores the list of files for the abundance of contigs (one per sample)
open(infile2,$contigcov) || die "Cannot find contig coverage file $contigcov\n";
my $currsample;
my %tcontigs;
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	if(!$currsample || ($k[$#k] ne $currsample)) {	#-- If we finished reading one sample, write the data
		foreach my $k(sort keys %tcontigs) { print outfile2 "$k\t0\n"; }
		close outfile2;
		$currsample=$k[$#k];
		print outfile1 "$dirbin/$currsample.abund\n";
		open(outfile2,">$dirbin/$currsample.abund") || die;	#-- Stores the abundances for current sample
		%tcontigs=%allcontigs;
  
		}
	print outfile2 "$k[0]\t$k[1]\n";
	delete $tcontigs{$k[0]};
	}
close infile2;	   
foreach my $k(sort keys %tcontigs) { print outfile2 "$k\t0\n"; }
close outfile2;
close outfile1;

my $command="perl $maxbin_soft -thread $numthreads -contig $contigsfna -abund_list $abundlist -out $dirbin/maxbin";
print "Now running Maxbin: $command\n";
system $command;
