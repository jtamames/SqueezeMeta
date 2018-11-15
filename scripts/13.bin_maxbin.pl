#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs MaxBin for binning

use strict;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($databasepath,$contigsfna,%bindirs,$contigcov,$maxbin_soft,$alllog,$tempdir,$numthreads);

my $maxchimerism=0.1;	#-- Threshold for excluding chimeric contigs
my $mingenes=1;		#-- Threshold for excluding small contigs (few genes than this)
my $smallnoannot=1;	#-- For excluding contigs with just one gene an no annotation

my %allcontigs;

	#-- Reading contigs

#open(infile1,$contigsfna) || die;
#while(<infile1>) {
#	chomp;
#	if($_=~/^\>([^ ]+)/) { $allcontigs{$1}=1; }
#	}
#close infile1;

print "Reading from $alllog\n";
open(infile1,$alllog) || die;
while(<infile1>) { 
	chomp;
	next if !$_;
	my @r=split(/\t/,$_);
	my($chimlevel,$numgenes);
	if($r[3]=~/Chimerism level\: (.*)/) { $chimlevel=$1; }
	if($r[4]=~/Genes\: (.*)/) { $numgenes=$1; }
	if(!$numgenes) { $numgenes=0; } 
	if(($numgenes>=$mingenes) && ($chimlevel<=$maxchimerism)) { $allcontigs{$r[0]}=1;  }	
	if($smallnoannot && ($numgenes<=1) && ($r[1] eq "Unknown")) { delete $allcontigs{$r[0]}; }
	}
close infile1;

my $tempfasta="$tempdir/bincontigs.fasta";
open(outfile1,">$tempfasta") || die;
open(infile1,$contigsfna) || die;
my $ingood=0;
while(<infile1>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { 
		my $tc=$1;
		if($allcontigs{$tc}) { $ingood=1; } else { $ingood=0; } 
		}
	if($ingood) { print outfile1 "$_\n"; }
	}
close infile1;
close outfile1; 

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
	if($allcontigs{$k[0]}) { print outfile2 "$k[0]\t$k[1]\n"; }
	delete $tcontigs{$k[0]};
	}
close infile2;	   
foreach my $k(sort keys %tcontigs) { print outfile2 "$k\t0\n"; }
close outfile2;
close outfile1;

my $command="perl $maxbin_soft -thread $numthreads -contig $tempfasta -abund_list $abundlist -out $dirbin/maxbin -markerpath $databasepath/marker.hmm";
print "Now running Maxbin: $command\n";
system $command;
