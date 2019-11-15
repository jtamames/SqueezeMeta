#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs binning with Metabat2

use strict;
use Cwd;
use lib ".";

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

#-- Configuration variables from conf file

our($contigsfna,$contigcov,$metabat_soft,$alllog,$tempdir,$mappingfile,$maxchimerism15,$mingenes15,$smallnoannot15,%bindirs);
my %skip;

print "Reading samples from $mappingfile\n";   #-- We will exclude samples with the "noassembly" flag
open(infile0,$mappingfile) || die "Can't open $alllog\n";
while(<infile0>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	if($_=~/nobinning/) { $skip{$t[0]}=1; }
	}
close infile0;

	#-- Reading contigs

my @allcontigs;
my(%abun,%allsets,%contiglen,%sumaver,%allcontigs);

open(infile1,$alllog) || die "Can't open $alllog\n";
while(<infile1>) { 
	chomp;
	next if !$_;
	my @r=split(/\t/,$_);
	my($chimlevel,$numgenes);
	if($r[3]=~/Disparity\: (.*)/) { $chimlevel=$1; }
	if($r[4]=~/Genes\: (.*)/) { $numgenes=$1; } 
	if(!$numgenes) { $numgenes=0; } 
	if(($numgenes>=$mingenes15) && ($chimlevel<=$maxchimerism15)) { push(@allcontigs,$r[0]); $allcontigs{$r[0]}=1; }	
	if($smallnoannot15 && ($numgenes<=1) && ($r[1] eq "Unknown")) { delete $allcontigs{$r[0]}; }
	}
close infile1;

my $tempfasta="$tempdir/bincontigs.fasta";
open(outfile1,">$tempfasta") || die "Can't open $tempfasta for writing\n";
open(infile1,$contigsfna) || die "Can't open $contigsfna\n";
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

my $dirbin=$bindirs{metabat2};
if(-d $dirbin) {} else { system "mkdir $dirbin"; }

	#-- Reading contig abundances

open(infile2,$contigcov) || die "Can't find contig coverage file $contigcov\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	next if($skip{$k[$#k]});
	$abun{$k[0]}{$k[$#k]}=$k[1];
	$allsets{$k[$#k]}++;
	$contiglen{$k[0]}=$k[3];
	$sumaver{$k[0]}+=$k[1];
	}
close infile2;

	#-- Creating abundance file

my $depthfile="$dirbin/contigs.depth.txt";
open(outfile1,">$depthfile") || die "Can't open $depthfile for writing\n";
print outfile1 "contigName\tcontigLen\ttotalAvgDepth";
foreach my $dataset(sort keys %allsets) { print outfile1 "\t$dataset.bam\t$dataset.bam-var"; }
print outfile1 "\n";
foreach my $contig(@allcontigs) {
	printf outfile1 "$contig\t$contiglen{$contig}\t%.4f",$sumaver{$contig};
	foreach my $dataset(sort keys %allsets) { 
		my $dat=$abun{$contig}{$dataset} || "0";
		printf outfile1 "\t%.4f\t0",$dat;
		}
print outfile1 "\n";
				   }

close outfile1;

	#-- Running metabat

my $command="$metabat_soft -t 8 -i $tempfasta -a $depthfile -o $dirbin/metabat2 --saveTNF saved_1500.tnf --saveDistance saved_1500.dist";
print "Now running metabat2: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

