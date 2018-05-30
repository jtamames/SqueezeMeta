#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs binning with Metabat2

use strict;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($contigsfna,$contigcov,$metabat_soft,%bindirs);

	#-- Reading contigs

my @allcontigs;
my(%abun,%allsets,%contiglen,%sumaver);

open(infile1,$contigsfna) || die;
while(<infile1>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { push(@allcontigs,$1); }
            }
close infile1;

	#-- Creating binning directory

my $dirbin=$bindirs{metabat2};
if(-d $dirbin) {} else { system "mkdir $dirbin"; }

	#-- Reading contig abundances

open(infile2,$contigcov) || die "Cannot find contig coverage file $contigcov\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	$abun{$k[0]}{$k[$#k]}=$k[1];
	$allsets{$k[$#k]}++;
	$contiglen{$k[0]}=$k[3];
	$sumaver{$k[0]}+=$k[1];
	}
close infile2;

	#-- Creating abundance file

my $depthfile="$dirbin/contigs.depth.txt";
open(outfile1,">$depthfile") || die;
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

my $command="$metabat_soft -t 8 -i $contigsfna -a $depthfile -o $dirbin/metabat2 --saveTNF saved_1500.tnf --saveDistance saved_1500.dist";
print "Now running metabat2: $command\n";
system $command;
