#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage for genes by counting bases mapped to each feature

use strict;

my $bedfile=$ARGV[0];	# htseq|bedtools file
my $gff_file=$ARGV[1];	# gff file
my $sample=$ARGV[2];	# Sample ID

my(%long,%ids);

#-- Read gff file to get gen positions

open(infile1,$gff_file) || die;
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @k=split(/\t/,$_);
	my($id,$newid);
	if($_=~/ID\=([^;]+)/) { 
		$id=$1;
		my @listh=split(/\_/,$id);
		# ($base,$contigid,$num)=split(/\_/,$id);
		$newid="$k[0]\_$listh[$#listh]"; 
	}

	my $length=$k[4]-$k[3]+1;
	$long{$id}=$length;
	$ids{$id}=$newid;
}
close infile1;

#-- Read mapping file (bedcount file) to count bases to features 

my(%abun,%long);
open(infile2,$bedfile) || die "Cannot open htseq file $bedfile\n";
while(<infile2>) {
	chomp;
	next if !$_;
	my @fd=split(/\t/,$_);
	my $gn=$fd[3]; 
	my $gcount=$fd[5];                      #-- Bedtools (OJO! v<0.24)
	my $oid=$ids{$gn};
	# next if(!$oid);
	$abun{$oid}+=$gcount;
	$long{$oid}=$fd[4];
}
close infile2;

#-- Creating output file

print "# Created by $0 from $bedfile, ",scalar localtime,"\n";
print "# Gen	Bases	Coverage	Sample	Gen length\n";
foreach my $oid(sort keys %abun) {
	my $coverage=$abun{$oid}/$long{$oid};
	printf "$oid\t$abun{$oid}\t%.3f\t$sample\t$long{$oid}\n",$coverage;
	}
