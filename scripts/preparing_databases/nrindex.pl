#!/usr/bin/env perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- From the nr database, extracts the taxonomic origin of the entries
#-- Needed to create LCA database. Creates nr.taxlist.db, containing the species corresponding to each nr entry

use strict;

#-- Output file. MUST POINT TO DATABASE DIRECTORY

my $databasedir=$ARGV[0];
my $nrfasta="$databasedir/nr.faa";	#-- Previously downloaded

my $resdir="$databasedir/LCA_tax";
if(-d $resdir) {} else { system("mkdir $resdir"); }
my $outfile="$resdir/nr.taxlist.db";
open(outfile1,">$outfile") || die;
open(infile1,$nrfasta) || die "Cannot open $nrfasta\n";
my %accum;
my $specname;
while(<infile1>) {
	chomp;
	next if(!$_ || ($_!~/^\>/));
	my ($id,$rest)=split(/\s+/,$_);
	my @f=split(/\|/,$id);
	my @spec=($_=~/\[.*?\]/g);
	foreach my $k(@spec) {
		$k=~s/\[|\]//g;
		my @wd=split(/\s+/,$k);
		next if(!$wd[1]);
		if($k=~/virus|phage|symbiont/i) { 
			$specname=$k;
			$specname=~s/\(.*//g;
			$specname=~s/ strain //g;
			$specname=~s/ genotype //g;
			$specname=~s/\s+$//;
			$specname=~s/\s+/ /g;
			 }
		elsif($wd[0] eq "Candidatus") { $specname="$wd[0] $wd[1] $wd[2]"; } 
		else  { $specname="$wd[0] $wd[1]"; }  
		$accum{$specname}++;
		}
	# print outfile1 "$f[1]\t$f[3]\t";	#-- Old nr, with GIs
	$id=~s/^\>//;				#-- New nr, unique identifier
	print outfile1 "$id\t";
	my $string="";		      
	foreach my $pr(sort keys %accum) { $string.="$pr;"; }
	chop $string;
	print outfile1 "$string\n";
	%accum=();
	}
	
close outfile1;
close infile1;
print "File created: $outfile\n";	    
