#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Counts the sizes for all taxa in the contiglog file created by summarycontigs3

use strict;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/SqueezeMeta_conf.pl";

	#-- Configuration variables from conf file

our($datapath,$resultpath,$contigslen,$alllog,$taxlist,$contigcov,$mcountfile);

my(%lon,%taxa,%abund,%abundreads,%samples,%accum,%accumbases,%accumreads,%taxcorr);

	#-- Read contig lengths

open(infile1,$contigslen);
while(<infile1>) {  
	chomp;
	next if !$_;
	my ($node,$len)=split(/\t/,$_);
	$lon{$node}=$len;
	}
close infile1;

	#-- Read contiglog file to get taxonomic assignment for contigs

open(infile2,$alllog) || die;
while(<infile2>) {
	chomp;
	next if !$_;
	my($node,$tax,$rest)=split(/\t/,$_);
	if($tax eq "No consensus") { $tax="Unknown"; }
	if(!$tax) { $tax="Unknown"; }
	$taxa{$node}=$tax;
	}
close infile2;

	#-- Read contigcov file to get abundances of each contig

open(infile3,$contigcov) || die;
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @f=split(/\t/,$_);
	my $sample=$f[$#f];
	$abund{$f[0]}{$sample}=$f[5];
	$abundreads{$f[0]}{$sample}=$f[4];
	$samples{$sample}=1; 
	my $node=$f[0];
	my $tlong=$lon{$node};
	my $tax=$taxa{$node};
	if(!$tax) { $tax="Unknown"; }
	my @tx=split(/\;/,$tax);
	my $string="";
	
	#-- For all the ranks of the current contig, add the contig size to the corresponding taxon
	
	for(my $n=0; $n<=$#tx; $n++) {
		$string.="$tx[$n];";
		$accum{$string}+=$tlong;
		
		#-- Add also bases and reads
		
		foreach my $samp(keys %samples) { 
			$accumbases{$string}{$samp}+=$abund{$node}{$samp};
			$accumreads{$string}{$samp}+=$abundreads{$node}{$samp}; 
 			}
		}
 
	#print "$f[0] $sample $f[5]\n";
	}
close infile3;

	#-- Read the equivalence between taxa and rank

open(infile4,$taxlist) || die;
my $nname;
while(<infile4>) {
	chomp;
	next if !$_;
	my($id,$ttax,$trank)=split(/\t/,$_);
	$taxcorr{$ttax}=$trank;
	if($trank eq "species") {
		my @wd=split(/\s+/,$ttax);
		if($wd[0] eq "Candidatus") { $nname="$wd[0] $wd[1] $wd[2]"; } else { $nname="$wd[0] $wd[1]"; }
		$taxcorr{$nname}=$trank;
		}
	}
close infile4;
 
	#-- Write the output file 
 
open(outfile1,">$mcountfile") || die "Cannot open $mcountfile\n";
print outfile1 "Rank\tTaxon\tAccumulated contig size";
foreach my $samp(sort keys %samples) { print outfile1 "\t$samp reads\t$samp bases"; }
print outfile1 "\n";
foreach my $kk(sort { $accum{$b}<=>$accum{$a}; } keys %accum) { 
	my $k=$kk;
	$k=~s/\;$//;
	my @l=split(/\;/,$k);
	my($rank,$tn)=split(/\:/,$l[$#l]);
	if(!$rank) { $rank="Unknown"; }				  
	print outfile1 "$rank\t$k\t$accum{$kk}"; 
	foreach my $samp(sort keys %samples) { print outfile1 "\t$accumreads{$kk}{$samp}\t$accumbases{$kk}{$samp}"; }
	print outfile1 "\n";
	}
close outfile1;
