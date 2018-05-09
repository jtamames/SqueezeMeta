#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Adds taxonomic assignments for bins, according to the consensus of the contigs belonging to it

$|=1;

use strict;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

	#-- Configuration variables from conf file

our($datapath,$alllog,$bintax,%bindirs);

	#-- Some configuration values for the algorithm
	
my $mincontigs=1;  		#-- Minimum number of genes for the contig 
my $minconsperc_asig=0.7;  	#-- Ratio contigs for the taxon/sum(genes all taxa). Therefore it only considers assigned contigs
my $minconsperc_total=0.5;  	#-- Ratio contigs for the taxon/number of contigs. Therefore it considers all (assigned+unassigned) contigs

my @ranks=('superkingdom','phylum','class','order','family','genus','species');
my(%tax,%taxlist);

	#-- Read taxonomic assignments for contigs

my $input=$alllog;
open(infile1,$input) || die "Cannot open $input\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my($contigid,$atax,$rest)=split(/\t/,$_);
	$tax{$contigid}=$atax;
	$atax=~s/\"//g;
	my @tf=split(/\;/,$atax);
	foreach my $uc(@tf) { 
		my($rank,$tax)=split(/\:/,$uc);
		if($rank ne "no rank") { $taxlist{$contigid}{$rank}=$tax; }
	}
}
close infile1;


open(outfile1,">$bintax") || die;
foreach my $binmethod(sort keys %bindirs) {		#-- For the current binning method
	my $bindir=$bindirs{$binmethod};
	print "Looking for $binmethod bins in $bindir\n";

	#-- Reading bin directories

	opendir(indir1,$bindir);
	my @files=grep(/fa$|fasta$/,readdir indir1);
	closedir indir1;
	
	#-- Working with each bin
	
	foreach my $k(@files) {
		my $tf="$bindir/$k";
		my $outf="$bindir/$k.tax";
		
		#-- We will create an output file for bin, with the contig names and consensus
		
		open(outfile2,">$outf") || die;			
		open(infile2,$tf) || die "Cannot open $tf\n";
		my $size=0;
		my %store=();
		
		#-- Read the contigs in the bin, with assignments and sizes
		
		while(<infile2>) {
			chomp;
			if($_=~/\>(.*)/) {
				my $contigid=$1;
				my $taxid=$tax{$contigid};
				print outfile2 "$_ $taxid\n";
				$store{$contigid}=1;
				}
			else { chomp; $size+=length $_; }
		}
		close infile2;
		
		#-- Call consensus() to fin the consensus taxon
		
		my ($constax,$rankc,$maxchim)=consensus(%store);
		
		#-- Write output
		
		printf outfile2 "Consensus: $constax\tTotal size: $size\tChimerism: %.3f\n",$maxchim;
		printf outfile1 "$binmethod\t$k\tConsensus: $constax\tTotal size: $size\tChimerism: %.3f\n",$maxchim;
		close outfile2;

 	}
}
close outfile1;


#----------------------- Calculate taxonomic consensus for variable ranks

sub consensus {

	my %store=shift;
	my($sep,$lasttax,$strg,$fulltax,$cattax)="";
	my $chimerism=0;
	
	#-- Loop for all ranks
	
	foreach my $rank(@ranks) { 
		my($totalcount,$totalas,$consf)=0;
		my %accumtax=();
		
		#-- For all contigs in the bin, count the taxa belonging to that rank
		
		foreach my $contig(keys %store) { 
			my $ttax=$taxlist{$contig}{$rank}; 
			$totalcount++;
			next if(!$ttax);
			$accumtax{$ttax}++;
			$totalas++; 
			# print "   $contig $ttax $accumtax{$ttax}\n";
			}
		next if(!$totalas);		#-- Exit if no contigs with assignment	
		
		#-- Consider the most abundant taxa, and look if it fulfills the criteria for consensus
							      
		my @listtax=sort { $accumtax{$b}<=>$accumtax{$a}; } keys %accumtax;
		my $mtax=$listtax[0];	
		my $times=$accumtax{$mtax};
		my($mtax2,$times2);
		if($listtax[1]) {	
			$mtax2=$listtax[1];	
			$times2=$accumtax{$mtax2};
			}
		else { $times2=0; }	
		my $percas=$times/$totalas;
		my $perctotal=$times/$totalcount;
		
		#-- If it does, store the consensus tax and rank
		
		if(($percas>=$minconsperc_asig) && ($perctotal>=$minconsperc_total) && ($totalcount>=$mincontigs) && ($times>$times2)) { 
			my $chimerism=($totalas-$times)/$totalas;
			 # print "***$mtax $times $percas $perctotal $totalas $chimerism\n";
			$cattax.="$mtax;";
			$fulltax.="$rank:$mtax;";
			$consf=1;
			$strg=$rank;
			# if($totalas!=$times) { print "***$contig $totalas $times\n"; }
			}
		last if(!$consf);
		}
	if(!$fulltax) { $fulltax="No consensus"; $strg="Unknown"; }
	
	#-- And return it
		
	return($fulltax,$strg,$chimerism);		
}
