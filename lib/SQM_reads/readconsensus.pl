#!/usr/bin/env perl

use strict;
use Tie::IxHash;
use Cwd;
use lib ".";

$|=1;

my $projectdir=$ARGV[0];	#-- directory for the file
my $input=$ARGV[1];		#-- wranks file
my $filter=$ARGV[2];		#-- With/without idfilters

my @ranks=('superkingdom','phylum','class','order','family','genus','species');
my @ranksabb=('k','p','c','o','f','g','s');

my $mingenes9=1;			     #-- STEP9: Minimum number of genes for the contig
my $minconsperc_asig9=0.7;      #-- STEP9: Ratio genes for the taxon/sum(genes all taxa). Therefore it only considers assigned genes
my $minconsperc_total9=0.1;     #-- STEP9: Ratio genes for the taxon/number of genes. Therefore it considers all (assigned+unassigned) genes

my @k=split(/\//,$projectdir);
my $installpath = "/media/disk5/tamames/SqueezeMeta";
my $extdatapath  = "$installpath/data";
my $databasepath = "/media/disk7/fer/SqueezeMeta/db";


#-- Reading taxonomic infomation (extracted from NCBI's taxonomy)

my(%taxcorr,%numorfs,%allcontigs,%orfs);

my $taxlist   = "$extdatapath/alltaxlist.txt";
open(infile1,$taxlist) || die "Can't open $taxlist\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @y=split(/\t/,$_);
	$taxcorr{$y[1]}=$y[2];
	}
close infile1;

my %parents;
open(infile0,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile0>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	$tax=~s/\[|\]//g;
	$parents{$tax}{wranks}=$par;
	my @m=split(/\;/,$par);
	foreach my $y(@m) {
		my($rt,$gtax)=split(/\:/,$y);
		$parents{$tax}{noranks}.="$gtax;"; 
		}
	chop $parents{$tax}{noranks};
	}
close infile0;


#-- Reading taxonomic information for the genes (wrank file)

my %taxlist;
print "  Reading $input\n";
open(infile3,$input) || die "Can't open $input\n";
while(<infile3>) {		#-- Looping on the ORFs
	my $contigid;
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my($node,$atax)=split(/\t/,$_);		#-- $node contains ORF name
	$atax=~s/\"//g;
	my @tf=split(/\;/,$atax);		#-- This contains rank:taxon pairs for the ORF
	my @id=split(/\_/,$node);
	$contigid=$id[0];	
	$numorfs{$contigid}++;

	#-- Stores rank and taxa for all the ORFs in the contig

	foreach my $uc(@tf) { 
		my ($rank,$tax)=split(/\_/,$uc);
		if($rank ne "n") { $taxlist{$contigid}{$rank}{$node}=$tax;  }
		}
	}
close infile3;

my($outputlong,$outputshort);
if($filter eq "idfilter") {
	$outputlong="$projectdir/readconsensus.log";
	$outputshort="$projectdir/readconsensus.txt";
	}
else {
	$outputlong="$projectdir/readconsensus\_noidfilter.log";
	$outputshort="$projectdir/readconsensus\_noidfilter.txt";
	}

print "  Output in $outputshort\n";
open(outfile1,">$outputlong") || die "Can't open $outputlong for writing\n";
open(outfile2,">$outputshort") || die "Can't open $outputshort for writing\n";

foreach my $contig(keys %numorfs) {
	my ($sep,$lasttax,$strg,$cattax,$fulltax,$lasttax)="";
	my ($consensus,$schim,$chimerism)=0;
	my(%consensus,%chimeracheck)=();
	my $totalcount=$numorfs{$contig}; 
	next if(!$totalcount);	
	print outfile1 "$contig\t$totalcount ORF";
	if($totalcount>1) { print outfile1 "s"; }
	print outfile1 "\n";
	
	#-- We loop on ranks
	
	foreach my $rank(@ranksabb) { 
		my $totalas=0;
		my %accumtax=();
		
		#-- And count how many times each taxon appears
		
		foreach my $orf(sort keys %{ $taxlist{$contig}{$rank} }) {  
			my $ttax=$taxlist{$contig}{$rank}{$orf}; 
			$accumtax{$ttax}++;
			if($ttax) { $totalas++; }
			}
		next if(!$totalas);	#-- Don't proceed if no ORF has been assigned at this current rank


		#-- We calculate if the most abundant taxon fulfills the criteria for assignment
								      
		my @listtax=sort { $accumtax{$b}<=>$accumtax{$a}; } keys %accumtax;
		my $mtax=$listtax[0];	
		my $times=$accumtax{$mtax};
		my($mtax2,$times2);
		if($listtax[1]) {	
			$mtax2=$listtax[1];	
			$times2=$accumtax{$mtax2};	
			}
		if(!$totalcount) { print "  Warning: Contig $contig has zero genes\n"; }
		my $percas=$accumtax{$mtax}/$totalas;
		my $perctotal=$accumtax{$mtax}/$totalcount;
		if(!$times2) { $times2=0; }
		#-- And if it does, we also calculate the disparity index of the contig for this rank

		if(($percas>=$minconsperc_asig9) && ($perctotal>=$minconsperc_total9) && ($totalcount>=$mingenes9) && ($times>$times2)) { 
	
			#-- Calculation of disparity for this rank
			my($chimera,$nonchimera,$unknown)=0;
			foreach my $orf(sort keys %{ $orfs{$contig} }) { 
				my $ttax=$taxlist{$contig}{$rank}{$orf}; 
				foreach my $orf2(sort keys %{ $orfs{$contig} }) { 
					my $ttax2=$taxlist{$contig}{$rank}{$orf2}; 
					next if($orf ge $orf2);
					if($chimeracheck{$orf}{$orf2} eq "chimera") {	#-- If it was a chimera in previous ranks, it is a chimera now
						$chimera++;
						next;
					}	
					if($chimeracheck{$orf}{$orf2} eq "unknown") {	#-- If it was an unknown in previous ranks, it is a unknown now
						$unknown++;
						next;
					}						
					if(($ttax && (!$ttax2)) || ($ttax2 && (!$ttax)) || ((!$ttax2) && (!$ttax))) { $chimeracheck{$orf}{$orf2}="unknown"; } #-- Unknown when one of the ORFs has no classification at this rank
					elsif($ttax eq $ttax2) { $chimeracheck{$orf}{$orf2}="nochimera"; $nonchimera++; }
					else { $chimeracheck{$orf}{$orf2}="chimera"; $chimera++; }
					}
			}


			# $chimerism=($totalas-$times)/$totalas;
			my $totch=$chimera+$nonchimera;
			if($totch) { $chimerism=$chimera/($chimera+$nonchimera); } else { $chimerism=0; }
			$consensus{$rank}=$mtax; 
			printf outfile1 "   $rank: $totalas\t$mtax: $times\tDisparity: %.3f\n",$chimerism; 

			#-- Then we add the new taxon to the taxa found for previous ranks

			$cattax.="$mtax;";
			$fulltax.="$rank\_$mtax;";
			$lasttax=$mtax;
			$strg=$rank;
			$schim=$chimerism;
			$consensus=1;
			}	
								      
		# if($totalas!=$times) { print "***$contig $totalas $times\n"; }
		last if(!$consensus{$rank});	#-- And exit if there is no consensus for this rank
		}

	#-- Finally, write the output
		
	if(!$consensus) { $cattax="No consensus"; $strg="Unknown"; }		
	my $abb=$parents{$lasttax}{wranks};	
	$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
	printf outfile2 "$contig\t$abb\t$strg\tDisparity: %.3f\tORFs: $numorfs{$contig}\n",$chimerism;
                                 }
close outfile1;
close outfile2;

