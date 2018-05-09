#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Assigns a consensus annotation for the whole contigs using the annotations for the genes it contains
#-- Also calculates chimerism index for each contig

use strict;
use Tie::IxHash;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$fun3tax,$taxlist,$aafile);

#-- Some local configuration variables

my @ranks=('superkingdom','phylum','class','order','family','genus','species');
my $mingenes=1;			#-- Minimum number of genes for the contig 
my $minconsperc_asig=0.7;	#-- Ratio genes for the taxon/sum(genes all taxa). Therefore it only considers assigned genes
my $minconsperc_total=0.5;	#-- Ratio genes for the taxon/number of genes. Therefore it considers all (assigned+unassigned) genes

#-- Input and output files

my $input="$fun3tax.wranks";
my $outputlong="$resultpath/08.$project.fulllog";
my $outputshort="$resultpath/08.$project.contiglog";

#-- Reading taxonomic infomation (extracted from NCBI's taxonomy)

my(%taxcorr,%numorfs,%allcontigs);

open(infile1,$taxlist) || die;
while(<infile1>) {
	chomp;
	next if !$_;
	my @y=split(/\t/,$_);
	$taxcorr{$y[1]}=$y[2];
	}
close infile1;

#-- Reading genes in contigs

open(infile2,$aafile) || die;
while(<infile2>) {
	chomp;
	next if($_!~/^\>/);
	if($_=~/^\>([^ ]+)/) { 
		my $contig=$1;
		$contig=~s/\_\d+$//;
		$numorfs{$contig}++; 
		$allcontigs{$contig}=1;
		}
	}
close infile2;

#-- Reading taxonomic information for the genes (wrank file)

my %taxlist;
print "Reading $input\n";
open(infile3,$input) || die "Cannot open $input\n";
while(<infile3>) {		#-- Looping on the ORFs
	my $contigid;
	chomp;
	next if !$_;
	my($node,$atax)=split(/\t/,$_);		#-- $node contains ORF name
	$atax=~s/\"//g;
	my @tf=split(/\;/,$atax);		#-- This contains rank:taxon pairs for the ORF
	my @id=split(/\_/,$node);
	# @id=split(/\|/,$node);
	
	#-- Different nomenclature for contigs in spades and other assemblers
	
	if($id[0]=~/NODE/) { $contigid="$id[0]\_$id[1]"; } else { $contigid="$id[0]"; }	

	$contigid=$node;
	$contigid=~s/\_\d+$//;
	$contigid=~s/\|\d+$//; 		#-- This is the contig name the current ORF belongs to

	#-- Stores rank and taxa for all the ORFs in the contig

	foreach my $uc(@tf) { 
		my ($rank,$tax)=split(/\:/,$uc);
		if($rank ne "no rank") { $taxlist{$contigid}{$rank}{$node}=$tax; }
		}
	}
close infile3;

#-- Preparing output files

print "Writing output to $outputshort\n";
open(outfile1,">$outputlong") || die "Cannot open $outputlong\n";
open(outfile2,">$outputshort") || die "Cannot open $outputshort\n";
print outfile1 "#- Created by $0 with data from $input, mingenes=$mingenes, minconsperc_asig=$minconsperc_asig, minconsperc_total=$minconsperc_total, ",scalar localtime,"\n";
print outfile2 "#- Created by $0 with data from $input, mingenes=$mingenes, minconsperc_asig=$minconsperc_asig, minconsperc_total=$minconsperc_total, ",scalar localtime,"\n";

foreach my $contig(keys %allcontigs) {
	my ($sep,$lasttax,$strg,$cattax,$fulltax)="";
	my ($consensus,$schim,$chimerism)=0;
	my %consensus=();
	my $totalcount=$numorfs{$contig}; 	
	print outfile1 "$contig\t$totalcount genes\n";
	
	#-- We loop on ranks
	
	foreach my $rank(@ranks) {
		my $totalas=0;
		my %accumtax=();
		
		#-- And count how many times each taxon appears
		
		foreach my $orf(keys %{ $taxlist{$contig}{$rank} }) { 
			my @y=split(/\_/,$orf);
			my $tpos=$y[$#y];
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
		my $percas=$accumtax{$mtax}/$totalas;
		my $perctotal=$accumtax{$mtax}/$totalcount;
		if(!$times2) { $times2=0; }
		# if($contig eq "k119_42524") {  printf "  $mtax $rank $accumtax{$mtax} $totalas $totalcount %.2f %.2f\n",$percas,$perctotal; }	      

		#-- And if it does, we also calculate the chimerism index of the contig for this rank

		if(($percas>=$minconsperc_asig) && ($perctotal>=$minconsperc_total) && ($totalcount>=$mingenes) && ($times>$times2)) { 
			$chimerism=($totalas-$times)/$totalas;
			# if($contig eq "k119_9990") {  printf "  $mtax $rank $accumtax{$mtax} $totalas $totalcount %.2f %.2f $chimerism\n",$percas,$perctotal; }	      
			$consensus{$rank}=$mtax; 
			printf outfile1 "   $rank: $totalas\t$mtax: $times\tChimerism: %.3f\n",$chimerism; 

			#-- Then we add the new taxon to the taxa found for previous ranks

			$cattax.="$mtax;";
			$fulltax.="$rank:$mtax;";
			$strg=$rank;
			$schim=$chimerism;
			$consensus=1;
			}	
								      
		# if($totalas!=$times) { print "***$contig $totalas $times\n"; }
		last if(!$consensus{$rank});	#-- And exit if there is no consensus for this rank
		}

	#-- Finally, write the output
		
	if(!$consensus) { $cattax="No consensus"; $strg="Unknown"; }			
	printf outfile2 "$contig\t$fulltax\t$strg\tChimerism level: %.3f\n",$chimerism;
                                 }
close outfile1;
close outfile2;
