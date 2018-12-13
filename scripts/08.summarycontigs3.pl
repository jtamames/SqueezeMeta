#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Assigns a consensus annotation for the whole contigs using the annotations for the genes it contains
#-- Also calculates disparity index for each contig
#-- 23/05/2018, corrected calculation of disparity

use strict;
use Tie::IxHash;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$fun3tax,$taxlist,$aafile,$contigsfna,$rnafile);

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

my(%taxcorr,%numorfs,%allcontigs,%orfs);

open(infile1,$taxlist) || die;
while(<infile1>) {
	chomp;
	next if !$_;
	my @y=split(/\t/,$_);
	$taxcorr{$y[1]}=$y[2];
	}
close infile1;

#-- Reading genes in contigs

tie %allcontigs,"Tie::IxHash";

open(infile2,$contigsfna) || die;
while(<infile2>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { $allcontigs{$1}=1; }
	}
close infile2;


open(infile2,$aafile) || die;
while(<infile2>) {
	chomp;
	next if($_!~/^\>/);
	if($_=~/^\>([^ ]+)/) { 
		my $orf=$1;
		my $contig=$orf;
		$contig=~s/\_\d+$//;
		$numorfs{$contig}++; 
		# $allcontigs{$contig}=1;
		$orfs{$contig}{$orf}=1;
		}
	}
close infile2;

open(infile2,$rnafile) || die;
while(<infile2>) {
	chomp;
	next if($_!~/^\>/);
	if($_=~/^\>([\w]+)/) { 
		my $orf=$1;
		my $contig=$orf; 
		$contig=~s/\_RNA\d+$//;
		$numorfs{$contig}++;
		$orfs{$contig}{$orf}=1;
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
	next if(!$_ || ($_=~/^\#/));
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
		if($rank ne "no rank") { $taxlist{$contigid}{$rank}{$node}=$tax;  }
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
	my(%consensus,%chimeracheck)=();
	my $totalcount=$numorfs{$contig}; 	
	print outfile1 "$contig\t$totalcount genes\n";
	
	#-- We loop on ranks
	
	foreach my $rank(@ranks) {
		my $totalas=0;
		my %accumtax=();
		
		#-- And count how many times each taxon appears
		
		foreach my $orf(sort keys %{ $taxlist{$contig}{$rank} }) { 
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

		#-- And if it does, we also calculate the disparity index of the contig for this rank

		if(($percas>=$minconsperc_asig) && ($perctotal>=$minconsperc_total) && ($totalcount>=$mingenes) && ($times>$times2)) { 
	
			#-- Calculation of disparity for this rank
			my($chimera,$nonchimera,$unknown)=0;
			foreach my $orf(sort keys %{ $orfs{$contig} }) { 
				my $ttax=$taxlist{$contig}{$rank}{$orf}; 
				foreach my $orf2(sort keys %{ $orfs{$contig} }) { 
					my $ttax2=$taxlist{$contig}{$rank}{$orf2}; 
					next if($orf ge $orf2);
					if($chimeracheck{$orf}{$orf2} eq "chimera") {	#-- If it was a chimera in previous ranks, it is a chimera now
					#	if($contig eq "3539") { print "$rank $orf $orf2 -> Previous chimera\n"; }
						$chimera++;
						next;
					}	
					if($chimeracheck{$orf}{$orf2} eq "unknown") {	#-- If it was an unknown in previous ranks, it is a unknown now
					#	 if($contig eq "3539") { print "$rank $orf $orf2 -> Previous unknown\n"; }
						$unknown++;
						next;
					}						
					if(($ttax && (!$ttax2)) || ($ttax2 && (!$ttax)) || ((!$ttax2) && (!$ttax))) { $chimeracheck{$orf}{$orf2}="unknown"; } #-- Unknown when one of the ORFs has no classification at this rank
					elsif($ttax eq $ttax2) { $chimeracheck{$orf}{$orf2}="nochimera"; $nonchimera++; }
					else { $chimeracheck{$orf}{$orf2}="chimera"; $chimera++; }
					# if($contig eq "3539") { print "$rank $orf $orf2 -> $ttax $ttax2 -> $chimeracheck{$orf}{$orf2}\n"; }
					}
			}


			# $chimerism=($totalas-$times)/$totalas;
			my $totch=$chimera+$nonchimera;
			if($totch) { $chimerism=$chimera/($chimera+$nonchimera); } else { $chimerism=0; }
			# print "$chimerism*****\n";
			# if($contig eq "k119_9990") {  printf "  $mtax $rank $accumtax{$mtax} $totalas $totalcount %.2f %.2f $chimerism\n",$percas,$perctotal; }	      
			$consensus{$rank}=$mtax; 
			printf outfile1 "   $rank: $totalas\t$mtax: $times\tDisparity: %.3f\n",$chimerism; 

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
	printf outfile2 "$contig\t$fulltax\t$strg\tDisparity: %.3f\tGenes: $numorfs{$contig}\n",$chimerism;
                                 }
close outfile1;
close outfile2;
