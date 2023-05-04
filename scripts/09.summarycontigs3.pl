#!/usr/bin/env perl

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
my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$databasepath,$interdir,$tempdir,$mode,$fun3tax,$taxlist,$aafile,$alllog,$allorfs,$contigsfna,$rnafile,$fna_blastx,$doublepass,$euknofilter,$fun3tax_blastx,$mingenes9,$minconsperc_asig9,$minconsperc_total9,$consensus,$syslogfile);

#-- Some local configuration variables


my @ranks=('superkingdom','phylum','class','order','family','genus','species');
my @ranksabb=('k','p','c','o','f','g','s');
my(%validranks,%incluster);
map{ $validranks{$_}=1; } @ranksabb;

if($consensus) { $minconsperc_total9=$consensus; }  #-- Overrides the parameter if user-set or mode minion

open(syslogfile,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Reading taxonomic infomation (extracted from NCBI's taxonomy)

my(%taxcorr,%numorfs,%allcontigs,%orfs);
if(!$euknofilter) { $euknofilter="0"; }

open(infile1,$taxlist) || die "Can't open $taxlist\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @y=split(/\t/,$_);
	$taxcorr{$y[1]}=$y[2];
	}
close infile1;

if($mode eq "clustered") {	#-- Reading reporesentatives if in clustering mode
	my $clusterfile="$tempdir/allfaa.clus.rep.tsv";
	print "  Reading representatives from $clusterfile\n";
	open(inclus,$clusterfile) || die;
	while(<inclus>) {
		chomp;
		next if !$_;
		my($repr,$member)=split(/\t/,$_);
		$incluster{$repr}{$member}=1;
		}
	close(inclus);
	}


#-- Reading genes in contigs

tie %allcontigs,"Tie::IxHash";

print "  Reading contigs from $contigsfna\n";
open(infile2,$contigsfna) || die "Can't open $contigsfna\n";
while(<infile2>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { $allcontigs{$1}=1; }
	}
close infile2;


open(infile2,$aafile) || die "Can't open $aafile for writing\n";
while(<infile2>) {
	chomp;
	next if($_!~/^\>/);
	if($_=~/^\>([^ ]+)/) { 
		my $orf=$1;
		my $contig=$orf;
		$contig=~s/\_\d+\-\d+$//;
		$numorfs{$contig}++; 
		# $allcontigs{$contig}=1;
		$orfs{$contig}{$orf}=1;
		}
	}
close infile2;

if($doublepass) {
	open(infile2,$fna_blastx) || die "Can't open $fna_blastx\n";
	while(<infile2>) {
		chomp;
		next if($_!~/^\>/);
		if($_=~/^\>([^ ]+)/) { 
			my $orf=$1;
			my $contig=$orf;
			$contig=~s/\_\d+\-\d+$//;
			$numorfs{$contig}++; 
			# $allcontigs{$contig}=1;
			$orfs{$contig}{$orf}=1;
			}
		}
	close infile2;
	}

open(infile2,$rnafile) || warn "Can't open $rnafile\n";
while(<infile2>) {
	chomp;
	next if($_!~/^\>/);
	if($_=~/^\>([\w]+)/) { 
		my $orf=$1;
		my $contig=$orf; 
		$contig=~s/\_\d+\-\d+$//;
		$numorfs{$contig}++;
		$orfs{$contig}{$orf}=1;
		}
	}
close infile2;

my %parents;
open(infile0,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile0>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	#$tax=~s/\[|\]//g;
	$parents{$tax}{wranks}=$par;
	my @m=split(/\;/,$par);
	foreach my $y(@m) {
		my($rt,$gtax)=split(/\:/,$y);
		$parents{$tax}{noranks}.="$gtax;"; 
		}
	chop $parents{$tax}{noranks};
	}
close infile0;

#-- Input and output files

my $input;
my @taxfiles=('idfilter','noidfilter');
foreach my $tfile(@taxfiles) {
	my $prefix;
	if($tfile eq 'noidfilter') { $prefix=".noidfilter"; }
	if($doublepass) { $input="$fun3tax_blastx$prefix.wranks"; } else { $input="$fun3tax$prefix.wranks"; }
	my $outputlong=$allorfs.$prefix;
	my $outputshort=$alllog.$prefix;
	my $outputall="$interdir/09.$projectname.contiglog_allranks$prefix";


	#-- Reading taxonomic information for the genes (wrank file)

	my %taxlist;
	print "  Reading $input\n";
	print syslogfile "  Reading taxa for genes from $input\n";
	open(infile3,$input) || die "Can't open $input\n";
	while(<infile3>) {		#-- Looping on the ORFs
		my $contigid;
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my($node,$atax)=split(/\t/,$_);		#-- $node contains ORF name
		my @arep;
		print "--$node \n";
		foreach my $i(keys %{ $incluster{$node} }) { print "  -> $i\n"; }
		if($mode eq "clustered") { @arep=keys %{ $incluster{$node} }; } else { push(@arep,$node); }
		foreach my $clusrep(@arep) { 
			$atax=~s/\"//g;
			my @tf=split(/\;/,$atax);		#-- This contains rank:taxon pairs for the ORF
			my @id=split(/\_/,$clusrep);
			$contigid=$clusrep;
			$contigid=~s/\_\d+\-\d+$//;
			$contigid=~s/\|\d+$//; 		#-- This is the contig name the current ORF belongs to
		
			#-- Stores rank and taxa for all the ORFs in the contig

			foreach my $uc(@tf) { 
				 my ($rank,$tax)=split(/\_/,$uc);
				 next if(!$validranks{$rank});
				if($rank ne "n") { $taxlist{$contigid}{$rank}{$clusrep}=$tax; }
				}
			}
		}
	close infile3;

	if($euknofilter) {	#-- Remove filters for Eukaryotes
		my $eukinput=$input;
		$prefix=".noidfilter";
		if($doublepass) { $eukinput="$fun3tax_blastx$prefix.wranks"; } else { $eukinput="$fun3tax$prefix.wranks"; }

		#if($eukinput=~/noidfilter/) {	#-- Will only read the no id filter file for skipping filters for eukaryotes
			# $eukinput=~s/\.wranks/\.noidfilter\.wranks/;
			print syslogfile "  Reading results without eukaryotic filter from $eukinput\n";
			open(infile3,$eukinput) || die "Can't open $eukinput\n";
			while(<infile3>) {		#-- Looping on the ORFs
				my $contigid;
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my($node,$atax)=split(/\t/,$_);		#-- $node contains ORF name
				next if($atax!~/k\_Eukaryota/);
				$atax=~s/\"//g;
				my @arep;
				if($mode eq "clustered") { @arep=keys %{ $incluster{$node} }; } else { push(@arep,$node); }
				foreach my $clusrep(@arep) {
					my @tf=split(/\;/,$atax);		#-- This contains rank:taxon pairs for the ORF
					my @id=split(/\_/,$clusrep);
					# @id=split(/\|/,$node);
	
					#-- Different nomenclature for contigs in spades and other assemblers
	
					if($id[0]=~/NODE/) { $contigid="$id[0]\_$id[1]"; } else { $contigid="$id[0]"; }	

					$contigid=$clusrep;
					$contigid=~s/\_\d+\-\d+$//;
					$contigid=~s/\|\d+$//; 		#-- This is the contig name the current ORF belongs to
					#-- Stores rank and taxa for all the ORFs in the contig

					foreach my $uc(@tf) { 
						my ($rank,$tax)=split(/\_/,$uc);
						if($rank ne "n") { $taxlist{$contigid}{$rank}{$clusrep}=$tax;  }
						}
					}
				}	
			close infile3;
		}


	#-- Preparing output files

	print "  Writing output to $outputshort\n";
	print syslogfile "  Writing output to $outputlong and $outputshort\n";
	open(outfile1,">$outputlong") || die "Can't open $outputlong for writing\n";
	open(outfile2,">$outputshort") || die "Can't open $outputshort for writing\n";
	open(outfile3,">$outputall") || die "Can't open $outputall for writing\n";
	print outfile1 "#- Created by $0 with data from $input, mingenes=$mingenes9, minconsperc_asig=$minconsperc_asig9, minconsperc_total=$minconsperc_total9, euknofilter=$euknofilter, ",scalar localtime,"\n";
	print outfile2 "#- Created by $0 with data from $input, mingenes=$mingenes9, minconsperc_asig=$minconsperc_asig9, minconsperc_total=$minconsperc_total9, euknofilter=$euknofilter,",scalar localtime,"\n";
	print outfile3 "#- Created by $0 with data from $input, euknofilter=$euknofilter,",scalar localtime,"\n";
	print outfile3 "# Contig\tRank\tPerc_assigned\tPerc_total\tGene number\tDisparity\tTax\n";

	foreach my $contig(keys %allcontigs) { 
		# print "***$contig\n";
		my ($sep,$lasttax,$strg,$cattax,$fulltax,$lasttax)="";
		my ($consensus,$schim,$chimerism)=0;
		my(%consensus,%chimeracheck,%consensusall)=();
		my $totalcount=$numorfs{$contig}; 
		next if(!$totalcount);	
		print outfile1 "$contig\t$totalcount genes\n";
	
		#-- We loop on ranks
	
		foreach my $rank(@ranksabb) {
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
			if(!$totalcount) { print "  Warning: Contig $contig has zero genes\n"; }
			my $percas=$accumtax{$mtax}/$totalas;
			my $perctotal=$accumtax{$mtax}/$totalcount;
			if(!$times2) { $times2=0; }
			 # if($contig=~/278124/) {  printf "  $mtax $rank $accumtax{$mtax} $totalas $totalcount %.2f %.2f\n",$percas,$perctotal; }	      

			#-- And if it does, we also calculate the disparity index of the contig for this rank

			# print "$contig\t$rank\t$mtax\t$percas\t$perctotal\t$totalcount\n";

			if(($percas>=$minconsperc_asig9) && ($perctotal>=$minconsperc_total9) && ($totalcount>=$mingenes9) && ($times>$times2)) { 
	
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
						}
				}


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
			$consensusall{$rank}{tax}=$mtax; 
			$consensusall{$rank}{percas}=$percas; 
			$consensusall{$rank}{perctotal}=$perctotal;
			$consensusall{$rank}{totalcount}=$totalcount;
			$consensusall{$rank}{disparity}=$chimerism;
		
								      
			# if($totalas!=$times) { print "***$contig $totalas $times\n"; }
			# last if(!$consensus{$rank});	#-- And exit if there is no consensus for this rank
			}

		#-- Finally, write the output
		
		if(!$consensus) { $cattax="No consensus"; $strg="Unknown"; }		
		my $abb=$parents{$lasttax}{wranks};	
		$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
		# printf outfile2 "$contig\t$fulltax\t$strg\tDisparity: %.3f\tGenes: $numorfs{$contig}\n",$chimerism;
		printf outfile2 "$contig\t$abb\t$strg\tDisparity: %.3f\tGenes: $numorfs{$contig}\n",$chimerism;
		foreach my $rank(@ranksabb) {
			my $ttax=$consensusall{$rank}{tax};
			if($ttax) { 
				my $abb=$parents{$ttax}{wranks};	
				$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
				printf outfile3 "$contig\t$rank\t%.1f\t%.1f\t$consensusall{$rank}{totalcount}\t%.1f\t$abb\n",$consensusall{$rank}{percas},$consensusall{$rank}{perctotal},$consensusall{$rank}{disparity};
				}
			}
		}
	close outfile1;
	close outfile2;
	close outfile3;
}
close syslogfile;
