#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Adds taxonomic assignments for bins, according to the consensus of the contigs belonging to it

$|=1;

use strict;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;

do "$project/SqueezeMeta_conf.pl";

	#-- Configuration variables from conf file

our($datapath,$alllog,$bintax,%bindirs,%dasdir);

	#-- Some configuration values for the algorithm
	
my $mincontigs=1;  		#-- Minimum number of contigs for the bin 
my $minconsperc_asig=0.6;  	#-- Ratio contigs for the taxon/sum(genes all taxa). Therefore it only considers assigned contigs
my $minconsperc_total=0.3;  	#-- Ratio contigs for the taxon/number of contigs. Therefore it considers all (assigned+unassigned) contigs

my $verbose=0;

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
foreach my $binmethod(sort keys %dasdir) {		#-- For the current binning method
	my $bindir=$dasdir{$binmethod};
	print "Looking for $binmethod bins in $bindir\n";

	#-- Reading bin directories

	opendir(indir1,$bindir);
	my @files=grep(/fa$|fasta$/,readdir indir1);
	closedir indir1;
	
	#-- Working with each bin
	
	foreach my $k(@files) { 
		my $tf="$bindir/$k";
		my $outf="$bindir/$k.tax";
 		# next if($k ne "maxbin.002.fasta");
		
		#-- We will create an output file for bin, with the contig names and consensus
		
		open(outfile2,">$outf") || die;			
		open(infile2,$tf) || die "Cannot open $tf\n";
		my $size=0;
		my %store=();
		
		#-- Read the contigs in the bin, with assignments and sizes
		
		while(<infile2>) {
			chomp;
			if($_=~/\>([^ ]+)/) {
				my $contigid=$1;
				my $taxid=$tax{$contigid};
				print outfile2 "$_ $taxid\n";
				$store{$contigid}=1;
				}
			else { chomp; $size+=length $_; }
		}
		close infile2;
		
		#-- Call consensus() to find the consensus taxon
		
		my %chimeracheck;
		my($sep,$lasttax,$strg,$fulltax,$cattax)="";
		my $chimerism=0;
	
		#-- Loop for all ranks
	
		foreach my $rank(@ranks) { 
			print ">>$rank\n" if $verbose;
			my($totalcount,$totalas,$consf)=0;
			my %accumtax=();
		
			#-- For all contigs in the bin, count the taxa belonging to that rank
		
			foreach my $contig(keys %store) { 
				my $ttax=$taxlist{$contig}{$rank}; 
				$totalcount++;
				next if(!$ttax);
				$accumtax{$ttax}++;
				$totalas++; 
				# if($k eq "maxbin.002.fasta") { print "   $contig $ttax $accumtax{$ttax}\n"; }
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
				#-- Calculation of disparity for this rank
				my($chimera,$nonchimera,$unknown)=0;
				foreach my $contig(sort keys %store) { 
					my $ttax=$taxlist{$contig}{$rank};
					foreach my $contig2(sort keys %store) { 
						my $ttax2=$taxlist{$contig2}{$rank}; 
						next if($contig ge $contig2);
						if($chimeracheck{$contig}{$contig2} eq "chimera") {	#-- If it was a chimera in previous ranks, it is a chimera now
							$chimera++;
							next;
						}	
						if($chimeracheck{$contig}{$contig2} eq "unknown") {	#-- If it was an unknown in previous ranks, it is a unknown now
							$unknown++;
							next;
						}						
						if(($ttax && (!$ttax2)) || ($ttax2 && (!$ttax)) || ((!$ttax2) && (!$ttax))) { $chimeracheck{$contig}{$contig2}="unknown"; } #-- Unknown when one of the ORFs has no classification at this rank
						elsif($ttax eq $ttax2) { $chimeracheck{$contig}{$contig2}="nochimera"; $nonchimera++; }
						else { $chimeracheck{$contig}{$contig2}="chimera"; $chimera++; }
						# if($contig eq "3539") { print "$rank $orf $orf2 -> $ttax $ttax2 -> $chimeracheck{$orf}{$orf2}\n"; }
						}
				}
				my $totch=$chimera+$nonchimera;
				if($totch) { $chimerism=$chimera/($chimera+$nonchimera); } else { $chimerism=0; }
				print "***$mtax $times $percas $perctotal $totalas $chimerism\n" if $verbose;
				$cattax.="$mtax;";
				$fulltax.="$rank:$mtax;";
				$consf=1;
				$strg=$rank;
			#	 if($consf && ($k eq "maxbin.002.fasta")) { print "**$fulltax $percas $perctotal -- $times $totalas $totalcount\n"; }
				 print "**$fulltax $percas $perctotal -- $times $totalas $totalcount -- $consf\n" if $verbose; 
				 if($totalas!=$times) { print "***$k$totalas $times\n" if $verbose; }
				}
			 last if(!$consf);
			}
		if(!$fulltax) { $fulltax="No consensus"; $strg="Unknown"; }
		
		
		
		
		
		
		#-- Write output
		
		printf outfile2 "Consensus: $fulltax\tTotal size: $size\tDisparity: %.3f\n",$chimerism;
		printf outfile1 "$binmethod\t$k\tConsensus: $fulltax\tTotal size: $size\tDisparity: %.3f\n",$chimerism;
		close outfile2;

 	}
}
close outfile1;


