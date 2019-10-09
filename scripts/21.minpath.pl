#!/usr/bin/env perl

use strict;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;

if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

#-- Configuration variables from conf file

our($extdatapath,$contigsinbins,$mergedfile,$tempdir,$resultpath,$minpath_soft,$bintable,$minfraction21,%bindirs,%dasdir);
my(%pathid,%ec,%ecs,%kegg,%inbin,%bintax);

open(infile1,"$extdatapath/metacyc_pathways_onto.txt") || die "Can't open $extdatapath/metacyc_pathways_onto.txt\n";
while(<infile1>) { 
	chomp;
	next if !$_;
	my @k=split(/\t/,$_);
	my($pathw,$onto)=split(/ -> /,$k[1]); 
	$pathid{$k[0]}=$pathw;
	}
close infile1; 

open(infile2,"$extdatapath/kegg2ec.txt") || die "Can't open $extdatapath/kegg2ec.txt\n";
while(<infile2>) {
	chomp;
	next if !$_;
	my @k=split(/\t/,$_);
	$k[1]=~s/^EC\://;
	$ec{$k[0]}{$k[1]}=1;
	}
close infile2;

open(infile3,$contigsinbins) || die "Can't open $contigsinbins\n";
while(<infile3>) {
	chomp;
	next if !$_;
	my @k=split(/\t/,$_);
	$inbin{$k[0]}{$k[2]}=$k[1];
	}
close infile3;

open(infile4,$bintable) || die "Can't open $bintable\n";
while(<infile4>) {
	chomp;
	next if !$_;
	my @k=split(/\t/,$_);
	$bintax{$k[0]}=$k[2];
	}
close infile4;

my($header,$keggpos);
open(infile5,$mergedfile)  || die "Can't open $mergedfile\n";
while(<infile5>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	if(!$header) {
		$header=$_;
		my @ht=split(/\t/,$header);
		for(my $cpos=0;$cpos<=$#ht; $cpos++) { 
			if($ht[$cpos] eq "KEGG ID") { $keggpos=$cpos; last; }
			}
		next;
	}
	my @u=split(/\t/,$_);
	my $keggid=$u[$keggpos];
	$keggid=~s/\*//g; # Bug fix, removing * from Kegg annotation in ORF table.
	my $contigid=$u[1];
	foreach my $ibin(keys %{ $inbin{$contigid} }) { $kegg{$ibin}{$keggid}++; }
	foreach my $thisec(keys %{ $ec{$keggid} }) {
		foreach my $ibin(keys %{ $inbin{$contigid} }) { $ecs{$ibin}{$thisec}++; } 
		}
	     }
close infile5;

my(%pathways,%allpaths,%totalpathways,%pathgenes);
kegg();
outres("kegg");

my(%pathways,%allpaths,%totalpathways,%pathgenes);
metacyc();
outres("metacyc");

sub outres {
	my $clas=shift;
	my $totpath;
	open(outfile5,">$resultpath/21.$project.$clas.pathways") || die "Can't open $resultpath/21.$project.$clas.pathways for writing\n";
	print outfile5 "Bin\tTax\tPathways found";
	foreach my $pt(sort keys %allpaths) { print outfile5 "\t$pt"; $totpath++; }
	print outfile5 "\n";
	print outfile5 "\t\t$totpath";
	foreach my $pt(sort keys %allpaths) { print outfile5 "\t$pathgenes{total}{$pt}"; }
	print outfile5 "\n";
	foreach my $mt(sort keys %pathways) {
		print outfile5 "$mt\t$bintax{$mt}\t$totalpathways{$mt}";
		foreach my $pt(sort keys %allpaths) { 
			my $fraction=$pathgenes{$mt}{$pt}/$pathgenes{total}{$pt};
			my $pv;
			if((($fraction>=$minfraction21) || ($pathgenes{$mt}{$pt}>=5)) && ($pathways{$mt}{$pt})) { $pv=$pathgenes{$mt}{$pt} } else { $pv="NF"; }
			print outfile5 "\t$pv"; 
			}
		print outfile5 "\n";
		}
	close outfile5;
	}
 
sub metacyc {
	foreach my $kbin(sort keys %ecs) {
		my $outec="$tempdir/$kbin.minpath.temp";
		open(outfile1,">$outec") || die "Can't open $outec for writing\n";
		my $id=0;
		foreach my $ecbin(sort keys %{ $ecs{$kbin} }) {
		       next if !$ecbin;
		       next if($ecbin=~/\-/);
		       $ecbin=~s/\*//g;
		       $id++;
		       print outfile1 "read$id\t$ecbin\n";
		       }
		close outfile1; 
		print "Running MinPath for metacyc: $kbin         \r";
		my $command="$minpath_soft -any $outec -map ec2path -report $tempdir/$kbin.minpath.temp.report -details $tempdir/$kbin.metacyc.details  > /dev/null";
		my $ecode = system $command;
 		if($ecode!=0) {
			print "WARNING: Error running command:    $command\n";
			print "This is normally due to a small bin having no annotations recognizable by MinPath.\n"
			}

		open(infile5,"$tempdir/$kbin.minpath.temp.report") || next;
		my %accum=();
		while(<infile5>) {
			chomp;
			if($_=~/minpath 1/) {
				my @k=split(/\s+/,$_);
				my $thisp=$k[$#k];
				my $thisonto=$pathid{$thisp};
				if($thisonto) { $accum{$thisonto}=1; }
				}
			}
		close infile5;  
		open(infile7,"$tempdir/$kbin.metacyc.details") || next;
			while(<infile7>) {
			chomp;
				if($_=~/^path.*fam0 (\d+) fam-found (\d+) \# (.*)/) {
				my $pname=$pathid{$3};
				$pathgenes{total}{$pname}=$1;
				$pathgenes{$kbin}{$pname}=$2;
				}
			}
		close infile7;
	
	
		open(outfile4,">$tempdir/$kbin.metacyc.pathways") || die "Can't open $tempdir/$kbin.metacyc.pathways for writing\n";
		foreach my $konto(sort keys %accum) { 
			print outfile4 "$konto\n"; 
			$pathways{$kbin}{$konto}=1; 
			$allpaths{$konto}=1;
			$totalpathways{$kbin}++;
			}  
		close outfile4;     
		}
	print "\n";				 
	}

sub kegg {
	foreach my $kbin(sort keys %kegg) {
		my $outkegg="$tempdir/$kbin.minpath.temp.kegg";
		my($binmethod,$rest)=split(/\./,$kbin);
		my $outdir="$resultpath/$binmethod";
		open(outfile3,">$outkegg") || die "Can't open $outkegg for writing\n";
		my $id=0;
		foreach my $keggbin(sort keys %{ $kegg{$kbin} }) {
			next if !$keggbin;
			$keggbin=~s/\*//g;
			$id++;
			print outfile3 "read$id\t$keggbin\n";
			}
		close outfile3;	
		print "Running MinPath for kegg: $kbin         \r";
		my $command="$minpath_soft -ko $outkegg -map ec2path -report $tempdir/$kbin.minpath.temp.report -details $outdir/$kbin.kegg.details > /dev/null";
		my $ecode = system $command;
 		if($ecode!=0) {
			print "WARNING: Error running command:    $command\n";
			print "This is normally due to a small bin having no annotations recognizable by MinPath.\n"
			}
		open(infile6,"$tempdir/$kbin.minpath.temp.report") || next;
		my %accum=();
		my $pathname;
		while(<infile6>) {
			chomp;
			if($_=~/minpath 1/) {
				if($_=~/name\s+(.*)/) { $pathname=$1; }
				my @k=split(/\s+/,$_);
				my $thisp=$k[$#k];
				my $thisonto=$pathid{$thisp};
				# $accum{$thisonto}=1;
				$accum{$pathname}=1;
				}
			}
		close infile6;	
	
		open(infile7,"$outdir/$kbin.kegg.details") || next;
		while(<infile7>) {
			chomp;
			if($_=~/^path.*fam0 (\d+) fam-found (\d+) \# (.*)/) {
				$pathgenes{total}{$3}=$1;
				$pathgenes{$kbin}{$3}=$2;
				}
			}
		close infile7;
	
	
		open(outfile4,">$outdir/$kbin.kegg.pathways") || die "Can't open $outdir/$kbin.kegg.pathways for writing\n";
		foreach my $konto(sort keys %accum) { 
			print outfile4 "$konto\n"; 
			$pathways{$kbin}{$konto}=1; 
			$allpaths{$konto}=1;
			$totalpathways{$kbin}++;
			}  
		close outfile4;     
		}
	print "\n";
	}				 



