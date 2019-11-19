#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

#-- This program prepares sqm output for being loaded in pavian

$|=1;

local $; = '#';
use strict;
use Cwd;
use lib ".";

my $pwd=cwd();
my $feature=$ARGV[1];
my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

our($resultpath,$databasepath,$mcountfile);

my %rankequival=('k','D','p','P','c','C','o','O','f','F','g','G','s','S');
my %spaces=('D',2,'P',4,'C',6,'O',8,'F',10,'G',12,'S',14);
my(%totalreads,%ncbitax);

if(!$feature) { $feature="reads"; }
my $parentfile="$databasepath/LCA_tax/parents.txt";
open(in,$parentfile) || die "Cannot open $parentfile\n";
while(<in>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	$ncbitax{$k[0]}=$k[2];
	}
close in;

open(in,"$resultpath/10.$project.mappingstat") || die;
while(<in>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	if($feature eq "reads") { $totalreads{$k[0]}=$k[1]; }
	elsif($feature eq "bases") { $totalreads{$k[0]}=$k[4]; }
	else { die "Bad feature: Use reads or bases\n"; }
	print "*$k[0]*$totalreads{$k[0]}*\n";
	}
close in;

open(in,$mcountfile) || die;
my($header,$sample,$r);
my(%store,%abundance,%accum);
my @header;
while(<in>) {
	chomp;
	next if !$_;
	if(!$header) { 
		$header=$_;
		@header=split(/\t/,$header);
		next;
		}
	my @fields=split(/\t/,$_);
	my $rank=$fields[0];
	my $tax=$fields[1];
	for(my $pos=2; $pos<=$#fields; $pos++) {
		if($header[$pos]=~/ $feature$/) {
			($sample,$r)=split(/\s+/,$header[$pos]);	 
			my $abun=$fields[$pos];
			$store{$sample}{$tax}=1;
			$abundance{$sample}{$tax}=$abun;
			my @tx=split(/\;/,$tax);
			pop @tx;
			my $parent=join(";",@tx);
			$accum{$sample}{$parent}+=$abun;
			}
		}    
	}
close in;

foreach my $insample(sort keys %store) {
	my $outfile="$insample.pavian";
	open(outfile1,">$outfile") || die;
	my $tabun_unclas=$abundance{$insample}{'Unknown'};
	my $tperc=($tabun_unclas/$totalreads{$insample})*100;
	next if($tperc<0.01);
	# printf outfile1 " %.2f\%\t$tabun_unclas\t$tabun_unclas\tU\t0\tunclassified\n",$tperc;
	my $tabun_root=$totalreads{$insample}-$tabun_unclas;
	my $tperc=($tabun_root/$totalreads{$insample})*100;
	next if($tperc<0.01);
	printf outfile1 " %.2f\%\t$tabun_root\t$tabun_root\tR\t1\troot\n",$tperc;


	my($lastblank,$thisblank);
	foreach my $p(sort keys %{ $store{$insample} }) { 
		my @tx=split(/\;/,$p);
		my $tt=$tx[$#tx];
		my($trank,$thistax)=split(/\_/,$tt);
		next if($tt eq "Unknown");
		my $tabun=$abundance{$insample}{$p};
		my $inthisnode=$abundance{$insample}{$p}-$accum{$insample}{$p};
		my $tperc=($tabun/$totalreads{$insample})*100;
		next if($tperc<0.01);
		# print "$trank\t$thistax\t$tabun\t$inthisnode\t$tperc\n"; 
		my $nrank=$rankequival{$trank};
		my $tax_id=$ncbitax{$thistax};
		if(!$tax_id) { $tax_id=1; }
		$thisblank=$spaces{$nrank};	#-- All this stupid stuff is needed to accomodate the suboptimal native format for input data (kraken report format)
		# print "*$thistax*$thisblank*$lastblank*\n";
		if($thisblank>$lastblank+2) { $thisblank=$lastblank+2; }
		my $spblank=" " x $thisblank;
		$lastblank=$thisblank;
		printf outfile1 " %.2f\%\t$tabun\t$inthisnode\t$nrank\t$tax_id\t$spblank$thistax\n",$tperc;   
		}
	print "Output in $feature for sample $insample in $outfile\n";	
	}
print "You can now use these data in your R pavian or in the web app https://fbreitwieser.shinyapps.io/pavian\n";
