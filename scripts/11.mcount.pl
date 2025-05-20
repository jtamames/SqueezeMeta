#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Counts the sizes for all taxa in the contiglog file created by summarycontigs3

use strict;
use Cwd;
use lib ".";

use File::Basename;
use Cwd 'abs_path';
our $scriptdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $scriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $scriptdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$scriptdir/..");

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

	#-- Configuration variables from conf file

our($datapath,$resultpath,$contigslen,$alllog,$taxlist,$contigcov,$mcountfile,$syslogfile);

my(%lon,%taxa,%abund,%abundreads,%samples,%accum,%accumbases,%accumreads,%taxcorr,%cseen);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

	#-- Read contig lengths

print outsyslog "Reading contig length from $contigslen\n";
open(infile1,$contigslen) || die "Can't open $contigslen for reading\n";
while(<infile1>) {  
	chomp;
	next if !$_;
	my ($node,$len)=split(/\t/,$_);
	$lon{$node}=$len;
	}
close infile1;

	#-- Read contiglog file to get taxonomic assignment for contigs

print outsyslog "Reading contig taxa from $alllog\n";
open(infile2,$alllog) || die "Can't open $alllog for writing\n";
while(<infile2>) {
	chomp;
	next if !$_;
	my($node,$tax,$rest)=split(/\t/,$_);
	if($tax eq "No consensus") { $tax="Unknown"; }
	if(!$tax) { $tax="Unknown"; }
	$taxa{$node}=$tax;
	my $string="";
	my @tx=split(/\;/,$tax);
	for(my $n=0; $n<=$#tx; $n++) {
		$string.="$tx[$n];";
		$accum{$string}+=$lon{$node};
		}
	$cseen{$node}=1;	
	}
close infile2;

		#-- Contigs with no genes are not in the contiglog file, they must be counted apart

foreach my $acg(keys %lon) {
	if(!$cseen{$acg}) { $accum{"Unknown;"}+=$lon{$acg}; }
	}

	#-- Read contigcov file to get abundances of each contig

print outsyslog "Reading contig coverages from $contigcov\n";
open(infile3,$contigcov) || die "Can't open $contigcov for writing\n";
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @f=split(/\t/,$_);
	my $sample=$f[$#f];
	$abund{$f[0]}{$sample}=$f[6];
	$abundreads{$f[0]}{$sample}=$f[5];
	$samples{$sample}=1; 
	my $node=$f[0];
	my $tlong=$lon{$node};
	my $tax=$taxa{$node};
	if(!$tax) { $tax="n"; }
	my @tx=split(/\;/,$tax);
	my $string="";
	
	#-- For all the ranks of the current contig, add the contig size to the corresponding taxon
	
	for(my $n=0; $n<=$#tx; $n++) {
		$string.="$tx[$n];";
		# $accum{$string}+=$tlong;
		
		#-- Add also bases and reads
		
		$accumbases{$string}{$sample}+=$abund{$node}{$sample};
		$accumreads{$string}{$sample}+=$abundreads{$node}{$sample}; 
 		}
 
	#print "$f[0] $sample $f[5]\n";
	}
close infile3;

	#-- Read the equivalence between taxa and rank

open(infile4,$taxlist) || die "Can't open $taxlist for writing\n";
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

print outsyslog "Writing output to $mcountfile\n"; 
open(outfile1,">$mcountfile") || die "Can't open $mcountfile for writing\n";
print outfile1 "Rank\tTaxon\tAccumulated contig size";
foreach my $samp(sort keys %samples) { print outfile1 "\t$samp reads\t$samp bases"; }
print outfile1 "\n";
foreach my $kk(sort { $accum{$b}<=>$accum{$a}; } keys %accum) { 
	my $k=$kk;
	$k=~s/\;$//;
	my @l=split(/\;/,$k);
	my($rank,$tn)=split(/\_/,$l[$#l]);
	if((!$rank) || ($rank eq "Unknown")) { $rank="n"; }				  
	print outfile1 "$rank\t$k\t$accum{$kk}"; 
	foreach my $samp(sort keys %samples) { print outfile1 "\t$accumreads{$kk}{$samp}\t$accumbases{$kk}{$samp}"; }
	print outfile1 "\n";
	}
close outfile1;
close outsyslog;
