#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs MaxBin for binning

use strict;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($databasepath,$contigsfna,%bindirs,$contigcov,$maxbin_soft,$alllog,$contigslen,$interdir,$singletons,$tempdir,$numthreads,$mappingfile,$methodsfile,$syslogfile);

my $maxchimerism=0.1;	#-- Threshold for excluding chimeric contigs
my $mingenes=1;		#-- Threshold for excluding small contigs (few genes than this)
my $smallnoannot=1;	#-- For excluding contigs with just one gene an no annotation

my(%allcontigs,%skip);

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


my %singletonlist;
if($singletons) {		#-- Excluding singleton raw reads from binning
	my $singletonlist="$interdir/01.$projectname.singletons";
	print "  Excluding singleton reads from $singletonlist\n";
	print outsyslog "  Excluding singleton reads from $singletonlist\n";
	open(infile0,$singletonlist) || die "Cannot open singleton list in $singletonlist\n";
	while(<infile0>) {
		chomp;
		next if !$_;
		my @y=split(/\t/,$_);
		$singletonlist{$y[0]}=1;
		}
	close infile0;
	}

print "  Reading samples from $mappingfile\n";   #-- We will exclude samples with the "noassembly" flag
open(infile0,$mappingfile) || die "Can't open $mappingfile\n";
while(<infile0>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	if($_=~/nobinning/) { $skip{$t[0]}=1; }
	}
close infile0;


if(-e $alllog) {
        open(infile1,$alllog) || die "Can't open $alllog\n";
        while(<infile1>) {
                chomp;
                next if !$_;
                my @r=split(/\t/,$_);
                next if($singletonlist{$r[0]});
                my($chimlevel,$numgenes);
                if($r[3]=~/Disparity\: (.*)/) { $chimlevel=$1; }
                if($r[4]=~/Genes\: (.*)/) { $numgenes=$1; }
                if(!$numgenes) { $numgenes=0; }
                if(($numgenes>=$mingenes) && ($chimlevel<=$maxchimerism)) { $allcontigs{$r[0]}=1; }
                if($smallnoannot && ($numgenes<=1) && ($r[1] eq "Unknown")) { delete $allcontigs{$r[0]}; }
                }
        close infile1;
} else {
        # if we don't have alllog (bc we did not run annotation) then get the list of contigs from $contigslen, which is generated at assembly
        #  so we don't do any filtering by disparity or numgenes
        print outsyslog "\nWARNING: $alllog was not found. We will not filter contigs by number of genes or chimerism!\n\n";
        open(infile1,$contigslen) || die "Can't open $contigslen\n";
        while(<infile1>) {
                chomp;
                next if !$_;
                my @r=split(/\t/,$_);
                next if($singletonlist{$r[0]}); # we can still avoid singletons so we do it
                $allcontigs{$r[0]}=1;
                }
        }



my $tempfasta="$tempdir/bincontigs.fasta";
open(outfile1,">$tempfasta") || die "Can't open $tempfasta for writing\n";
open(infile1,$contigsfna) || die "Can't open $contigsfna for writing\n";
my $ingood=0;
while(<infile1>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { 
		my $tc=$1;
		if($allcontigs{$tc}) { $ingood=1; } else { $ingood=0; } 
		if($singletonlist{$tc}) { $ingood=0; }
		}
	if($ingood) { print outfile1 "$_\n"; }
	}
close infile1;
close outfile1; 

	#-- Creating binning directory

my $dirbin="$interdir/binners/maxbin";
if(-d $dirbin) {} else { system "mkdir $dirbin"; }

	#-- Creating abundance file

my $abundlist="$dirbin/abund.list";
print outsyslog "Creating abundance file in $abundlist\n";
open(outfile1,">$abundlist") || die "Can't open $abundlist for writing\n";	#-- Stores the list of files for the abundance of contigs (one per sample)
open(infile2,$contigcov) || die "Can't open contig coverage file $contigcov\n";
my $currsample;
my %tcontigs;
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	next if($skip{$k[$#k]});
	if(!$currsample || ($k[$#k] ne $currsample)) {	#-- If we finished reading one sample, write the data
		foreach my $k(sort keys %tcontigs) { print outfile2 "$k\t0\n"; }
		close outfile2;
		$currsample=$k[$#k];
		print outfile1 "$dirbin/$currsample.abund\n";
		open(outfile2,">$dirbin/$currsample.abund") || die "Can't open $dirbin/$currsample.abund for writing\n"; #-- Stores the abundances for the current sample
		%tcontigs=%allcontigs;
  
		}
	if($allcontigs{$k[0]}) { print outfile2 "$k[0]\t$k[1]\n"; }
	delete $tcontigs{$k[0]};
	}
close infile2;	   
foreach my $k(sort keys %tcontigs) { print outfile2 "$k\t0\n"; }
close outfile2;
close outfile1;

my $command="perl $maxbin_soft -thread $numthreads -contig $tempfasta -abund_list $abundlist -out $dirbin/maxbin -markerpath $databasepath/marker.hmm";
print outsyslog "Running Maxbin: $command\n";
print "  Running Maxbin (Wu et al 2016, Bioinformatics 32(4), 605-7)\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

print outmet "Binning was done using MaxBin2 (Wu et al 2016, Bioinformatics 32(4), 605-7)\n";
close outmet;
close outsyslog;
