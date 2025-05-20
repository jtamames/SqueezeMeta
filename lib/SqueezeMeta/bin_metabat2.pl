#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs binning with Metabat2

use strict;
use Cwd;
use lib ".";

use File::Basename;
use Cwd 'abs_path';
our $sqmlibdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $sqmlibdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $sqmlibdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$sqmlibdir/../..");

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;
do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($contigsfna,$contigcov,$bindir,$metabat_soft,$jgi_summ_soft,$alllog,$datapath,$tempdir,$interdir,$singletons,$contigslen,$mappingfile,$methodsfile,$maxchimerism14,$mingenes14,$smallnoannot14,%bindirs,$syslogfile,$numthreads);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


#-- Creating binning directory

my $dirbin="$interdir/binners/metabat2";
if(-d $dirbin) {} else { system "mkdir $dirbin"; }


#-- Calculate average coverages and variances for each contig

# Exclude samples with the nobinning flag
my %samples;
my %skip;
print "  Reading samples from $mappingfile\n";
open(infile0,$mappingfile) || die "Can't open $mappingfile\n";
while(<infile0>) {
        chomp;
        next if !$_;
        my @t=split(/\t/,$_);
	$samples{$t[0]}=1;
        if($_=~/nobinning/) { $skip{$t[0]}=1; }
        }
close infile0;

my $bamlist = "";
foreach my $sample (keys %samples) {
	next if($skip{$sample});
	$bamlist = "$bamlist $datapath/bam/$projectname.$sample.bam";
}
if(!$bamlist) { die "All samples have the \"nobinning\" flag so there are no valid BAM files. Please check your samples file"; }

# Call jgi_summarize_bam_contig_depths
my $depthfile="$dirbin/contigs.depth.txt";
my $depthfileprov="$depthfile.PROV";
print "  Creating coverage file in $depthfile\n";
print outsyslog "Creating coverage file in $depthfile\n";
my $command = "$jgi_summ_soft $bamlist --outputDepth $depthfileprov >> $syslogfile 2>&1";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }


#-- Exclude singleton raw reads from binning

my %singletonlist;
if($singletons) {               #-- Excluding singleton raw reads from binning
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


#-- Reading contigs
my @allcontigs;
my %allcontigs;

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
		if(($numgenes>=$mingenes14) && ($chimlevel<=$maxchimerism14)) { push(@allcontigs,$r[0]); $allcontigs{$r[0]}=1; }	
		if($smallnoannot14 && ($numgenes<=1) && ($r[1] eq "Unknown")) { delete $allcontigs{$r[0]}; }
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
		push(@allcontigs,$r[0]); $allcontigs{$r[0]}=1;
		}
	}



#-- Create a temp fasta excluding the singleton contigs
my $tempfasta="$tempdir/bincontigs.fasta";
open(outfile1,">$tempfasta") || die "Can't open $tempfasta for writing\n";
open(infile1,$contigsfna) || die "Can't open $contigsfna\n";
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


#-- Create a contig depth file excluding the singleton contigs
open(outfile2,">$depthfile") || die "Can't open $depthfile for writing\n";
open(infile2,$depthfileprov) || die "Can't open $depthfileprov\n";
while(<infile2>) {
	chomp;
	next if !$_;
	my @z=split(/\t/,$_);
	unless($allcontigs{$z[0]} or $z[0] eq "contigName") { next; } # Skip if the contigs is low quality (except if this is the first line)
	unless($singletonlist{$z[0]}) { print outfile2 "$_\n"; }      # Write only if this is not a singleton
}


#-- Run metabat2
my $command="$metabat_soft -t $numthreads -i $tempfasta -a $depthfile -o $dirbin/metabat2";
print outsyslog "Running metabat2 : $command\n";
print "  Running metabat2 (Kang et al 2019, PeerJ 7, e7359)\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "Binning was done using Metabat2 (Kang et al 2019, PeerJ 7, e7359)\n";
close outmet;


