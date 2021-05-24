#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use Linux::MemInfo;
use Term::ANSIColor qw(:constants);
use lib ".";
use strict;
use threads;

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my($symlinkpath,$symlinkdest,$utilsdir);
if(-l __FILE__)
        {
        $symlinkpath = dirname(__FILE__);
        $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $utilsdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$utilsdir/..");
my $scriptdir = "$installpath/scripts";
my $auxdir = "$installpath/lib/SQM_reads";
###

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

my($numthreads,$outdir,$aadir,$project,$samplesfile,$notax,$nocog,$nokegg,$blastmode,$blocksize,$hel,$blockoption,$printversion);

my $start_run = time();
print BOLD "\nSQM_annot v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;

my $helpshort="Usage: SQM_mapper.pl -r <reference file> -s <samples file> -f <raw fastq dir> -g <GFF file> [options]\n";

my $helptext = <<END_MESSAGE;
Usage: SQM_annot.pl -f <faa directory> -s <samples file> -o <output directory> [options]

Mandatory parameters:
   -f: Directory for sequence files (REQUIRED)
   -o: Output directory (REQUIRED)
   -s: Samples file (REQUIRED)

 Options:
   -t: Number of threads (Default: 12)
   -b: Diamond block size
   -blastmode: blastp (for aminoacid sequences) or blastx (for nucleotide) (Default: blastp)
   --notax: Skip taxonomic annotation
   --nocog: Skip COGs annotation
   --nofun: Skip KEGG annotation
   -version: Print version
   -h: this help
   
END_MESSAGE

my $result = GetOptions ("t=i" => \$numthreads,
                     "s=s" => \$samplesfile,
                     "f=s" => \$aadir,
		     "o=s" => \$outdir,
		     "b=i" => \$blocksize,
		     "notax" => \$notax,
		     "nocog" => \$nocog,
		     "nokegg" => \$nokegg,
		     "blastmode=s" => \$blastmode,
		     "v|version" => \$printversion,
		     "h" => \$hel
		    );

if($hel) { print "$helptext\n"; exit; }
if($printversion) { exit; }
if(!$samplesfile) { die "Please specify a file of samples\n"; }
if(!$outdir) { $outdir="."; }
if(!$aadir) { $aadir="."; }
if($blocksize) { $blockoption="-b $blocksize"; }
if($blastmode=~/blastx/) { $blastmode="blastx"; } else { $blastmode="blastp"; }

my %files;
tie %files,"Tie::IxHash";
open(in,$samplesfile) || die;
my $tempsample="temp.samples";  #-- We need a samples file to trick SqueezeMeta
open(out,">$tempsample") || die;
while(<in>) {
	chomp;
	my @e=split(/\t/,$_);
	$files{$e[0]}=$e[1];
	print out "$_\tpair1\n";
	}
close out;

my $command="perl $scriptdir/SqueezeMeta.pl  -s temp.samples -f $aadir -m sequential $blockoption --nopfam --empty";
system($command);
foreach my $thissample(keys %files) {
	print "Working with $thissample\n";
	my $aafile=$files{$thissample};
	$project=$thissample;
	system("cp $aadir/$aafile $project/results/03.$project.faa");
	$command="perl $scriptdir/04.rundiamond.pl $project $notax $blastmode";
	system $command;

	if(!$notax) { $command="perl $scriptdir/06.lca.pl $project"; system $command; }
	if((!$nocog) || (!$nokegg)) { $command="perl $scriptdir/07.fun3assign.pl $project"; system $command; }
	}
	
