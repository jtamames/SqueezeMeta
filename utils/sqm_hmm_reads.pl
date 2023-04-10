#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use lib ".";
use strict;
use Term::ANSIColor qw(:constants);


###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my $pwd=cwd();
my $utilsdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $utilsdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$utilsdir/..");
my $scriptdir = "$installpath/scripts";
my $auxdir = "$installpath/lib/SQM_reads";
my $shortpair_soft = abs_path("$utilsdir/../bin/Short-Pair/Short-Pair.py");

###

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

my($numthreads,$pair1,$pair2,$pfam,$outfile,$dietext,$tempfile1,$tempfile2,$hel,$printversion);

my $start_run = time();

my $helptext = <<END_MESSAGE;
Search for specific PFAMs in a collection of short reads using the Short-Pair tool (PMID: 29072140).

Usage: sqm_hmm_reads.pl -pfam <PFAM list> -pair1 <pair1 fasta file> -pair2 <pair2 fasta file> [options]

Arguments:

 Mandatory parameters:
   -pfam: List of PFAM ids to retrieve, comma-separated (eg: PF00069,PF00070) (REQUIRED)
   -pair1: Fasta file for pair1 (REQUIRED)
   -pair2: Fasta file for pair2 (REQUIRED)
   
 Optional:
   -t: Number of threads (Default: 12)
   -output: Name of the output file (Default: sqm_pfam.out)
   -v|version: Print version
   -h: this help


END_MESSAGE

my $result = GetOptions ("t=i" => \$numthreads,
			"pair1=s" => \$pair1,
			"pair2=s" => \$pair2,
			"pfam=s" => \$pfam,
			"output=s" => \$outfile,
		        "v|version" => \$printversion,
		    	"h" => \$hel
			);
			

print "\nsqm_hmm_reads.pl v$version- (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n
The Short-Pair tool is published in Techa-Angkoon, Sun & Lei J, BMC bioinformatics 18:414 (2017). doi: https://doi.org/10.1186/s12859-017-1826-2\n\n";
if($printversion) { exit; }

if($hel) { die "$helptext\n"; } 
if(!$pair1) { $dietext.="MISSING PAIR1\n"; }
if(!$pair2) { $dietext.="MISSING PAIR2\n"; }
if(!$pfam) { $dietext.="MISSING PFAM\n"; }
if($dietext) { die "$dietext\n$helptext\n"; }
if($pair1!~/\.fasta\b|\.fa\b/) { print RED; print "WARNING: This script requires FASTA files, and yours are not named that way\n"; print RESET; print "I will proceed assuming these are fasta files, but if I crash it will likely because of this (and you should feel very bad)\n\n"; }

if(!$outfile) { $outfile="sqm_pfam.out"; }
my $pfamfile="pfam.hmm";
my $pfamseed="pfam.seed";

get_pfam();
rewrite_files();
run_short_pairs();

sub run_short_pairs {
	system("mkdir out1");
	print "\nRunning ShortPair for PFAM $pfam\n  (This can take a while, please be patient)\n";
	my $command="$shortpair_soft -m $pfamfile -s $pfamseed -x $tempfile1 -y $tempfile2 -o $outfile";
	system($command);
	system("rm *alldomains.allframe*; rm fragment_length*; rm hmms.sav; rm pfam.seed.sav; rm $pfamfile; rm $pfamseed");
	system("rm $tempfile1; rm $tempfile2");
	system("rm -r Protein; rm -r Out_extracted; rm -r out1; rm -r HMMs; rm -r faaSP; rm -r fastaSP");
	my $currtime=timediff();
	print "\n[",$currtime->pretty,"]: Output created in $outfile\n";
	}

sub rewrite_files {
	my($hcount,$acount,$swc);
	$tempfile1=$pair1."temp.1.fasta";
	$tempfile2=$pair2."temp.2.fasta";
	if($pair1=~/gz$/) { open(infile1,"zcat $pair1 |") || die "Cannot open input file 1 in $pair1\n"; }
	else { open(infile1,$pair1) || die "Cannot open input file 1 in $pair1\n"; }
	while(<infile1>) {		#-- To know if the naming schema has .1, \1, _1 at the end of read names
		chomp;			#-- We evaluate the first 10 headers to check if all follow the same pattern
		next if(!$_);
		if($_=~/^\>([^ ]+)/) { 
			$hcount++;
			my $nameseq=$1;
			if($nameseq=~/\W1$/) { $acount++; }
			last if($hcount>=10);
			}
		}
	close infile1;
	if($acount==$hcount) { $swc=1; }	#-- If they do, we will remove that .1, \1, _1 and replace it for .1, which ShortPair wants
	if($pair1=~/gz$/) { open(infile1,"zcat $pair1 |") || die "Cannot open input file 1 in $pair1\n"; }
	else { open(infile1,$pair1) || die "Cannot open input file 1 in $pair1\n"; }
	open(outfile1,">$tempfile1") || die "Cannot open pair1 OUTPUT file in $tempfile1, check if you have permissions in that directory\n";
	while(<infile1>) {
		chomp;
		next if(!$_);
		if($_=~/^\>([^ ]+)/) { 
			my $nameseq=$1;
			if($swc) { $nameseq=~s/\W1$/\.1/; } else { $nameseq.=".1"; }	#-- Substitution occurs here
			print outfile1 ">$nameseq\n"; 
			}
		else { print outfile1 "$_\n"; }
		}
	close infile1;
	close outfile1;
	if($pair2=~/gz$/) { open(infile2,"zcat $pair2 |") || die "Cannot open input file 2 in $pair2\n"; }
	else { open(infile2,$pair2) || die "Cannot open input file 2 in $pair2\n"; }
	open(outfile2,">$tempfile2") || die "Cannot open pair2 OUTPUT file in $tempfile2, check if you have permissions in that directory\n";
	while(<infile2>) {
		chomp;
		next if(!$_);
		if($_=~/^\>([^ ]+)/) { 
			my $nameseq=$1;
			if($swc) { $nameseq=~s/\W2$/\.2/; } else { $nameseq.=".2"; }	#-- Same for the second pair
			print outfile2 ">$nameseq\n"; 
			}
		else { print outfile2 "$_\n"; }
		}
	close infile2;
	close outfile2;
	}	


sub get_pfam {		#-- Retrieving corresponding pfam from PFAM
	my @plist=split(/\,/,$pfam);
	foreach my $tpfam(@plist) {
		print "Getting Pfam and alignment for $tpfam\n";
		
		#-- WARNING: Pfam has been DISCONTINUED and it is now included in UniProt, which means that pfam.xfam.org no longer exists.
		#-- for now, the pfam-legacy.xfam.org URL still works, but can cause trouble in the future.
		
		my $hmm_query="wget -nv http://pfam-legacy.xfam.org/family/$tpfam/hmm -O $tpfam.hmm";
		system $hmm_query;
		if(-e $pfamfile) { system("cat $pfamfile $tpfam.hmm > outh; mv outh $pfamfile; rm $tpfam.hmm"); }
		else { system("mv $tpfam.hmm $pfamfile"); }
		my $ali_query="wget -nv http://pfam-legacy.xfam.org/family/$tpfam/alignment/seed -O $tpfam.seed";
		system $ali_query;
		if(-e $pfamseed) { system("cat $pfamseed $tpfam.seed > outh; mv outh $pfamseed; rm $tpfam.seed"); }
		else { system("mv $tpfam.seed $pfamseed"); }
		}
	}
	

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}
