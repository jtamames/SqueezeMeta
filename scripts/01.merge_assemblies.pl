#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Merges individual assemblies using minimus2, for merged mode. It also uses cd-hit for excluding identical contigs

use strict;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

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

my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

#-- Configuration variables from conf file

our($resultpath,$interdir,$tempdir,$scriptdir,$cdhit_soft,$contigid,$extassembly,$minimus2_soft,$toamos_soft,$prinseq_soft,$numthreads,$methodsfile,$syslogfile,$singletons,$force_overwrite);

#-- Merges the assemblies in a single dataset

open(out_met,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my $finalcontigs="$resultpath/01.$project.fasta";
my($ecode,$command);
if($extassembly) { 
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	system("cp $extassembly $finalcontigs");
	}
#elsif((-e $finalcontigs) && (!$force_overwrite)) { print "  Assembly results already present in file $finalcontigs, skipping\n"; exit; }
else {
	if(-e $finalcontigs) { system("rm $finalcontigs"); }
	my $merged="$tempdir/mergedassemblies.$project.fasta";
	$command="cat $interdir/01*fasta > $merged";
	print outsyslog "Merging assemblies: $command\n";
	system $command;
	if(-z $merged) { die "$merged is empty\n"; }

	#-- Uses cd-hit to identify and remove contigs contained in others

	my $merged_clustered="$tempdir/mergedassemblies.$project.99.fasta";
	$command="$cdhit_soft -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9 > /dev/null 2>&1";
	print "  Running cd-hit-est for removing redundant contigs\n";
	print outsyslog "Running cd-hit-est for removing redundant contigs: $command\n";
	$ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	if(-z $merged_clustered) { die "$merged_clustered is empty\n"; }
	print out_met "Redundant contigs were removed using cd-hit(Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";

	#-- Uses Amos to chage format to afg (for minimus2)

	my $afg_format="$tempdir/mergedassemblies.$project.99.afg";
	$command="$toamos_soft -s $merged_clustered -o $afg_format > /dev/null 2>&1 ";
	print "  Transforming to afg format\n";
	print outsyslog "Transforming to afg format: $command\n";
	$ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	if(-z $afg_format) { die "$afg_format is empty\n"; }

	#-- Uses minimus2 to assemble overlapping contigs

	$command="$minimus2_soft $tempdir/mergedassemblies.$project.99 -D OVERLAP=100 -D MINID=95 -D THREADS=$numthreads >> $syslogfile 2>&1";
	print "  Merging with minimus2\n";
	print outsyslog "Merging with minimus2: $command";
	$ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	if(-z $afg_format) { die "$afg_format is empty\n"; }
	print out_met "Contigs were merged using Minimus2(Treangen et al 2011, Curr Protoc Bioinfomatics 11)\n";

	#-- Create the final result (overlapping contigs plus singletons)

	$command="cat $tempdir/mergedassemblies.$project.99.fasta $tempdir/mergedassemblies.$project.99.singletons.seq > $finalcontigs.prov";
	print outsyslog "Joining contigs and singletons: $command\n";
	system($command);

	#-- Renaming contigs

	my $cocount;
	if(!$contigid) { $contigid="Merged"; }
	open(outfile0,">$finalcontigs") || die;
	open(infile0,"$finalcontigs.prov") || die;
	while(<infile0>) {
		chomp;
		if($_=~/^\>/) {
			$cocount++;
			my $newname="$contigid\_$cocount";
			print outfile0 ">$newname\n"; 
			}
		else { print outfile0 "$_\n"; }
		}
	close infile0;
	close outfile0;
	system("rm $finalcontigs.prov");
			

	#-- Remove files from temp

	 system("rm -r $tempdir/mergedassemblies*");
	}

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $finalcontigs -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats ";
print outsyslog "Using prinseq for statistics: $command\n";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print out_met "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	

#-- Count length of contigs (needed later)

my $contigslen="$interdir/01.$project.lon";
print "  Counting contig lengths\n";
open(outfile1,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile1,$finalcontigs) || die "Can't open $finalcontigs\n";
my($thisname,$contigname,$numc,$seq);
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
	$numc++;
	$thisname=$1;
	if($contigname) {
		my $len=length $seq;
		print outfile1 "$contigname\t$len\n"; 
		}
	$seq="";
	$contigname=$thisname;
	}
	else { $seq.=$_; }
}
close infile1;
if($contigname) { my $len=length $seq; print outfile1 "$contigname\t$len\n"; }
close outfile1;
print "  Number of contigs: $numc\n";
print outsyslog "  Number of contigs: $numc\n";

close outsyslog;
close out_met;

