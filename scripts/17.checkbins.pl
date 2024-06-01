#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs checkM for evaluating bins
#

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

our($installpath,$datapath,$taxlist,$binresultsdir,$checkm2_soft,$alllog,$resultpath,$tempdir,$minsize17,$numthreads,$interdir,$methodsfile,$syslogfile,$checkmfile,$gtdbtk,$gtdbtk_data_path,$gtdbtkfile);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

print "  Evaluating bins with CheckM2 (Chklovski et al 2023, Nat Met 20, 1203-12)\n\n";

my $markerdir="$datapath/checkm_markers";
my $checktemp="$interdir/checkm2";
my $gtdbtktemp="$interdir/gtdbtk";
my $command;

my $binmethod="DAS";

if (-d $checktemp) { system "rm -r $checktemp"; }

$command = "$checkm2_soft predict -i $binresultsdir -o $checktemp -t $numthreads -x fa >> $syslogfile 2>&1";
print outsyslog "$command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

print "\n  Storing results for $binmethod in $checkmfile\n";
$command = "cp $checktemp/quality_report.tsv $checkmfile";
print outsyslog "$command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "Bin statistics were computed using CheckM2 (Chklovski et al, Nat Met 20, 1203-12)\n";

if($gtdbtk) {
	print "\n  Running GTDB-Tk to classify the bins\n";
	my $command = "GTDBTK_DATA_PATH=$gtdbtk_data_path gtdbtk classify_wf --genome_dir $binresultsdir --out_dir $gtdbtktemp -x fa --cpus $numthreads --mash_db $gtdbtk_data_path >> $syslogfile 2>&1";
	print outsyslog "$command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	my @files = glob( $gtdbtktemp . '/*summary.tsv' );
	my $nfiles = scalar(@files);
	if(!$nfiles) { die "No GTDB-Tk results found"; }
	if($nfiles>2) { die "More than two GTDB-Tk results found (Bacteria, Archaea, ??)"; }
	my $nf = 0;
	open(my $outfile, ">", $gtdbtkfile) || die "Can't open $gtdbtkfile for writing";
	foreach my $f (@files) {
		open(my $infile, "<", $f) || die "Can't open $f for reading";
		my $nl = 0;
        	while(<$infile>) {
			if($nf == 0 or $nl > 0) { print $outfile $_; }
			$nl++;
			}
		close $infile;
 		$nf++;
		}
	close $outfile;
	print "\n  GTDB-Tk results can be found in $gtdbtkfile\n";
	print outmet "Bins were classified using GTDB-Tk v2 (Chaumeil et al 2022, Bioinformatics 38, 5315-16)\n";
}

close outmet;
close outsyslog;

