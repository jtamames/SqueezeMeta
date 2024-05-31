#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 27/01/2020 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs MaxBin for binning

use strict;
use Cwd;
use Term::ANSIColor qw(:constants);
use lib ".";

$|=1;

my $pwd=cwd();
my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

#-- Configuration variables from conf file

our($installpath,$databasepath,$installpath,$resultpath,$interdir,$contigsfna,%binscripts,$contigcov,$maxbin_soft,$alllog,$tempdir,$numthreads,$mappingfile,$binners,$methodsfile,$syslogfile);

my $finaltrace;
my @binner=split(/\,/,$binners);

foreach my $tbinner(@binner) { #-- For all the specified binners
	my @binfiles;
	my $wsize=0;
	my $firstfile;
	my $dirbin="$interdir/binners/$tbinner";
	if(-d $dirbin) {	#-- If the result directory exists, don't create it, and check if bins are already there
		opendir(indir1,$dirbin) || die "Can't open $dirbin directory\n";
		@binfiles=grep(/fasta$|fa$/,readdir indir1);
		closedir indir1;
		$firstfile="$dirbin/$binfiles[0]";
		# $wsize=checksize($firstfile); #commented out since this was no longer used
		}
	else { system("mkdir $dirbin"); }
	
	#-- Skip the binning if results are already present
		
       	# if($wsize>2)         { print "Binning result $firstfile already found for binner $tbinner, skipping\n"; next; }		
	
	#-- Run the binner
	
	my $scriptname=$binscripts{$tbinner};
	if(!$scriptname) { print RED; print "WARNING in STEP14 -> No binner found for $tbinner\n"; print RESET; $finaltrace.="WARNING in STEP15: No binner found for $tbinner\n"; next; }
	print "  Running $tbinner from $scriptname\n";
	my $ecode = system("LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installpath/lib perl $scriptname $projectpath >> $tempdir/$project.log");
	if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
	
	#-- Check the bins, to verify that all is correct (there are at least some bins)
	opendir(indir1,$dirbin) || die "Can't open $dirbin directory\n";
	@binfiles=grep(/fasta$|fa$/,readdir indir1);
	closedir indir1;
	$firstfile="$dirbin/$binfiles[0]";
	if(-e $firstfile) { $wsize=checksize($firstfile); }
	else { $wsize=0; }
	if($wsize<2) { print RED; print "WARNING in STEP14 -> $scriptname. No $tbinner results!\n"; print RESET; $finaltrace.="WARNING in STEP14: No $tbinner results!\n"; }
}


sub checksize {
	my $tfile=shift;
	my $wc=qx(wc -l $tfile);
        my($wsize,$rest)=split(/\s+/,$wc);
	return $wsize;
	}
