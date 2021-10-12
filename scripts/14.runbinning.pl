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

foreach my $tbinner(@binner) { 
	my $scriptname=$binscripts{$tbinner};
	if(!$scriptname) { print RED; print "WARNING in STEP14 -> No binner found for $tbinner\n"; print RESET; $finaltrace.="WARNING in STEP15: No binner found for $tbinner\n"; next; }
	print "  Running $tbinner from $scriptname\n";
	my $ecode = system("perl $scriptname $projectpath >> $tempdir/$project.log");
	if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
	my $dirbin="$interdir/binners/$tbinner";
	if(-d $dirbin) {} else { system("mkdir $dirbin"); }
	my @binfiles;
	opendir(indir1,$dirbin) || die "Can't open $dirbin directory\n";
	@binfiles=grep(/fasta$|fa$/,readdir indir1);
	closedir indir1;
	my $firstfile="$dirbin/$binfiles[0]";
	my ($wsize,$rest);
	if(-e $firstfile) {
		my $wc=qx(wc -l $firstfile);
		($wsize,$rest)=split(/\s+/,$wc);
		}
	else { $wsize==0; }
	if($wsize<2) { print RED; print "WARNING in STEP14 -> $scriptname. No $tbinner results!\n"; print RESET; $finaltrace.="WARNING in STEP14: No $tbinner results!\n"; }
}

