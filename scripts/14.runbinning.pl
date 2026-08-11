#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 27/01/2020 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs binning

use strict;
use Cwd;
use Term::ANSIColor qw(:constants);
use lib ".";
use Getopt::Long;

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

my $pwd=cwd();
my $projectpath = shift @ARGV;

if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

#-- Configuration variables from conf file

our($databasepath,$resultpath,$interdir,$contigsfna,%binscripts,$contigcov,$maxbin_soft,$alllog,$tempdir,$numthreads,$mappingfile,$binners,$methodsfile,$syslogfile);

#-- Handle positional args
GetOptions( 'threads=i' => \my $numthreads_override
          , 'force_overwrite' => \my $force_overwrite
          );

#-- Override numthreads if requested

if($numthreads_override) { $numthreads = $numthreads_override; }

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my @binner=split(/\,/,$binners);

foreach my $tbinner(@binner) { #-- For all the specified binners

	#-- Check whether we have results and skip in that case
	my @binfiles;
	my $wsize=0;
	my $firstfile;
	my $dirbin="$interdir/binners/$tbinner";
	opendir(indir1,$dirbin);
	@binfiles=grep(/fasta$|fa$/,readdir indir1);
	closedir indir1;
	$firstfile="$dirbin/$binfiles[0]";
	# Check that there is at least one file and that it has content
	if(scalar @binfiles) { $wsize=checksize($firstfile); }
	if(($wsize>=2) && (!$force_overwrite)) {
		print "  Binning result $firstfile already found for binner $tbinner, skipping\n";
		print outsyslog "  Binning result $firstfile already found for binner $tbinner, skipping\n";
		next; 
		}
	#-- If we are going to run the binner, remove whatever intermediate result is left from previous runs
	if(-d $dirbin) { system("rm -r $dirbin; mkdir $dirbin"); }
	
	#-- Run the binner
	
	my $scriptname=$binscripts{$tbinner};
	if(!$scriptname) {
		print RED; print "WARNING in STEP14 -> No binner found for $tbinner\n"; print RESET;
		print outsyslog "WARNING in STEP14 -> No binner found for $tbinner\n";
		next;
		}
	print "  Running $tbinner from $scriptname\n";
	my $ecode = system("LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installpath/lib perl $scriptname $projectpath -t $numthreads >> $tempdir/$project.log");
	if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
	
	#-- Check the bins, to verify that all is correct (there are at least some bins)
	opendir(indir1,$dirbin) || die "Can't open $dirbin directory\n";
	@binfiles=grep(/fasta$|fa$/,readdir indir1);
	closedir indir1;
	$firstfile="$dirbin/$binfiles[0]";
	if(-e $firstfile) { $wsize=checksize($firstfile); }
	else { $wsize=0; }
	if($wsize<2) {
		print RED; print "WARNING in STEP14 -> $scriptname. No $tbinner results!\n"; print RESET;
		print outsyslog "WARNING in STEP14 -> $scriptname. No $tbinner results!\n"; 
		}
	}

close outsyslog;

sub checksize {
        my $tfile=shift;
        my $wsize;
        if(-e $tfile) {
                $wsize=qx(grep -cv "^#" $tfile); # this excludes comments!
                }
        else { $wsize=0; }
        return $wsize;
        }

