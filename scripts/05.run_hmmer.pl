#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs hmmer searches against Pfam databases

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

our($hmmer_soft,$pfamhmmer,$numthreads,$pfam_db,$aafile,$evaluehmmer5,$methodsfile,$syslogfile);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

print "  Running HMMER3 (Eddy 2009, Genome Inform 23, 205-11) for Pfam\n";
my $command="$hmmer_soft --domtblout $pfamhmmer -E $evaluehmmer5 --cpu $numthreads $pfam_db $aafile > /dev/null 2>&1";
print outsyslog "Running HMMER3 for Pfam: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "HMM homology searches were done by HMMER3 (Eddy 2009, Genome Inform 23, 205-11) for the Pfam database (Finn et al 2016, Nucleic Acids Res 44, D279-85)\n";
close outmet;
close outsyslog;


