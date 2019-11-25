#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs hmmer searches against Pfam databases

use strict;
use Cwd;
use lib ".";

my $pwd=cwd();

my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

#-- Configuration variables from conf file

our($hmmer_soft,$pfamhmmer,$numthreads,$pfam_db,$aafile,$evaluehmmer5,$methodsfile,$syslogfile);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

print "Running HMMER3 (Eddy 2009, Genome Inform 23, 205-11) for Pfam\n";
my $command="$hmmer_soft --domtblout $pfamhmmer -E $evaluehmmer5 --cpu $numthreads $pfam_db $aafile > /dev/null 2>&1";
print outsyslog "Running HMMER3 for Pfam: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "HMM homology searches were done by HMMER3 (Eddy 2009, Genome Inform 23, 205-11) for the Pfam database (Finn et al 2016, Nucleic Acids Res 44, D279-85)\n";
close outmet;
close outsyslog;


