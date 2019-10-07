#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs hmmer searches against Pfam databases

use strict;
use warnings;
use Cwd;
use lib ".";

my $pwd=cwd();

my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

#-- Configuration variables from conf file

our($hmmer_soft,$pfamhmmer,$numthreads,$pfam_db,$aafile,$evaluehmmer5);

print "Running hmmer search\n";
my $command="$hmmer_soft --domtblout $pfamhmmer -E $evaluehmmer5 --cpu $numthreads $pfam_db $aafile > /dev/null";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

