#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs hmmer searches against Pfam databases

use strict;
use warnings;
use Cwd;

my $pwd=cwd();

my $project=$ARGV[0];
$project=~s/\/$//; 

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($hmmer_soft,$pfamhmmer,$numthreads,$pfam_db,$aafile);

print "Running hmmer search\n";
my $command="$hmmer_soft --domtblout $pfamhmmer -E 1e-10 --cpu $numthreads $pfam_db $aafile > /dev/null";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

