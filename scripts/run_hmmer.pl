#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs hmmer searches against Pfam databases

use strict;
use warnings;
use Cwd;

my $pwd=cwd();

my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($hmmer_soft,$pfamhmmer,$numthreads,$pfam_db,$aafile);

print "Running hmmer search\n";
my $command="$hmmer_soft --domtblout $pfamhmmer -E 1e-10 --cpu $numthreads $pfam_db $aafile > /dev/null";
system $command;

