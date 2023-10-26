#!/usr/bin/env perl

$|=1;

use strict;
use Cwd;
use Linux::MemInfo;
use lib ".";

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
my $sample=$ARGV[1];
my $par1name=$ARGV[2];
my $par2name=$ARGV[3];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

#-- Configuration variables from conf file

our($installpath);

# Just pass things to assembly_spades.pl, which can deal with different spades modes.
my $command="perl $installpath/lib/SqueezeMeta/assembly_spades.pl $projectdir $sample $par1name $par2name";
my $ecode=system $command;
if($ecode) { die "Error running command     $command";}
