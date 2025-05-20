#!/usr/bin/env perl

use strict;
use Cwd;
use Linux::MemInfo;
use lib ".";

$|=1;

use File::Basename;
use Cwd 'abs_path';
our $sqmlibdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $sqmlibdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $sqmlibdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$sqmlibdir/../..");

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

# Just pass things to assembly_spades.pl, which can deal with different spades modes.
my $command="perl $installpath/lib/SqueezeMeta/assembly_spades.pl $projectdir $sample $par1name $par2name";
my $ecode=system $command;
if($ecode) { die "Error running command     $command";}
