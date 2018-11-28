#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Prodigal software for predicting ORFs

use strict;
use warnings;
use Cwd;

my $pwd=cwd();

my $project=$ARGV[0];
$project=~s/\/$//; 

do "$project/SqueezeMeta_conf.pl";

our($resultpath,$tempdir,$aafile,$ntfile,$gff_file,$prodigal_soft);

#-- Runs prodigal and cat the gff file with the RNA's one coming from barrnap (previous step)

my $tempgff="$tempdir/02.$project.cds.gff";
my $maskedcontigs="$resultpath/02.$project.maskedrna.fasta";
my $command="$prodigal_soft -q -m -p meta -i $maskedcontigs -a $aafile -d $ntfile -f gff -o $tempgff";
print "Running prodigal: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
system("cat $tempgff $tempdir/02.$project.rna.gff > $gff_file");
