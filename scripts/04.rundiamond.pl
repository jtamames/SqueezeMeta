#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Diamond for homology searches, against nr, COGs and KEGG databases

use strict;
use warnings;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($aafile,$numthreads,$diamond_soft,$nocog,$nokegg,$cog_db,$kegg_db,$nr_db,$evalue,$miniden,$cogdiamond,$keggdiamond,$taxdiamond);
my $command;

#-- COG database

if(!$nocog) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $cog_db -e $evalue --id $miniden -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $cogdiamond";
	print "Running Diamond for COGS: $command\n";
	system $command;
}

#-- KEGG database

if(!$nokegg) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $kegg_db -e $evalue --id $miniden -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $keggdiamond";
	print "Running Diamond for KEGG: $command\n";
	system $command;
}

#-- nr database

$command="$diamond_soft blastp -q $aafile -p $numthreads -d $nr_db -e $evalue -f tab -b 8 -o $taxdiamond";
print "Running Diamond for taxa: $command\n";
system $command;
