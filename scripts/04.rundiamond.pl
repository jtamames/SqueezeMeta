#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Diamond for homology searches, against nr, COGs and KEGG databases

use strict;
use warnings;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

#-- Configuration variables from conf file

our($aafile,$numthreads,$diamond_soft,$nocog,$nokegg,$interdir,$cog_db,$kegg_db,$nr_db,$blocksize,$evaluetax4,$minidentax4,$evaluefun4,$minidenfun4,$cogdiamond,$keggdiamond,$taxdiamond,$opt_db,$resultpath);
my $command;

#-- COG database

if(!$nocog) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $cog_db -e $evaluefun4 --id $minidenfun4 --quiet -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $cogdiamond";
	print "Running Diamond for COGS (This can take a while, please be patient)\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
}

#-- KEGG database

if(!$nokegg) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $kegg_db -e $evaluefun4 --id $minidenfun4 --quiet -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $keggdiamond";
	print "Running Diamond for KEGG (This can take a while, please be patient)\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
}

#-- Optional databases

if($opt_db) {
	open(infile1,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($dbname,$extdb,$dblist)=split(/\t/,$_);
		my $outdb="$interdir/04.$project.$dbname.diamond";
		$command="$diamond_soft blastp -q $aafile -p $numthreads -d $extdb -e $evaluefun4 --id $minidenfun4 --quiet -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outdb";
		print "Running Diamond for optional database $dbname\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		}
}

#-- nr database

$command="$diamond_soft blastp -q $aafile -p $numthreads -d $nr_db -e $evaluetax4 --id $minidentax4 -f tab -b $blocksize --quiet -o $taxdiamond";
print "Running Diamond for taxa (This can take a long while, please be even more patient)\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
