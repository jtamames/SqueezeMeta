#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Diamond for homology searches, against nr, COGs and KEGG databases

use strict;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();

my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

#-- Configuration variables from conf file

our($aafile,$numthreads,$diamond_soft,$nocog,$nokegg,$interdir,$cog_db,$kegg_db,$nr_db,$blocksize,$evaluetax4,$minidentax4,$evaluefun4,$minidenfun4,$cogdiamond,$keggdiamond,$taxdiamond,$opt_db,$resultpath,$methodsfile,$syslogfile);
my $command;

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";
print outmet "Similarity searches for ";

print "Running Diamond (Buchfink et al 2015, Nat Methods 12, 59-60) for";
#-- COG database

if(!$nocog) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $cog_db -e $evaluefun4 --id $minidenfun4 --quiet -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $cogdiamond";
	print " COGS";
	print outsyslog "Running Diamond for COGs: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "eggNOG (Huerta-Cepas et al 2016, Nucleic Acids Res 44, D286-93), ";
}

#-- KEGG database

if(!$nokegg) {
	$command="$diamond_soft blastp -q $aafile -p $numthreads -d $kegg_db -e $evaluefun4 --id $minidenfun4 --quiet -b 8 -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $keggdiamond";
	print " KEGG";
	print outsyslog "Running Diamond for KEGG: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "KEGG (Kanehisa and Goto 2000, Nucleic Acids Res 28, 27-30), ";
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
		print " $dbname";
		print outsyslog "Running Diamond for $dbname: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		}
}

#-- nr database

$command="$diamond_soft blastp -q $aafile -p $numthreads -d $nr_db -e $evaluetax4 --id $minidentax4 -f tab -b $blocksize --quiet -o $taxdiamond";
print " taxa\n";
print outsyslog "Running Diamond for taxa: $command\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "GenBank (Clark et al 2016, Nucleic Acids Res 44, D67-D72), ";
print outmet "were done using Diamond (Buchfink et al 2015, Nat Methods 12, 59-60)\n";
close outmet;
close outsyslog;

