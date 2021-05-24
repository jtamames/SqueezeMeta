#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Diamond for homology searches, against nr, COGs and KEGG databases

use strict;
use Cwd;
use Linux::MemInfo;
use lib ".";

$|=1;

my $pwd=cwd();

my $projectdir=$ARGV[0];
my $notax=$ARGV[1];
my $blastmode=$ARGV[2];

if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($aafile,$databasepath,$numthreads,$diamond_soft,$nodiamond,$nocog,$nokegg,$interdir,$cog_db,$kegg_db,$nr_db,$blocksize,$evaluetax4,$minidentax4,$evaluefun4,$minidenfun4,$cogdiamond,$keggdiamond,$taxdiamond,$opt_db,$resultpath,$methodsfile,$syslogfile);
my $command;

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Setting block size for Diamond
	
if($blocksize eq "NF") {
	print "  Setting block size for Diamond\n";
	my %mem=get_mem_info;
	my $ram=$mem{"MemAvailable"}/(1024*1024);
	my $ramstr=sprintf('%.2f',$ram);
	my $block_size_set=sprintf('%.1f',$ram/5);
	if($block_size_set>8) { $block_size_set=8; }	
	if($block_size_set<1) { $block_size_set=1; }
	print "  AVAILABLE (free) RAM memory: $ramstr Gb\n  We will set Diamond block size to $block_size_set (Gb RAM/5, Max 8).\n  You can override this setting using the -b option when starting the project, or changing\n  the \$blocksize variable in SqueezeMeta_conf.pl\n";
	print outsyslog "Diamond block size set to $block_size_set (Free Mem $ramstr Gb)\n";
	$blocksize=$block_size_set;
	}

my($taxfound,$cogfound,$keggfound);
if(-e $taxdiamond) { $taxfound=1; }
if(-e $cogdiamond) { $cogfound=1; }
if(-e $keggdiamond) { $keggfound=1; }
my($donediamond,$stringmethods);
if(!$blastmode) { $blastmode="blastp"; }
print outmet "Similarity searches for ";
$stringmethods="Similarity searches for ";

$command="cp $databasepath/DB_BUILD_DATE $interdir";
my $ecode = system $command;
if($ecode!=0) { warn "Error running command:     $command"; }

#-- nr database

if(!$notax) {
	$command="$diamond_soft $blastmode -q $aafile -p $numthreads -d $nr_db -e $evaluetax4 --id $minidentax4 -f tab -b $blocksize --quiet -o $taxdiamond";
	print " taxa";
	print outsyslog "Running Diamond for taxa: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "GenBank (Clark et al 2016, Nucleic Acids Res 44, D67-D72), ";
	}

#-- COG database

if(!$nocog) {
	$command="$diamond_soft $blastmode -q $aafile -p $numthreads -d $cog_db -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $cogdiamond";
	print " COGS";
	print outsyslog "Running Diamond for COGs: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "eggNOG (Huerta-Cepas et al 2016, Nucleic Acids Res 44, D286-93), ";
}

#-- KEGG database

if(!$nokegg) {
	if((!$nodiamond) || ($nodiamond && !$keggfound)) {
		$command="$diamond_soft blastp -q $aafile -p $numthreads -d $kegg_db -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $keggdiamond";
		print "   Running Diamond (Buchfink et al 2015, Nat Methods 12, 59-60) for KEGG\n";
		print outsyslog "Running Diamond for KEGG: $command\n";
		my $ecode = system $command;
		$donediamond=1;
		if($ecode!=0) { die "Error running command:    $command"; }
		$stringmethods.="KEGG (Kanehisa and Goto 2000, Nucleic Acids Res 28, 27-30), ";
		}
	else { print "  Found --nodiamond flag and Diamond result for KEGG ($keggdiamond): skipping\n"; }	
	}

#-- Optional databases

my $outfound;
if($opt_db) {
	open(infile1,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($dbname,$extdb,$dblist)=split(/\t/,$_);
		my $outdb="$interdir/04.$project.$dbname.diamond";
		$command="$diamond_soft $blastmode -q $aafile -p $numthreads -d $extdb -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outdb";
		print " $dbname";
		print outsyslog "Running Diamond for $dbname: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		print outmet "$dbname, ";
		}
}

$stringmethods.=" were done using Diamond (Buchfink et al 2015, Nat Methods 12, 59-60)\n";
print "\n";
if($donediamond) { print outmet $stringmethods; } else { print outsyslog "Skipping Diamond runs because of --nodiamond flag\n"; }
close outsyslog;
close outmet;

