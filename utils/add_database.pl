#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Diamond for homology searches, against nr, COGs and KEGG databases

use strict;
use Cwd;
use Linux::MemInfo;
use Term::ANSIColor qw(:constants);
use Tie::IxHash;
use lib ".";

$|=1;

my $pwd=cwd();

my $projectdir=$ARGV[0];
my $new_opt_db=$ARGV[1];

my $helpshort="Usage: add_database.pl <project> <database file>\n";

my $helptext = <<END_MESSAGE;
Usage: add_database.pl <project> <database file>

Mandatory parameters:
   project : Project name (a project must exist already)
   database file: File containing the information of the databases to add. Please refer to the manual for details.
	
END_MESSAGE

if(!$projectdir) { die "$helptext\nPlease provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "$helpshort\nCan't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($diamond_soft);
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($aafile,$databasepath,$numthreads,$diamond_soft,$nodiamond,$nocog,$nokegg,$nopfam,$opt_db,$scriptdir,$interdir,$tempdir,$cog_db,$kegg_db,$nr_db,$blocksize,$evaluetax4,$minidentax4,$evaluefun4,$minidenfun4,$cogdiamond,$keggdiamond,$taxdiamond,$resultpath,$methodsfile,$syslogfile);
my $command;
open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


print BOLD "\nadd_database.pl  - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;

my $blastmode="blastp";
my $noprev;
if($opt_db) { $noprev=1; }
if($new_opt_db) { $opt_db=$new_opt_db; }
if(!$opt_db) { die "$helpshort\nList of external databases needs to be specified. Please see the manual for knowing the exact format of that file\n"; }

#-- Running Diamond on newc optional databases

my $outfound;
open(infile1,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my($dbname,$extdb,$dblist)=split(/\t/,$_);
	my $outdb="$interdir/04.$project.$dbname.diamond";
	if(-e $outdb) { print "  File $outdb for database $dbname found, skipping run\n"; print outsyslog "  File $outdb for database $dbname found, skipping run\n"; }
	else {
		$command="$diamond_soft $blastmode -q $aafile -p $numthreads -d $extdb -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outdb";
		print " $dbname";
		print outsyslog "Running Diamond for $dbname: $command\n";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		}
	}

print "\n";
close outsyslog;
close outmet;

#-- Change SqueezeMeta_conf.pl for adding ext_db (if not present already)

if($new_opt_db && (!$noprev)) {
	my $oldconf="$projectdir/SqueezeMeta_conf.pl";
	my $provline="$tempdir/provline.txt";
	my $provconf="$tempdir/provconf.txt";
	open(out,">$provline") || die "Cannot open provisional conf file in $provline\n";
	print out "\n#-- New external database\n\n\$opt_db        = \"$new_opt_db\";\n";
	close out;
	system("cat $oldconf $provline > $provconf");
	system("mv $provconf $oldconf");
	system("rm $provline");
	}

#-- Calling other scripts for creating and updating annotations	
	
system("$scriptdir/07.fun3assign.pl $projectdir");
system("$scriptdir/12.funcover.pl $projectdir");
system("$scriptdir/13.mergeannot2.pl $projectdir");
system("$scriptdir/21.stats.pl $projectdir");
		
		
		
	
	

