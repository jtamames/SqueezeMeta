#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 30/01/2020 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Concoct for binning

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

our($databasepath,$contigsfna,$contigcov,$alllog,$tempdir,$interdir,$datapath,$numthreads,$mappingfile,$methodsfile,$syslogfile);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my $concoctdir="/home/jtamames/software/CONCOCT-1.1.0";
my $bindir="$interdir/binners/concoct";
print outsyslog "\nRUNNING concoct\n";
if(-d $bindir) {} else { print "  Creating $bindir directory\n"; print outsyslog "  Creating $bindir directory\n"; system("mkdir $bindir"); }

my $samdir="$datapath/sam";
opendir(indir,$samdir) || die;
my @samfiles=grep(/sam$/,readdir indir);
closedir indir;

my $command;
foreach my $sam(@samfiles) {
	print "  Working for $sam\n";
	print outsyslog "  Working for $sam\n";
	my $samfile="$samdir/$sam";
	my $bamfile="$tempdir/$sam";
	$bamfile=~s/sam$/bam/;
	$command="samtools view -b -\@$numthreads -S $samfile -o $bamfile > /dev/null 2>&1";
	print outsyslog "  Creating BAM for $sam: $command\n";
	print "  Creating BAM for $sam\n";
	system $command;
	my $sorted=$bamfile;
	$sorted=~s/\.bam$/\.sorted/;
	$command="samtools sort $bamfile $sorted -\@$numthreads > /dev/null 2>&1";
	print "  Sorting BAM\n";
	print outsyslog "  Sorting BAM: $command\n";
	system $command;
	$command="samtools index $sorted.bam > /dev/null 2>&1";
	print outsyslog "  Indexing sorted BAM: $command\n";
	print "  Indexing sorted BAM\n";
	system $command;
	}
	
print "\n  Cutting contigs in pieces!\n";
my $bedfile="$tempdir/$project.contigs.bed";
my $choppedcontigs="$tempdir/$project.choppedcontigs.fasta";
$command="$concoctdir/scripts/cut_up_fasta.py $contigsfna -c 10000 -o 0 --merge_last -b $bedfile > $choppedcontigs";
print outsyslog "\n  Cutting contigs in pieces!: $command\n";	
system $command;	
print "  Creating abundance table\n";	
$command="$concoctdir/scripts/concoct_coverage_table.py $bedfile $tempdir/*.sorted.bam > $bindir/coverage_table.tsv";
print outsyslog "  Creating abundance table: $command\n";
system($command);
print "  Running concoct\n";
$command="$concoctdir/bin/concoct --composition_file $choppedcontigs --coverage_file $bindir/coverage_table.tsv --threads $numthreads  -b $bindir/concoct_int/ > /dev/null 2>&1";
print outsyslog "  Running concoct: $command\n";
system($command);
print "  Merging clusters\n";
$command="$concoctdir/scripts/merge_cutup_clustering.py $bindir/concoct_int/clustering_gt1000.csv > $bindir/concoct_int/clustering_merged.csv";
print outsyslog "  Merging clusters: $command\n";
system($command);
print "  Extracting final bins\n";
$command="$concoctdir/scripts/extract_fasta_bins.py $contigsfna $bindir/concoct_int/clustering_merged.csv --output_path $bindir > /dev/null 2>&1";
print outsyslog "  Extracting final bins: $command\n";
system($command);
print "  Renaming bins\n";
opendir(indir1,$bindir) || die "Cannot open bin directory $bindir\n";
my @mfiles=grep(/\.fa/,readdir indir1);
closedir(indir1);
foreach my $renam(@mfiles) {
	my $newname="concoct.".$renam;
	system("mv $bindir/$renam $bindir/$newname");
	}

close outsyslog;
