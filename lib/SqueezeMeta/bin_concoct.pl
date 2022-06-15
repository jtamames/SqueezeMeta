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

our($installpath, $samtools_soft, $concoct_dir,$databasepath,$contigsfna,$singletons,$contigcov,$alllog,$tempdir,$interdir,$mappingfile,$datapath,$numthreads,$mappingfile,$methodsfile,$syslogfile);

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my $bindir="$interdir/binners/concoct";
print outsyslog "\nRUNNING concoct\n";
if(-d $bindir) {} else { print "  Creating $bindir directory\n"; print outsyslog "  Creating $bindir directory\n"; system("mkdir $bindir"); }

my %skip;
print "  Reading samples from $mappingfile\n";   #-- We will exclude samples with the "nobinning" flag
open(infile0,$mappingfile) || die "Can't open $alllog\n";
while(<infile0>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	if($_=~/nobinning/) { $skip{$t[0]}=1; }
	}
close infile0;

my %singletonlist;
my $outcontig="$tempdir/contigs.nosingle.fasta";
if($singletons) {		#-- Excluding singleton raw reads from binning
	my $singletonlist="$interdir/01.$projectname.singletons";
	print "  Excluding singleton reads from $singletonlist\n";
	print outsyslog "  Excluding singleton reads from $singletonlist\n";
	open(infile0,$singletonlist) || die "Cannot open singleton list in $singletonlist\n";
	while(<infile0>) {
		chomp;
		next if !$_;
		my @y=split(/\t/,$_);
		$singletonlist{$y[0]}=1;
		}
	close infile0;
	open(provfasta,">$outcontig") || die;
	open(infilec,$contigsfna) || die;
	my $gd=1;
	while(<infilec>) {
		if($_=~/^>([^ ]+)/) {
			my $thiscontig=$1;
			if($singletonlist{$thiscontig}) { $gd=0; } else { $gd=1; }
			}
		if($gd) { print provfasta $_; }
		}
	close infilec;
	close provfasta;
	$contigsfna=$outcontig;	
	}

my $samdir="$datapath/sam";
opendir(indir,$samdir) || die;
my @samfiles=grep(/sam$/,readdir indir);
closedir indir;

my $command;
foreach my $sam(@samfiles) {
	my @w=split(/\./,$sam);
	if($skip{$w[1]}) { print "  Sample $w[1] flagged for nobinning, skipping\n"; next; }
	print "  Working for $sam\n";
	print outsyslog "  Working for $sam\n";
	my $samfile="$samdir/$sam";
	my $sorted="$tempdir/$sam";
	$sorted=~s/sam$/sorted.bam/;
	$command="$samtools_soft sort $samfile -o $sorted -\@$numthreads > /dev/null 2>&1";
	print outsyslog "  Creating sorted BAM for $sam: $command\n";
	print "  Creating sorted BAM for $sam\n";
	system $command;
	$command="$samtools_soft index $sorted > /dev/null 2>&1";
	print outsyslog "  Indexing sorted BAM: $command\n";
	print "  Indexing sorted BAM\n";
	system $command;
	}
	
print "\n  Cutting contigs in pieces!\n";
my $bedfile="$tempdir/$project.contigs.bed";
my $choppedcontigs="$tempdir/$project.choppedcontigs.fasta";
$command="python3 $concoct_dir/scripts/cut_up_fasta.py $contigsfna -c 10000 -o 0 --merge_last -b $bedfile > $choppedcontigs";
print outsyslog "\n  Cutting contigs in pieces!: $command\n";	
system $command;	
print "  Creating abundance table\n";	
$command="PATH=$installpath/bin:\$PATH python3 $concoct_dir/scripts/concoct_coverage_table.py $bedfile $tempdir/*.sorted.bam > $bindir/coverage_table.tsv";
print outsyslog "  Creating abundance table: $command\n";
system($command);
print "  Running concoct\n";
$command="OMP_THREAD_LIMIT=$numthreads python3 $concoct_dir/bin/concoct --composition_file $choppedcontigs --coverage_file $bindir/coverage_table.tsv --threads $numthreads  -b $bindir/concoct_int/ > /dev/null 2>&1";
print outsyslog "  Running concoct: $command\n";
system($command);
print "  Merging clusters\n";
$command="python3 $concoct_dir/scripts/merge_cutup_clustering.py $bindir/concoct_int/clustering_gt1000.csv > $bindir/concoct_int/clustering_merged.csv";
print outsyslog "  Merging clusters: $command\n";
system($command);
print "  Extracting final bins\n";
$command="python3 $concoct_dir/scripts/extract_fasta_bins.py $contigsfna $bindir/concoct_int/clustering_merged.csv --output_path $bindir > /dev/null 2>&1";
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
if(-e $outcontig) { system("rm $outcontig"); }
