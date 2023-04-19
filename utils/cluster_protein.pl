#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs barrnap for RNA prediction. Sequentially runs model for bacteria, archaea, eukaryote and mitochondria
#-- Runs aragorn for tRNA/tmRNA prediction (21/11/2019)

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($installpath,$resultpath,$interdir,$tempdir,$mode,$aafile,$methodsfile,$syslogfile);
my $mmseqs_soft="/home/tamames/anaconda3/envs/SqueezeMeta1.6/bin/mmseqs";

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my @seqfiles;
opendir(indir,$interdir) || die;
@seqfiles=grep(/03.*.faa$/,readdir indir);
closedir indir;
# print "**@seqfiles**\n";

my $allfasta="$tempdir/allfaa.fasta";
open(out,">$allfasta") || die "Cannot open file $allfasta for writing\n";
foreach my $rfiles(@seqfiles) {
	my @k=split(/\./,$rfiles);
	my $project=$k[1];
	my $thisfile="$interdir/$rfiles";
	print "  Reading aa sequences from $thisfile\n";
	open(inf,$thisfile) || die "Cannot open $thisfile\n";
	while(<inf>) {
		chomp;
		$_=~s/^\>.*?\_/>$project\_/;
		$_=~s/\s+.*//;
		$_=~s/\*//;
		print out "$_\n";
		}
	close inf;
	}
print "  Sequences written to $allfasta\n";
close out;

print "\n  Now using mmseqs2 (Steinegger M., Soeding J. Nat. Biotech. 35, 1026-8, 2017) for clustering sequences\n";
my $dbfile="$tempdir/allfaa.db";
my $clusterfile="$tempdir/allfaa.clus";
my $repfile="$tempdir/allfaa.clus.rep";
my $repfasta="$tempdir/allfaa.clus.rep.fasta";
my $mapfile="$tempdir/allfaa.clus.rep.tsv";
print "  Creating database\n";
my $command="$mmseqs_soft createdb $allfasta $dbfile > /dev/null 2>&1";
system($command);
print "  Clustering sequences\n";
my $command="$mmseqs_soft linclust --cov-mode 0 -c 0.8 --min-seq-id 0.3 $dbfile $clusterfile $tempdir > /dev/null 2>&1";
system($command);
print "  Choosing representatives\n";
my $command="$mmseqs_soft createsubdb $clusterfile $dbfile $repfile  > /dev/null 2>&1";
system($command);
print "  Creating fasta for representatives\n";
my $command="$mmseqs_soft convert2fasta $repfile $repfasta > /dev/null 2>&1";
system($command);
print "  Creating mappings\n";
my $command="$mmseqs_soft createtsv $dbfile $dbfile $clusterfile $mapfile > /dev/null 2>&1";
system($command);

$command="mv $repfasta $aafile";
system($command);
