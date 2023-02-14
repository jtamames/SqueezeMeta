#!/usr/bin/env perl

$|=1;

use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
my $sample=$ARGV[1];
my $par1name=$ARGV[2];
my $par2name=$ARGV[3];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;
my $command;

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$mode,$megahit_soft,$assembler_options,$extassembly,$contigid,$numthreads,$spades_soft,$flye_soft,$prinseq_soft,$mappingfile,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$scriptdir,$singletons,$methodsfile,$syslogfile,$norename);

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

system("rm -r $datapath/canu > /dev/null 2>&1");
$outassembly="$datapath/canu/contigs.fasta";
if($canumem eq "NF") {
	print "  Setting available memory for Canu\n";
	my %mem=get_mem_info;
	my $ram=$mem{"MemAvailable"}/(1024*1024);
	my $ramstr=sprintf('%.2f',$ram);
	$canumem=sprintf('%.0f',int($ram));
	$canumem*=0.8;
	print "  AVAILABLE (free) RAM memory: $ramstr Gb. We will set canu to use $canumem Gb.\n  You can override this setting using the -canumem option when calling SqueezeMeta.pl\n";
	print outsyslog "canumem set to $canumem (Free Mem $ramstr Gb)\n";
	}
$command="$canu_soft  -p $project -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0  maxThreads=$numthreads maxMemory=$canumem $assembler_options -nanopore-raw  $par1name > $syslogfile 2>&1;"; 
$command.="mv $datapath/canu/$project.contigs.fasta $outassembly"; 

my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

if($mode=~/coassembly|sequential/) { system("mv $outassembly $contigsfna"); }
else {
	my $provname="$interdir/01.$sample.fasta";
	system("mv $outassembly $provname");
	}
pprint outmet "Assembly was done using Canu (Koren et al 2017, Genome Res 27(5):722-36)\n";
