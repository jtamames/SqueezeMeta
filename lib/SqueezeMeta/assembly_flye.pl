#!/usr/bin/perl

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

system("rm -r $datapath/flye > /dev/null 2>&1");
$outassembly="$datapath/flye/contigs.fasta";
$command="$flye_soft $assembler_options -o $datapath/flye --plasmids --meta --genome-size 2.6g --threads $numthreads --nano-raw $par1name > $syslogfile 2>&1; "; 
$command.="mv $datapath/flye/assembly.fasta $outassembly";

my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

if($mode=~/coassembly|sequential/) { system("mv $outassembly $contigsfna"); }
else {
	my $provname="$interdir/01.$sample.fasta";
	system("mv $outassembly $provname");
	}
	
print outmet "Assembly was done using Flye (Kolmogorov et al 2019, Nature Biotech 37, 540-546)\n";
