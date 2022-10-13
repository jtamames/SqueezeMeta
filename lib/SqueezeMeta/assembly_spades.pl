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

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$mode,$megahit_soft,$assembler_options,$extassembly,$contigid,$numthreads,$spades_soft,$flye_soft,$prinseq_soft,$mappingfile,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$scriptdir,$singletons,$methodsfile,$syslogfile,$norename);

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

system("rm -r $datapath/spades > /dev/null 2>&1");
if($assembler eq "spades") { $outassembly="$datapath/spades/contigs.fasta"; }
elsif($assembler eq "rnaspades") { $outassembly="$datapath/spades/transcripts.fasta"; }
my($command,$command_start);
if($assembler eq 'spades') { $command_start="$spades_soft --meta"; }
elsif($assembler eq 'rnaspades') { $command_start="$spades_soft --rna" }
if($par2name) { $command="$command_start --pe1-1 $par1name --pe1-2 $par2name -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile 2>&1"; }
else { $command="$command_start --s1 $par1name  -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile"; } #-- Support for single reads

my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

if($mode=~/coassembly|sequential/) { system("mv $outassembly $contigsfna"); }
else {
	my $provname="$interdir/01.$sample.fasta";
	system("mv $outassembly $provname");
	}

print outmet "Assembly was done using SPAdes (Bankevich et al 2012, J Comp Biol 19(5):455-77)\n";
