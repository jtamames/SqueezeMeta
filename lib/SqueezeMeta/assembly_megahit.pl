#!/usr/bin/env perl

use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 

$|=1;

use File::Basename;
use Cwd 'abs_path';
our $sqmlibdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $sqmlibdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $sqmlibdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$sqmlibdir/../..");

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

system("rm -r $datapath/megahit > /dev/null 2>&1"); 
$outassembly="$datapath/megahit/final.contigs.fa";
if($par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }
else {  $command="$megahit_soft $assembler_options -r $par1name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }  #-- Support for single reads

my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

if($mode=~/coassembly|sequential/) { system("mv $outassembly $contigsfna"); }
else {
	my $provname="$interdir/01.$sample.fasta";
	system("mv $outassembly $provname");
	}
print outmet "Assembly was done using Megahit (Li et al 2015, Bioinformatics 31(10):1674-6)\n";
