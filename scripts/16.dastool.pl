#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/11/2018 Version 0.3.1, (c) Javier Tamames, CNB-CSIC
#-- Runs DasTool for combining binning results

use strict;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

#-- Configuration variables from conf file

our($installpath,$datapath,$databasepath,$resultpath,$aafile,$contigsfna,%bindirs,$contigcov,$dastool_soft,$alllog,$tempdir,$score_tres16,$numthreads);

my $daspath="$resultpath/DAS";
system("mkdir $daspath");

#-- Creating contigs in bins tables

print "Creating tables of contigs in bins... ";
my($tables,$methods,$thiseq);
my @files;
foreach my $binmethod(sort keys %bindirs) {
	my $bindir=$bindirs{$binmethod};
	$tables.="$daspath/$binmethod.table,";
	$methods.="$binmethod,";
	opendir(indir1,$bindir) || die "Can't open $bindir directory\n";
	my @fastafiles=grep(/fasta$|fa$/,readdir indir1);
	closedir indir1;
	open(outfile1,">$daspath/$binmethod.table") || die "Can't open $daspath/$binmethod.table for writing\n";
	foreach my $tfil(@fastafiles) {
		my $bin=$tfil;
		$bin=~s/\.fa.tax|\.fasta.tax//g;
		open(infile1,"$bindir/$tfil") || die "Can't open $bindir/$tfil\n";
		while(<infile1>) { 
 			chomp;
			if($_=~/^>([^ ]+)/) { 
				$thiseq=$1; 
				if($thiseq) { print outfile1 "$thiseq\t$bin\n"; }
				}
			}
		close infile1;
		}
	close outfile1;
	}
chop $tables;
chop $methods;
print "done\n";

#-- Run DAS tool

#my $das_command="PATH=$installpath/bin:\$PATH $dastool_soft -i $tables -l $methods -c $contigsfna --write_bins 1 --proteins $aafile --score_threshold $score_tres16 --search_engine diamond -t $numthreads -o $resultpath/DAS/$project --db_directory $databasepath"; 
my $das_command="LD_LIBRARY_PATH=$installpath/lib PATH=$installpath/bin:\$PATH $dastool_soft -i $tables -l $methods -c $contigsfna --write_bins 1 --score_threshold $score_tres16 --search_engine diamond -t $numthreads -o $resultpath/DAS/$project --db_directory $databasepath";
 
print "Running DAS Tool for $methods\n";
 print "$das_command\n";
my $ecode = system $das_command;
if($ecode!=0) { die "Error running command:    $das_command"; }

