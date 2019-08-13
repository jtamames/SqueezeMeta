#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Merges individual assemblies using minimus2, for merged mode. It also uses cd-hit for excluding identical contigs

use strict;
use warnings;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

my $project=$ARGV[0];
$project=~s/\/$//;
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; } 
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($resultpath,$interdir,$tempdir,$cdhit_soft,$extassembly,$minimus2_soft,$toamos_soft,$prinseq_soft,$numthreads);

#-- Merges the assemblies in a single dataset

my $finalcontigs="$resultpath/01.$project.fasta";
my($ecode,$command,$mergestep,$merged);
opendir(indir0,$interdir);
my @indassemblies=grep(/fasta$/,readdir indir0);
my $numassem=$#indassemblies+1;
closedir indir0;

if($extassembly) { 
	print "External assembly provided: $extassembly. Overriding assembly\n";
	system("cp $extassembly $finalcontigs");
	}
else {
	while($numassem>1) {
		$mergestep++;	
		if(-e $finalcontigs) { system("rm $finalcontigs"); }
		$merged="$interdir/merged_$mergestep.$project.fasta";
		my $tomerge=system("/media/mcm/jtamames/kmerdist.pl $project");
		open(infile0,"$tempdir/$project.2merge") || die;
		$_=<infile0>;
		chomp;
		my($sample1,$sample2,$twosamples,$mdist)=split(/\t/,$_);
		$sample1.=".fasta";
		$sample2.=".fasta";
		close infile0;
		print "MERGE $mergestep, $sample1 and $sample2\n";
		#opendir(indir0,$interdir);
		#my @indassemblies=grep(/fasta$/,readdir indir0);
		#my @n=reverse sort @indassemblies;
		#@indassemblies=@n;
		#closedir indir0;
		#my $numassemblies=$#indassemblies+1;
		#for(my $posmix=1; $posmix<=($numassemblies-1); $posmix++) {
		#	if($posmix==1) { 
		#		$command="cat $interdir/$indassemblies[0]  $interdir/$indassemblies[1] > $merged"; 
		#		print "Merging assemblies 1 ($indassemblies[0]) and 2 (",$indassemblies[1],")\n";
		#		print "**$command\n";
		#		}
		#	else { 
		#		$command="cat $merged $interdir/$indassemblies[$posmix] > $tempdir/interm.$project.fasta; mv $tempdir/interm.$project.fasta $merged"; 
		#		print "Merging assembly ",$posmix+1," ($indassemblies[$posmix])\n";
		#		print "**$command\n";
		#		}
		$command="cat $interdir/$sample1  $interdir/$sample2 > $merged";
		system $command;
		if(-z $merged) { die "$merged is empty\n"; }

			#-- Uses cd-hit to identify and remove contigs contained in others

			my $merged_clustered="$tempdir/mergedassemblies.$project.99.fasta";
			$command="$cdhit_soft -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9 > /dev/null 2>&1";
			print "Running cd-hit-est: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $merged_clustered) { die "$merged_clustered is empty\n"; }

			#-- Uses Amos to chage format to afg (for minimus2)

			my $afg_format="$tempdir/mergedassemblies.$project.99.afg";
			$command="$toamos_soft -s $merged_clustered -o $afg_format";
			print "Transforming to afg format: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Uses minimus2 to assemble overlapping contigs

			$command="$minimus2_soft $tempdir/mergedassemblies.$project.99 -D OVERLAP=100 -D MINID=95 -D THREADS=$numthreads > /dev/null 2>&1";
			print "Merging with minimus2: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Create the final result (overlapping contigs plus singletons)

			system("cat $tempdir/mergedassemblies.$project.99.fasta $tempdir/mergedassemblies.$project.99.singletons.seq > $merged.prov");

			open(outfile0,">$merged") || die;
			open(infile0,"$merged.prov") || die;
			while(<infile0>) {
				chomp;
				if($_=~/^\>(\d+)/) {
					my $newname="Merged\_$1\_$mergestep";
					print outfile0 ">$newname\n"; 
					}
				else { print outfile0 "$_\n"; }
				}
			close infile0;
			close outfile0;
			
			system("mv $interdir/$sample1 $interdir/$sample1.orig");
			system("mv $interdir/$sample2 $interdir/$sample2.orig");
			opendir(indir0,$interdir);
			@indassemblies=grep(/fasta$/,readdir indir0);
			$numassem=$#indassemblies+1;
			}
				

		#-- Remove files from temp

		 system("rm -r $tempdir/mergedassemblies*");
	}
system("mv $merged $finalcontigs");

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $finalcontigs -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
	

#-- Count length of contigs (needed later)

my $contigslen="$interdir/01.$project.lon";
print "Counting lengths\n";
open(outfile1,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile1,$finalcontigs) || die "Can't open $finalcontigs\n";
my($thisname,$contigname,$seq);
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
	$thisname=$1;
	if($contigname) {
		my $len=length $seq;
		print outfile1 "$contigname\t$len\n"; 
		}
	$seq="";
	$contigname=$thisname;
	}
	else { $seq.=$_; }
}
close infile1;
close outfile1;
