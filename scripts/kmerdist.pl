#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. Compares metagenomes by their k-mer distance to provide a merging order 

use strict;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

my $project=$ARGV[0];
$project=~s/\/$//;
my $mergestep=$ARGV[1];
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; } 
do "$project/SqueezeMeta_conf.pl";

our($numthreads,$interdir,$tempdir,$resultpath,$kmerdb_soft);

  #-- Reading sequences

opendir(indir1,$interdir);
my @fastafiles=grep(/01.*?\.fasta$|merged.*?\.fasta$/,readdir indir1);
closedir indir1;

  #-- Writing samples file
  
my %toremove;
if($mergestep>1) {
	open(infile0,"$tempdir/mergelog") || die "Cannot open merge log file $tempdir/mergelog!\n";
	while(<infile0>) {
		chomp;
		next if !$_;
		$toremove{$_}=1;
		# print "*$_*\n";
		}
	close infile0;
	}

my $samples="$tempdir/samples.$project.txt";
open(out1,">$samples") || die;
foreach my $file(@fastafiles) { 
	# print "--$file--\n";
	next if($toremove{$file});
	# print "  ok\n";
	$file=~s/\.fasta.*//; 
	print out1 "$interdir/$file\n"; 
	}
close out1;

  #-- Running kmer-db
  
print "Calculating similarities between metagenomes using k-mer db\n";
my $command;
my $kmerdb="$tempdir/kmerdb.$project.txt";
$command="$kmerdb_soft build -t $numthreads $samples $kmerdb > /dev/null 2>&1";
#print "$command\n";
print "k-mer db: Building database\n";
my $ecode=system($command); 
if($ecode!=0) { die "Error running command:    $command"; }
my $kmertable="$tempdir/kmertable.$project.txt";
$command="$kmerdb_soft all2all -t $numthreads $kmerdb $kmertable > /dev/null 2>&1";
#print "$command\n";
print "k-mer db: Comparing metagenomes\n";
my $ecode=system($command);
if($ecode!=0) { die "Error running command:    $command"; }
my $disttable="$kmertable.jaccard";
 if(-e $disttable) { system("rm $disttable"); }
$command="$kmerdb_soft distance -t $numthreads $kmertable > /dev/null 2>&1";
#print "$command\n";
print "k-mer db: Calculating distances\n";
my $ecode=system($command);
if($ecode!=0) { die "Error running command:    $command"; }

  #-- Reading the distance file
  
open(infile1,$disttable) || die "Cannot open distance file in $disttable\n";
my $header;
my @fields;
my %distances;
while(<infile1>) {
	chomp;
	next if !$_;
	if(!$header) { 
		$header=$_;
		@fields=split(/\,/,$_);
		}
	else {
		my @t=split(/\,/,$_);
		my $s1=$t[0];
		for(my $pos=1; $pos<=$#t; $pos++) {
			my $s2=$fields[$pos];
			next if($s1 eq $s2);
			next if(!$t[$pos]);
			$distances{"$s1\t$s2"}=$t[$pos];
			}
		}
	}
close infile1;

  #-- Returns the maximum similarity value

my $tmerge="$tempdir/$project.2merge";
if(-e $tmerge) { system("rm $tmerge"); }

open(outfile2,">$tmerge") || die;
my @ordlist=sort { $distances{$b}<=>$distances{$a}; } keys %distances;
print outfile2 "$ordlist[0]\t$distances{$ordlist[0]}\n";
print "Maximum similarity for $ordlist[0]: $distances{$ordlist[0]}\n";
close outfile2;			
