#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes eggnog database, ready for diamond usage
#-- Requires Diamond for formatting

use strict;

$|=1;

my $databasedir="/media/mcm/jtamames/databases"; 	#-- THIS MUST POINT TO THE DATABASES DIRECTORY
my $software_dir="/home/jtamames/software/";		#-- THIS MUST POINT TO THE DIAMOND DIRECTORY

my $eggnogdb="http://eggnogdb.embl.de/download/eggnog_4.5";	#-- SITE FOR DOWLOADING DATA
my $outegg="$databasedir/eggnog4.fasta";

my $command;
print "Downloading eggnog data from $eggnogdb\n";
$command="wget $eggnogdb/data/NOG/NOG.members.tsv.gz -nc -P $databasedir";
system $command;
print "Downloading sequence data from $eggnogdb\n";
$command="wget $eggnogdb/eggnog4.proteins.all.fa.gz -nc -P $databasedir";
system $command;

my %cogs;
open(infile1,"zcat $databasedir/NOG.members.tsv.gz |") || die; 
while(<infile1>) {
	chomp;
	my @f=split(/\t/,$_);
	my @e=split(/\,/,$f[5]);
	map { $cogs{$_}{$f[1]}=1; } @e;
            }
close infile1;

print "Now creating database\n";
open(outfile1,">$outegg") || die;
open(infile2,"zcat $databasedir/eggnog4.proteins.all.fa.gz |") || die; 
while(<infile2>) {
	my($badid,$cogstr);
	if($_=~/^\>(.*)/) {
		my $id=$1;
		chomp $id;
		$cogstr="";
		$badid=0;
		foreach my $p(sort keys %{ $cogs{$id} }) { $cogstr.="$p;"; }
		if(!$cogstr) { $badid=1; next; }
		chop $cogstr;
		#  print ">$id|$cogstr\n";
		print outfile1 ">$id|$cogstr\n";
		}
	elsif(!$badid) { print outfile1 $_; }
	}
close infile2;
close outfile1;
 
print "Removing datafiles\n";
system("rm $databasedir/eggnog4.proteins.all.fa.gz; rm $databasedir/NOG.members.tsv.gz;");
print "Formatting database\n";
system("$software_dir/diamond makedb --in $outegg -d $databasedir/eggnog -p 8");
print "FINISHED: eggnog database created in $outegg\n";
