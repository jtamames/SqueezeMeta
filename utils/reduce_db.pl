#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 29/01/2019 Original version, (c) Javier Tamames, CNB-CSIC
#-- Last Common Ancestor (LCA) taxonomic assignment from a Diamond file. For blastx collapsed format 

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use strict;
use Tie::IxHash;
use Cwd;
use lib ".";
use threads;
use Getopt::Long;

my($numthreads,$listfile,$include,$exclude,$outfile,$download,$dbfile,$hel);

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my($symlinkpath,$symlinkdest,$utilsdir);
if(-l __FILE__)
        {
        $symlinkpath = dirname(__FILE__);
        $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $utilsdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$utilsdir/..");
my $scriptdir = "$installpath/scripts";
my $auxdir = "$installpath/lib/SQM_reads";
my $datadir = "$installpath/data/";
###

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

do "$scriptdir/SqueezeMeta_conf.pl";
our($databasepath);

my $start_run = time();
print BOLD "\nreduce_db v$version- (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;

my $helpshort="Usage: reduce_db.pl -i <nr database> -l <list of taxa> -o <output file> [options]\n";

my $helptext = <<END_MESSAGE;
Usage: reduce_db.pl -i <nr database> -l <list of taxa> -o <output file> [--include|--exclude]

Mandatory parameters:
   -i: Location of the nr database to use (REQUIRED)	
   -l: List of taxa to include or exclude (REQUIRED)
   -o: Output file for filtered database(REQUIRED)

 Options:
   --include: Use the list of taxa to INCLUDE entries corresponding to these taxa
   --exclude: Use the list of taxa to EXCLUDE entries corresponding to these taxa
   --download: Download a new nr database and store it in the file specified with the -i option
   -h: this help
   
END_MESSAGE

my $result = GetOptions (
			"l|list=s" => \$listfile,
			"include" => \$include,
			"exclude" => \$exclude,
			"download" => \$download,
			"o|output=s" => \$outfile,
			"i|input=s" => \$dbfile,
		     	"h" => \$hel
		  	  );



if($hel) { print "$helptext\n"; exit; }
if($include && $exclude) { die "$helpshort\nPlease select either --include or --exclude\n"; }
if(!$include) { $exclude=1; }
if(!$listfile) { die "$helpshort\nPlease specify a file containing taxa using -l option\n"; }
if(!$outfile) { die "$helpshort\nPlease provide a name for the output file\n"; }
if(!$dbfile) { die "$helpshort\nPlease provide a nr database\n"; }

if($download) {
	print "Downloading a new nr database\n";
	my $command="wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -nc -P $dbfile";
	system($command);
	}
elsif(-e $dbfile) {} 
else { die "$helpshort\nnr database does not exist in $dbfile, please provide a correct name\n"; }

my %taxlist;
my $ntax;
open(infile,$listfile) || die "Cannot open file $listfile contain taxa\n";
while(<infile>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	$taxlist{$_}=1;
	$ntax++;
	}
close infile; 
if(!$ntax) { die "File $listfile should contain taxa, but it is empty!\n"; }

my(%listd);
open(infile1,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my($tax,$par)=split(/\t/,$_);
	$par=~s/superkingdom\:/k_/; $par=~s/phylum\:/p_/; $par=~s/order\:/o_/; $par=~s/class\:/c_/; $par=~s/family\:/f_/; $par=~s/genus\:/g_/; $par=~s/species\:/s_/; $par=~s/no rank\:/n_/g; $par=~s/\w+\:/n_/g;
	foreach my $ilist(keys %taxlist) {
		if($par=~/$ilist/) { $listd{$tax}=$par; } #-- Taxa that falls in our list
		}
	}
close infile1;

print "Creating output file $outfile for ";
if($include) { print "INCLUDING "; } else { print "EXCLUDING "; }
print "entries in database $dbfile\n";
open(outf,">$outfile") || die "Cannot open output file $outfile\n";
my($good,$retained,$total);
if($dbfile=~/gz$/) { open(db,"zcat $dbfile |") || die "Cannot open input file $dbfile\n"; }
else { open(db,$dbfile) || die "Cannot open input file $dbfile\n"; }
while(<db>) {
	if($_=~/^\>/) {
		$good=0;
		$total++;
		my %accum=();
		my @spec=($_=~/\[.*?\]/g);
       		foreach my $k(@spec) {
                	$k=~s/\[|\]//g;
                	my @wd=split(/\s+/,$k);
                	next if(!$wd[1]);
			my $specname;
                	if($k=~/virus|phage|symbiont/i) { 
                        	$specname=$k;
                        	$specname=~s/\(.*//g;
                        	$specname=~s/ strain //g;
                        	$specname=~s/ genotype //g;
                        	$specname=~s/\s+$//;
                        	$specname=~s/\s+/ /g;
                         	}
                	elsif($k=~/Candidatus/) { $specname="$wd[0] $wd[1] $wd[2]"; } 
			elsif($k=~/bacterium$/) { $specname="$wd[0]"; }
                	else  { $specname="$wd[0] $wd[1]"; }  
                	$accum{$specname}++;
                }
		foreach my $u(keys %accum) { 
			if($listd{$u}) { 
				# print "$_ -> $u -> $listd{$u}\n"; 
				if($include) { $good=1; print outf $_; $retained++; last; }
				}
			elsif($exclude) { $good=1; print outf $_; $retained++; last; }				
			}
		}
		elsif($good) { print outf $_;  }
	if($total%10000==0) { print "Entries found: $total; Retained: $retained     \r"; }  	
	}
close db;
print "Entries found: $total; Retained: $retained          \n";		
print "Done! New database created in $outfile\n"; 
