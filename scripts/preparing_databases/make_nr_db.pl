#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes nr database (and nr reduced), ready for diamond usage
#-- Requires package ncbi-blast+ (for blastdbcmd)
#
# - 16-05-2018: Modify to use relative paths (Fernando Puente-Sánchez).

use strict;

$|=1;

my $reduce=0;		#-- Set to 1, make reduced database (not eukaryotes). Set to 0, full database

my $databasedir=$ARGV[0];			#-- THIS MUST POINT TO THE DATABASES DIRECTORY

# my $databasedir="/media/mcm/jtamames/temp/db";		#-- THIS MUST POINT TO THE DATABASES DIRECTORY
my $fastadb="$databasedir/nr.faa";			#-- Name of the fasta file to create
my $reducedfastadb="$databasedir/nr_reduced.fasta";	#-- Reduced database (No eukaryotes)
my $taxdict="$databasedir/../data/taxdict.txt";	#-- Resulting from parsing NCBI's taxonomy
my $bindir="$databasedir/../bin";
	#-- Getting the raw files from NCBI. This can take long and need 100 Gb disk space

my $command="wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr*gz -nc -P $databasedir";
system $command;

	#-- Read the downloaded files and changing its format

opendir(indir1,$databasedir) || die;
my @files=grep(/tar.gz/,readdir indir1);
closedir indir1;

foreach my $tfile(sort @files) {
	my @d=split(/\./,$tfile);
	my $rootname="$d[0]\.$d[1]";
	my $command="tar xvzf $databasedir/$tfile -C $databasedir; $bindir/blastdbcmd -entry all -outfmt %f -db $databasedir/$rootname -out $databasedir/$rootname.fasta"; 
	print "$command\n"; 
	system $command; 
	system("rm $databasedir/$rootname.p*");		#-- Deleting unnecessary files
	}
	
	#-- Join all files
	
system("cat $databasedir/*fasta > $fastadb");
if((-e $fastadb)  && (!(-z $fastadb))) { system("rm $databasedir/nr*tar.gz;rm $databasedir/nr*fasta"); } else { die "ERROR: No results in $fastadb\n"; }

	#-- Format the database

system("$bindir/diamond makedb --in $fastadb -d $databasedir/nr -p 8");

	#-- Make the reduced database

if($reduce) {
	print "\nMaking a REDUCED database\n";
	my(%euk,%prok);
	
	#-- Reading taxonomy for discriminating Eukaryotes (reduced database)

	open(infile1,$taxdict) || die;
	while(<infile1>) {
		chomp;
		my @t=split(/\t/,$_); 
		my @umm=split(/\s+/,$t[1]);
		my $wd=$umm[0];
		$wd=~s/\W+//g; 
		if($t[3] eq "E") { $euk{$wd}=1;  } else { $prok{$wd}=1;  }
		}
	close infile1;
	open(infile2,$fastadb) || die;
	
	#-- Writing sequences
	
	open(outfile1,">$reducedfastadb");
	while(<infile2>) {
		chomp;
		next if !$_;
		my $good;
		if($_=~/^\>/) {
			$good=0;
			my @mlist=($_=~/\[(.*?)\]/g);
			foreach my $taxname(@mlist) {
				my @words=split(/\W+/,$taxname);
				foreach my $tw(@words) { 
					next if($tw=~/unclassified\candidatus/i);
					next if($tw!~/^[A-Z]/); 
					if($prok{$tw}) { $good=1; last; }
				}
				last if $good;			 
			}
		}
		if($good) { print outfile1 "$_\n"; }
	}
	close infile2;
	close outfile1;
	
	system("$bindir/diamond makedb --in $reducedfastadb -d $databasedir/nr_reduced -p 8");
    }   
