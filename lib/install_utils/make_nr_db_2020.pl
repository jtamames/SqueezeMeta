#!/usr/bin/env perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes nr database (and nr reduced), ready for diamond usage
#-- Requires package ncbi-blast+ (for blastdbcmd)
#
# - 16-05-2018: Modify to use relative paths (Fernando Puente-Sánchez).
# - 06-06-2018: Add scriptdir patch FPS.
# - 10-02-2020: Directly download the fasta version of nr from NCBI. Drop the functionality for making a reduced database. FPS

use strict;
###scriptdir patch, Fernando Puente-Sánchez, 07-V-2018
use File::Basename;
our $dbscriptdir = dirname(__FILE__);
our $installpath = "$dbscriptdir/../..";
###

$|=1;


my $databasedir=$ARGV[0];			        #-- THIS MUST POINT TO THE DATABASES DIRECTORY

my $fastadb="$databasedir/nr.faa";			#-- Name of the fasta file to create
my $dbfile="$databasedir/nr.dmnd";
my $md5file="$databasedir/nr.md5";
my $bindir="$installpath/bin";

#-- Getting the raw files from NCBI. This can take long and need 100 Gb disk space
my $command="wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz -P $databasedir";
system $command;

#-- Format the database
system("gunzip $databasedir/nr.gz && mv $databasedir/nr $fastadb");

system("$bindir/diamond makedb --in $fastadb -d $databasedir/nr -p 8");
system("md5sum $dbfile > $md5file");
