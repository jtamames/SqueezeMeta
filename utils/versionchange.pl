#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 30/01/2020 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs Concoct for binning

use strict;
use Cwd;
use lib ".";
use Term::ANSIColor qw(:constants);
use File::Basename;

$|=1;

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

our $scriptdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $scriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $scriptdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$scriptdir/..");
###

my $pwd=cwd();
my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;
my $basedir="$projectpath/..";

our($version,$resultpath,$basedir,$interdir,$bintable,$bintax,$bincov,$contigsinbins,$contigtable,$mode,$scriptdir,%bindirs);
if(!$version) { $version="Unknown older version"; }
open(infile1,"$installpath/version.txt") || warn "Cannot find version number in $installpath/version.txt\n";
my $thisversion=<infile1>;
chomp $thisversion;
close infile1;

if($version ne $thisversion) {
	my $logfile="$projectname/changes.$projectname.txt";
	my $binners=join(",",sort keys %bindirs);
	print RED; 
	print "WARNING: The results you want to re-run happen to be created with an older version of SqueezeMeta ($version), different from the current one ($thisversion)\n";
	print "To make these results compatible with the current version, I would need to make some restructuring and name changes.\n";
	print "I will keep track of all these changes for you.\n";
	print RESET;
	my $answer;
	while(!$answer) {
		print "Do you want to proceed? (y/n)? ";
		$answer=<STDIN>;
		chomp($answer);
		}
	$answer=~tr/A-Z/a-z/;	
	if($answer eq "n") { print RED; print "Uncompatible versions\n"; print RESET; die; }
	open(changes,">$logfile") || die;
	print changes "Migrating results from version \"$version\" to \"$thisversion\", ",scalar localtime,"\n";
	if($thisversion eq "1.2.0") { 
		my $bindir="$interdir/binners";
		if(-d $bindir) {} else { system("mkdir $bindir"); }
		if(-d "$resultpath/maxbin") { 
			my $maxbindir="$bindir/maxbin";			
			print "  Moving $resultpath/maxbin to $maxbindir\n";
			print changes "Moving $resultpath/maxbin to $maxbindir\n";
			system("mv $resultpath/maxbin $bindir");
			}
		if(-d "$resultpath/metabat2") { 
			my $metabatdir="$bindir/metabat2";			
			print "  Moving $resultpath/metabat2 to $metabatdir\n";
			print changes "Moving $resultpath/metabat2 to $metabatdir\n";
			system("mv $resultpath/metabat2 $bindir");
			}		
		if(-d "$resultpath/DAS") { 
			my $dasdir="$bindir/DAS";
			my $finalbindir="$resultpath/bins";			
			print "  Moving $resultpath/DAS to $dasdir\n";
			print changes "Moving $resultpath/DAS to $dasdir\n";
			system("mv $resultpath/DAS $bindir");
			print "  Moving final bins to $finalbindir\n";
			print changes "  Moving final bins to $finalbindir\n";
			system("mv $dasdir/$projectname\_DASTool\_bins/* $finalbindir");
			}
		my $binresultdir="$resultpath/bins";			
		if(-d $binresultdir) {} else { print "  Creating directory $binresultdir\n"; print changes "Creating directory $binresultdir\n"; system("mkdir $binresultdir"); }}
		my $newfile="$interdir/16.$projectname.bintax";
		if(-e $bintax) { print "  Moving file $bintax to $newfile\n"; print changes "Moving file $bintax to $newfile\n"; system("mv $bintax $newfile"); }
		my $newfile="$interdir/17.$projectname.checkM";
		my $checkmfile="$interdir/18.$projectname.DASTool.checkM";
		if(-e $checkmfile) { print "  Moving file $checkmfile to $newfile\n"; print changes "Moving file $checkmfile to $newfile\n"; system("mv $checkmfile $newfile"); }
		my $newfile="$resultpath/18.$projectname.bintable";
		if(-e $bintable) { print "  Moving file $bintable to $newfile\n"; print changes "Moving file $bintable to $newfile\n"; system("mv $bintable $newfile"); }
		my $newfile="$interdir/18.$projectname.bincov";
		if(-e $bincov) { print "  Moving file $bincov to $newfile\n"; print changes "Moving file $bincov to $newfile\n"; system("mv $bincov $newfile"); }
		my $newfile="$interdir/18.$projectname.contigsinbins";
		if(-e $contigsinbins) { print "  Moving file $contigsinbins to $newfile\n"; print changes "Moving file $contigsinbins to $newfile\n"; system("mv $contigsinbins $newfile"); }
		my $newfile="$resultpath/19.$projectname.contigtable";
		if(-e $contigtable) { print "  Moving file $contigtable to $newfile\n"; print changes "Moving file $contigtable to $newfile\n"; system("mv $contigtable $newfile"); }
		my $newfile="$resultpath/20.$projectname.kegg.pathways";
		my $keggpathfile="$resultpath/21.$projectname.kegg.pathways";
		if(-e $keggpathfile) { print "  Moving file $keggpathfile to $newfile\n"; print changes "Moving file $keggpathfile to $newfile\n"; system("mv $keggpathfile $newfile"); }
		my $newfile="$resultpath/20.$projectname.metacyc.pathways";
		my $cycpathfile="$resultpath/21.$projectname.metacyc.pathways";
		if(-e $cycpathfile) { print "  Moving file $cycpathfile to $newfile\n"; print changes "Moving file $cycpathfile to $newfile\n"; system("mv $cycpathfile $newfile"); }
		my $newfile="$resultpath/21.$projectname.stats";
		my $statsfile="$resultpath/22.$projectname.stats";
		if(-e $statsfile) { print "  Moving file $statsfile to $newfile\n"; print changes "Moving file $statsfile to $newfile\n"; system("mv $statsfile $newfile"); }
	
		print "\n  Changes to files and directories stored in file $logfile\n";


		print "  Now changing file SqueezeMeta_conf.pl\n";
		system("mv $projectpath/SqueezeMeta_conf.pl $projectpath/SqueezeMeta_conf.pl.old");
		open(outfile1,">$projectpath/SqueezeMeta_conf.pl") || die;
		open(infile1,"$projectpath/SqueezeMeta_conf.pl.old") || die "Cannot open model $projectpath/SqueezeMeta_conf.pl.old\n";
		print outfile1 "\$version = \"$thisversion\";\n";
		print outfile1 "\$mode = \"$mode\";\n\n";
		print outfile1 "\$installpath = \"$installpath\";\n";
		while(<infile1>) {
			chomp;
			if($_=~/^\$projectname/) { print outfile1 "\$projectname = \"$projectname\";\n"; }
			elsif($_=~/^\%bindirs/) { print outfile1 "\$binresultdir = \"\$resultpath/bins\";\n"; }
			elsif($_=~/^\%dasdir/) { }
			elsif($_=~/^\$bintax/) { $_=~s/17\./16\./; print outfile1 "$_\n"; print outfile1 "\$checkmfile=\"\$interdir/17.\$projectname.checkM\";\n"; }
			elsif($_=~/^\$bincov|^\$bintable|^\$contigsinbins/) { $_=~s/19\./18\./; print outfile1 "$_\n"; }
			elsif($_=~/^\$contigtable/) { $_=~s/20\./19\./; print outfile1 "$_\n"; }
			elsif($_=~/^\$mapper/) { print outfile1 "$_\n\$binners=\"$binners\";\n"; }
			elsif($_=~/^\$binners/) { }
			else { print outfile1 "$_\n"; }
			}
		close infile1;
		close outfile1;
		print "  New SqueezeMeta_conf.pl written to $basedir/$projectname/SqueezeMeta_conf.pl\n"; 
		print changes "New SqueezeMeta_conf.pl written to $basedir/$projectname/SqueezeMeta_conf.pl\n";	
		close changes;
	
		}
