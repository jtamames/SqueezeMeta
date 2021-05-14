#!/usr/bin/env perl

use Cwd;
use Getopt::Long;
use strict;
use File::Basename;
use Cwd 'abs_path';
use lib "."; 

$|=1;

my($project,$taxon,$funclass,$hel,$reqfunctions,$colorreq,$outfile);

#-- Define help text

my $helptext = <<END_MESSAGE;
Usage: sqm2ipath.pl -p project <options>

Options:

   -classification [cog|kegg]: Functional classification to use (Default: kegg)
   -taxon [string]: Taxon to be considered
   -functions [file]: File containing the name of the functions to be considered
   -color [string]: Color to plot in ipath
   -out [file]: Output file (Default: ipath.out)
   -h: This help
     
END_MESSAGE

my $result = GetOptions ("p=s" => \$project,
                     "taxon=s" => \$taxon,
                     "classification=s" => \$funclass,
		     "h" => \$hel,
		     "functions=s" => \$reqfunctions,
		     "color=s" => \$colorreq,
		     "o|out=s" => \$outfile,
 		    );

if($hel) { die "$helptext\n"; } 
if(!$funclass) { $funclass="KEGG"; }
if(!$project) { die "Please provide a valid project name or project path\n"; }
my $projectname=$project;
$projectname=~s/.*\///g;
if(!$outfile) { $outfile="ipath.out"; }

our($installpath,$extdatapath,$contigsinbins,$mergedfile,$aafile,$tempdir,$resultpath,$minpath_soft,$bintable,$extpath);
my $resultpath="$project/results";

my(%requested);
if($reqfunctions) {
	open(infile1,$reqfunctions) || die "Cannot open requested functions file $reqfunctions\n";
	print "Reading requested functions from $reqfunctions\n";
	while(<infile1>) {
		chomp;
		my($id,$attr)=split(/\t|\s+/,$_,2);
		if(!$attr) { 
			if($colorreq) { $attr=$colorreq; }
			else { $attr="#ff0000"; }
			}
		$requested{$id}=$attr;
		}
	close infile1;
	}

my(%accumfun)=();
my($cogpos,$keggpos,$taxpos,$funpos);
my $funfile="$resultpath/13.$projectname.orftable";
print "Reading functions in $funfile\n";
open(infile1,$funfile) || die;
$_=<infile1>;
$_=<infile1>;
my $headerf=$_; 
chomp $headerf;
my @headerf=split(/\t/,$_);
for(my $pos=0; $pos<=$#headerf; $pos++) {
	if($headerf[$pos] eq "COG ID") { $cogpos=$pos; }
	elsif($headerf[$pos] eq "KEGG ID") { $keggpos=$pos; }
	elsif($headerf[$pos] eq "Tax") { $taxpos=$pos; }
	}	
if($funclass=~/COG/i) {
	if(!$cogpos) { die "No COG results in table $funfile!\n"; }
	else{ $funpos=$cogpos; }
	}
if($funclass=~/KEGG/i) {
	if(!$keggpos) { die "No KEGG results in table $funfile!\n"; }
	else{ $funpos=$keggpos; }
	}
print "Looking for $taxon\n";
while(<infile1>) {
	chomp;
	my @line=split(/\t/,$_);
	next if !$_;
	# print "$line[$taxpos]\n";
	if($taxon && ($line[$taxpos]=~/$taxon/i)) { 
		my $tfun=$line[$funpos];
		$tfun=~s/\*$//;     # <-- Add here treatment for besthit/bestaver
		next if(!$tfun);
		if(!$reqfunctions || ($requested{$tfun})) { $accumfun{$tfun}++; }
		}
	}
close infile;

open(out,">$outfile") || die "Cannot open output file $outfile\n";
foreach my $thisfun(sort keys %accumfun) {
	my $thiscolor=$requested{$thisfun};	#-- COLOR: 1) Specified in functions file; 2) Specified in command line; 3) Default: Red
	if(!$thiscolor) {
		if($colorreq) { $thiscolor=$colorreq; }
		else { $thiscolor="#ff0000"; }
		}
	print out "$thisfun $thiscolor\n";
	}
close out;
print "Output results in $outfile. Now load these into iPath (https://pathways.embl.de)\n";

	
		



