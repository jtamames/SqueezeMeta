#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/06/2022 Original version, (c) Javier Tamames, CNB-CSIC
#-- Filters a Pfam file selecting or removing the Pfams specified in a external file

$|=1;

use strict;
use Getopt::Long;

my $helptext = <<END_MESSAGE;
Usage: filter_pfam.pl [arguments]

Arguments:

  -pfam <pfam dir>: Pfam file to filter(MANDATORY)
  -filter <file name>: File with Pfam ids/names/accs to select (MANDATORY)
  -o|output <file name>: Resulting filtered Pfam file (OPTIONAL. Default: Pfam_filtered.hmm)
  -h: This help	

END_MESSAGE

my($pfamfile,$filterfile,$outputfile,$hel);
my $result = GetOptions ("pfam=s" => \$pfamfile,
		     "filter=s" => \$filterfile,
		     "o|output=s" => \$outputfile,
		     "h" => \$hel
		     );
		     
if($hel) { die "$helptext\n"; } 
if(!$filterfile) { die "Please provide a file with Pfam ids/names/accs to select\n"; }
if(!$pfamfile) { die "Please provide a Pfam file to filter\n"; }
if(!$outputfile) { $outputfile="Pfam_filtered.hmm"; }
	     
my %selected;
my($picked,$excluded);
open(inf,$filterfile) || die "Cannot open file $filterfile with Pfam ids/names/accs to select\n"; 
while(<inf>) {
	chomp;
	next if !$_;
	if($_=~/^\!(.*)/) { $selected{$1}=1; $excluded++; }
	else { $selected{$_}=1; $picked++; }
	}
close inf;

if($excluded && $picked) { die "Please select Pfam either to select OR to skip, not both\n"; }
if($excluded) { print "Excluding $excluded Pfams\n"; }
if($picked) { print "Selecting $picked Pfams\n"; }

open(outf,">$outputfile") || die "Cannot create output file $outputfile\n";
open(pf,$pfamfile) || die "Cannot open file $pfamfile to filter\n";
my($thisentry,$into,$count);
$into=0;
while(<pf>) {	
	$thisentry.=$_;
	chomp;
	next if !$_;
	if(!$into) { 
		if($_=~/^NAME\s+(.*)/) { $into=check($1); }
		if($_=~/^ACC\s+(.*)/) { $into=check($1); }
		if($_=~/^DESC\s+(.*)/) { $into=check($1); }
		}
	if($_=~/^\/\//) {
		if($into && $picked) { print outf $thisentry; $count++;  }
		if((!$into) && $excluded) { print outf $thisentry; $count++;  }
		$thisentry="";
		$into=0; 		
		}		
	}
close pf;

print "Resulting Pfams: $count\n";
print "Output file: $outputfile\n";

sub check {
	my $lf=shift;
	my $i=0;
	foreach my $ter(keys %selected) {
		if($lf=~/$ter/) { $i=1; }
		}
	return $i;	    
	}


