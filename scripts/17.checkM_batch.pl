#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs checkM for evaluating bins
#
#-- FPS (30-V-2018). Override $PATH for external dependencies of checkm.

use strict;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($installpath,$datapath,$taxlist,%bindirs,%dasdir,$checkm_soft,$alllog,$resultpath,$tempdir,$numthreads);

my $minsize=20000;  #-- Minimum size of a bin to be considered
my $markerdir="$datapath/checkm_markers";
my $checktemp="$tempdir/checkm_batch";

if(-d $markerdir) {} else { system "mkdir $markerdir"; }
if(-d $checktemp) {} else { system "mkdir $checktemp"; print "Creating $checktemp\n";  }

my(%tax,%bins,%consensus,%alltaxa);

#-- Read NCBI's taxonomy 

open(infile1,$taxlist) || die "Cannot find taxonomy list $taxlist\n";
print "Reading $taxlist\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	$tax{$t[1]}=$t[2];
	}
close infile1;

	#-- Read bin directories

foreach my $binmethod(sort keys %dasdir) {
	my $bindir=$dasdir{$binmethod};
	print "Looking for $binmethod bins in $bindir\n";
	
	#-- Read contigs in bins
	
	open(infile2,$alllog) || die;
	while(<infile2>) {
		chomp;
		next if !$_;
		my @t=split(/\t/,$_);
		$tax{$t[0]}=$t[1];
		}
	close infile2;

	#-- Read all bins

opendir(indir,$bindir);
my @files=grep(/tax$/,readdir indir);
closedir indir;


my $checkmfile="$resultpath/17.$project.$binmethod.checkM";	#-- From checkM_batch.pl, checkM results
if(-e $checkmfile) { system("rm $checkmfile"); }

	#-- Working for each bin

foreach my $m(@files) { 
	#if($exclude{$m}) { print "**Excluding $m\n"; next; }
	my $binname=$m;
	my $binname=~s/\.tax//g;
	my $thisfile="$bindir/$m";
	$bins{$thisfile}=$binname;
	print "Working for bin $thisfile\n";
 
	#-- Reading the consensus taxa for the bin
 
	open(infile3,$thisfile) || die;
	while(<infile3>) { 
		chomp;
		if($_=~/Consensus/) {
			my($cons,$size,$chim,$chimlev)=split(/\t/,$_);
			$cons=~s/Consensus\: //;
			$size=~s/Total size\: //g;
			if($size<$minsize) { print "Skipping bin $bins{$thisfile} because of low size\n"; next; }
			$consensus{$thisfile}=$cons;
			my @k=split(/\;/,$cons);
		
			#-- We store the full taxonomy for the bin because not all taxa have checkm markers
		
			foreach my $ftax(reverse @k) { 
				my($ntax,$rank);
				if($ftax!~/\:/) { $ntax=$ftax; } else { ($rank,$ntax)=split(/\:/,$ftax); }
				$ntax=~s/unclassified //gi;
				$ntax=~s/ \<.*\>//gi; 
				if($tax{$ntax} && ($tax{$ntax} ne "species")  && ($tax{$ntax} ne "no rank")) { 
				push( @{ $alltaxa{$thisfile} },"$tax{$ntax}:$ntax");
				#   print "$m\t$ntax\t$tax{$ntax}\n";
				}
			}
		}
	}
	close infile3;

	my $inloop=1;

	#-- We will find the deepest taxa with checkm markers
	while($inloop) {  
		my $taxf=shift(@{ $alltaxa{$thisfile}  }); 
		if(!$taxf) { last; $inloop=0; }
		my($rank,$tax)=split(/\:/,$taxf);
		$tax=~s/ \<.*//g;
		if($rank eq "superkingdom") { $rank="domain"; }
		print "Using profile for rank $rank: $tax for $thisfile\n";   
		my $marker="$markerdir/$tax.ms"; 
	
		#-- Use already existing tax profile or create it
	
		if(-e $marker) {} else { 
			my $command="export PATH=\"$installpath/bin/pplacer\":\$PATH; $checkm_soft taxon_set $rank $tax $marker > /dev/null"; #Override $PATH for external dependencies of checkm. (FPS).
                        my $ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			}
	
		#-- If it was not possible to create the profile, go for the upper rank
		
		if(-e $marker) {} else { next; }
	
		my $fastafile=$thisfile;
		$fastafile=~s/\.tax//;
		$fastafile=~s/.*\///;
		# print ">>> $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null\n";
		# system("$checkm_soft analyze -t $numthreads -x $bins{$thisfile} $marker $bindir $checktemp > /dev/null");
		my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null";
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }

		my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft qa -t $numthreads $marker $checktemp >> $checkmfile"; #Override $PATH for external dependencies of checkm. (FPS).
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
	#	system("rm -r $checktemp");
		$inloop=0;
		}
 	} 
}
