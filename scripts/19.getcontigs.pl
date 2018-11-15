#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes the contig table, gathering data from previous results

use strict;
use Cwd;
use Tie::IxHash;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$alllog,$contigsfna,$aafile,$contigcov,$contigsinbins,$nobins,$contigtable,%bindirs,%dasdir);

my(%contig,%allsamples);

	#-- Reading taxonomic assignment and chimerism for the contigs

open(infile1,$alllog) || warn "Cannot open contiglog file $alllog\n";
print "Reading taxa for contigs information...";
while(<infile1>) { 
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	$contig{$t[0]}{tax}=$t[1]; 
	if($t[3]=~/Chimerism level\: (.*)/i) { $contig{$t[0]}{chimerism}=$1; }
}
close infile1;

	#-- Reading GC content and length of the contigs
	
print "done!\nReading GC & length... ";
open(infile2,$contigsfna) || warn "Cannot open fasta file $contigsfna\n";
my($thisname,$contigname,$seq);
while(<infile2>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {		#-- If we are reading a new contig, store the data for the last one
		$thisname=$1;
		if($contigname) {
			my $gc=gc_count($seq);
			$contig{$contigname}{gc}=$gc;
			$contig{$contigname}{len}=length $seq;
			}
		$seq="";
		$contigname=$thisname;
		}
 else { $seq.=$_; }			#-- Otherwise store the sequence of the current	
            }
close infile2;
if($contigname) {
	my $gc=gc_count($seq);
	$contig{$contigname}{gc}=$gc;
	$contig{$contigname}{len}=length $seq;
}

	#-- Reading number of genes for the contigs

print "done!\nReading number of genes... ";
open(infile3,$aafile) || warn "Cannot open aa file $aafile\n";
while(<infile3>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
		my $contigname=$1;
		$contigname=~s/\_\d+$//; 
		$contig{$contigname}{numgenes}++;
	}
}
close infile3;

  #-- Reading contig coverages 
  
print "done!\nReading coverages... ";
open(infile4,$contigcov) || die;
while(<infile4>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @cc=split(/\t/,$_);
	$contig{$cc[0]}{coverage}{$cc[$#cc]}=$cc[1];	#-- Coverage values
	$contig{$cc[0]}{rpkm}{$cc[$#cc]}=$cc[2];	#-- RPKM values
	$contig{$cc[0]}{raw}{$cc[$#cc]}=$cc[4];		#-- Raw read counts
	$allsamples{$cc[$#cc]}=1;
}
close infile4;  

  #-- Reading bins (if any)

if(!$nobins) {				#-- Skip this step if no bins were requested  
	print "done!\nReading bins... ";
	open(infile5,$contigsinbins);
	while(<infile5>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @e=split(/\t/,$_);
		$contig{$e[0]}{bin}{$e[1]}=$e[2];
		}
	close infile5;
}

	#-- CREATING CONTIG TABLE
	
print "done!\nCreating contig table...";
open(outfile1,">$contigtable") || die;

	#-- Headers

print outfile1 "#Created by $0, ",scalar localtime,"\n";
print outfile1 "Contig ID\tTax\tChimerism\tGC perc\tLength\tNum genes\tBin ID";
foreach my $countfile(sort keys %allsamples) { print outfile1 "\tCoverage $countfile\tRPKM $countfile\tRaw $countfile"; }
print outfile1 "\n";

	#-- Contig data

foreach my $p(sort keys %contig) { 
	my $binfield;
	next if(!$contig{$p}{numgenes});

	#-- bins

	if(!$nobins) {
		my $ld=0;
		$binfield="{\"Bins\": [";
		foreach my $binmet(sort keys %dasdir) { 
			if($contig{$p}{bin}{$binmet}) { 
				if($ld) { $binfield.=","; }
				$binfield.="{ \"$binmet\":\"$contig{$p}{bin}{$binmet}\" }"; 
				$ld=1;
				}
			}
		$binfield.="] }";					 
	 if(!$ld) { $binfield=""; }
	}	 					 

	#-- Output

	printf outfile1 "$p\t$contig{$p}{tax}\t$contig{$p}{chimerism}\t%.2f\t$contig{$p}{len}\t$contig{$p}{numgenes}\t$binfield",$contig{$p}{gc}; 
	foreach my $countfile(sort keys %allsamples) { printf outfile1 "\t%.3f\t%.3f\t%d",$contig{$p}{coverage}{$countfile},$contig{$p}{rpkm}{$countfile},$contig{$p}{raw}{$countfile}; }
	print outfile1 "\n";
}
close outfile1;

print "done!\nTable created in $contigtable\n";


#----------------------------- GC counting

sub gc_count {
	my $seq=shift;
	my @m=($seq=~/G|C/gi);
	my $gc=(($#m+1)/length $seq)*100;
	return $gc;
}
