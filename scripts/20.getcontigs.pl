#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes the contig table, gathering data from previous results

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$alllog,$contigsfna,$aafile,$contigcov,$contigsinbins,$nobins,$contigtable,%bindirs,%dasdir);

my(%contig,%allsamples);
tie %allsamples,"Tie::IxHash";

	#-- Reading taxonomic assignment and disparity for the contigs

open(infile1,$alllog) || warn "Can't open contiglog file $alllog\n";
print "Reading taxa for contigs information...";
while(<infile1>) { 
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @t=split(/\t/,$_);
	$contig{$t[0]}{tax}=$t[1]; 
	if($t[3]=~/Disparity\: (.*)/i) { $contig{$t[0]}{chimerism}=$1; }
}
close infile1;

	#-- Reading GC content and length of the contigs
	
print "done!\nReading GC & length... ";
open(infile2,$contigsfna) || warn "Can't open fasta file $contigsfna\n";
my($thisname,$contigname,$seq);
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
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
open(infile3,$aafile) || warn "Can't open aa file $aafile\n";
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	if($_=~/^\>([^ ]+)/) {
		my $contigname=$1;
		$contigname=~s/\_\d+\-\d+$//; 
		$contig{$contigname}{numgenes}++;
	}
}
close infile3;

  #-- Reading contig coverages 
  
print "done!\nReading coverages... ";
open(infile4,$contigcov) || die "Can't open $contigcov\n";
while(<infile4>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @cc=split(/\t/,$_);
	$contig{$cc[0]}{coverage}{$cc[$#cc]}=$cc[1];	#-- Coverage values
	$contig{$cc[0]}{rpkm}{$cc[$#cc]}=$cc[2];	#-- RPKM values
	$contig{$cc[0]}{tpm}{$cc[$#cc]}=$cc[3];	#-- TPM values
	$contig{$cc[0]}{raw}{$cc[$#cc]}=$cc[5];		#-- Raw read counts
	$allsamples{$cc[$#cc]}=1;
}
close infile4;  

  #-- Reading bins (if any)

if(!$nobins) {				#-- Skip this step if no bins were requested  
	print "done!\nReading bins... ";
	open(infile5,$contigsinbins); # File will be missing if running in sequential mode
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
open(outfile1,">$contigtable") || die "Can't open $contigtable for writing\n";

	#-- Headers

print outfile1 "#Created by $0, ",scalar localtime,"\n";
print outfile1 "Contig ID\tTax\tDisparity\tGC perc\tLength\tNum genes\tBin ID";
foreach my $countfile(keys %allsamples) { print outfile1 "\tCoverage $countfile\tRPKM $countfile\tTPM $countfile\tRAW READ COUNT $countfile"; }
print outfile1 "\n";

	#-- Contig data

my (@listcontigs,@sortedcontigs);
foreach my $ctg(keys %contig) {
	my @y=split(/\_/,$ctg);
	push(@listcontigs,{'contig',=>$ctg,'number'=>$y[1]});
	}
@sortedcontigs=sort {
	$a->{'number'} <=> $b->{'number'}
	} @listcontigs;

foreach my $ctg(@sortedcontigs) { 
	my $p=$ctg->{'contig'};
	my $binfield;
	#next if(!$contig{$p}{numgenes});

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
	foreach my $countfile(keys %allsamples) { printf outfile1 "\t%.3f\t%.3f\t%.3f\t%d",$contig{$p}{coverage}{$countfile},$contig{$p}{rpkm}{$countfile},$contig{$p}{tpm}{$countfile},$contig{$p}{raw}{$countfile}; }
	print outfile1 "\n";
}
close outfile1;

print "done!\n";
print "============\nCONTIG TABLE CREATED: $contigtable\n============\n\n";


#----------------------------- GC counting

sub gc_count {
	my $seq=shift;
	my @m=($seq=~/G|C/gi);
	my $gc=(($#m+1)/length $seq)*100;
	return $gc;
}
