#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes the Bins table

use strict;
use Cwd;
use lib ".";
use Tie::IxHash;

$|=1;

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$bincov,$contigcov,%bindirs,%dasdir,$contigsinbins,$resultpath,$interdir,$bintable,$syslogfile);
my(%bins,%contigs,%allsamples,%mapped,%totalreadcount,%taxrna,%rinsample);
tie %allsamples,"Tie::IxHash";
open(syslogfile,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

	#-- Read 16S in contigs

my $rnafile="$resultpath/02.$project.16S.txt";
print syslogfile "  Reading 16S rRNA in contigs from $rnafile\n";
open(infile1,$rnafile) || warn "Can't open $rnafile\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @t=split(/\t/,$_);
	my $contigid=$t[0]; 
	$contigid=~s/\_RNA\d+$//;
	my @cidtemp = split("_", $contigid);
	pop @cidtemp;
	$contigid = join("_", @cidtemp);
	$taxrna{$contigid}{$t[4]}++;
	}
close infile1;
	
	#-- Create the bin coverage table and the contigsinbins file

if(-e $bincov) { system("rm $bincov"); }
print syslogfile "  Creating bin coverage table in $bincov\n";
open(outfile1,">>$bincov") || die "Can't open $bincov for writing\n";
print outfile1 "#--Created by $0,",scalar localtime,"\n";
print outfile1 "# Bin ID\tMethod\tCoverage\tRPKM\tTPM\tSample\n";

if(-e $contigsinbins) { system("rm $contigsinbins"); }
print syslogfile "  Creating contigs in bins table in $contigsinbins\n";
open(outfile2,">>$contigsinbins") || die "Can't open $contigsinbins for writing\n";
print outfile2 "#--Created by $0,",scalar localtime,"\n";
print outfile2 "# Contig\tMethod\tBin ID\n";

foreach my $binmethod(sort keys %dasdir) {

	#-- For all the binning methods

	print "  Method:$binmethod\n";
	my $bindir=$dasdir{$binmethod};
	my $checkmfile="$interdir/18.$project.$binmethod.checkM";
	print "  Reading checkM results in $checkmfile\n";
	
		#-- Read checkM results for each bin
	
	open(infile2,$checkmfile) || warn "Can't open checkM results in $checkmfile\n";
	while(<infile2>) {
		chomp;
		next if !$_;
		if(($_=~/^\s+/) && ($_!~/Bin Id/)) {
		$_=~s/^\s+//g;
		$_=~s/Bacteria/Bacteria ()/g;
		$_=~s/Archaea/Archaea ()/g;
		my @k=split(/\s+/,$_);
		my $bin=$k[0]; 
		my $complete=$k[12];
		$bins{$binmethod}{$bin}{complete}=$k[12];
		$bins{$binmethod}{$bin}{contamination}=$k[13];
		$bins{$binmethod}{$bin}{strain}=$k[14];  
		}
	}
	close infile2;

	#-- Read data for each bin (tax, size, disparity)

	opendir(indir1,$bindir) || die "Can't open $bindir directory\n";
	my @files=grep(/tax$/,readdir indir1);
	closedir indir1;

	foreach my $tfil(@files) {
	my $bin=$tfil;
	$bin=~s/\.fa.tax|\.fasta.tax//g;
	print "  Reading data for bin $bin           \r";
	open(infile3,"$bindir/$tfil") || die "Can't open $bindir/$tfil\n";
	while(<infile3>) {
 		chomp;
		if($_=~/^Consensus/) {
			my($consensus,$size,$chimerism)=split(/\t/,$_);
			$consensus=~s/Consensus\: //g;
			$size=~s/Total size\: //g;
			$chimerism=~s/Disparity\: //g;
			$bins{$binmethod}{$bin}{consensus}=$consensus;
			$bins{$binmethod}{$bin}{size}=$size;
			$bins{$binmethod}{$bin}{chimerism}=$chimerism;
			}
		elsif($_=~/^\>/) {
			$_=~s/\>//;
			my @t=split(/\s+/,$_);
			$contigs{$binmethod}{$t[0]}=$bin; 
			print outfile2 "$t[0]\t$binmethod\t$bin\n";
			$bins{$binmethod}{$bin}{contignum}++; 
			if($taxrna{$t[0]}) {
				foreach my $cid(sort keys %{ $taxrna{$t[0]} }) { $bins{$binmethod}{$bin}{rna}{$cid}++; }
				}
			}
		}
	close infile3;

	#-- Calculate GC for the bin

	my $fasta=$tfil;
	$fasta=~s/\.tax//g;
	my $seq;
	open(infile4,"$bindir/$fasta") || die "Can't open $bindir/$fasta\n";
	while(<infile4>) { 
		chomp;
		if($_!~/\>/) { $seq.=$_; } 
		}
	close infile4;
	my @m=($seq=~/G|C/gi);
	my $gc=(($#m+1)/length $seq)*100;
	$bins{$binmethod}{$bin}{gc}=$gc;
	$seq="";	      
	}

	#-- Count coverages for the bins

	print "\n  Calculating coverages\n";
	print syslogfile "  Calculating coverages for bins from $contigcov\n";
	open(infile5,$contigcov) || die "Can't open contig coverage file $contigcov\n";
	while(<infile5>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_); 
		my $bincorr=$contigs{$binmethod}{$k[0]}; 
		next if(!$bincorr);
		my $mappedbases=$k[6];
		my $sample=$k[$#k];
		$allsamples{$sample}=1;
		$mapped{$binmethod}{$bincorr}{$sample}{bases}+=$mappedbases;
		$mapped{$binmethod}{$bincorr}{$sample}{reads}+=$k[5];
		$rinsample{$binmethod}{$sample}{$bincorr}+=$k[5];
		$totalreadcount{$binmethod}{$sample}+=$k[5];
		}
	close infile5;

	my(%rpk,%accumrpk);
	foreach my $samps(sort keys %{ $rinsample{$binmethod} }) {
		foreach my $binst(sort keys %{ $rinsample{$binmethod}{$samps} }) {
		my $longt=$bins{$binmethod}{$binst}{size};
		$rpk{$binmethod}{$binst}{$samps}=$mapped{$binmethod}{$binst}{$samps}{reads}/$longt;
		$accumrpk{$samps}+=$rpk{$binmethod}{$binst}{$samps};
		}
	$accumrpk{$samps}/=1000000;
	}
	
	foreach my $binst(sort keys %{ $mapped{$binmethod} }) {
		foreach my $samps(sort keys %{ $mapped{$binmethod}{$binst} }) {
			my $coverage=$mapped{$binmethod}{$binst}{$samps}{bases}/$bins{$binmethod}{$binst}{size};
			my $rpkm=($mapped{$binmethod}{$binst}{$samps}{reads}*1000000000)/($bins{$binmethod}{$binst}{size}*$totalreadcount{$binmethod}{$samps});
			my $tpm=$rpk{$binmethod}{$binst}{$samps}/$accumrpk{$samps};
			$bins{$binmethod}{$binst}{tpm}{$samps}=$tpm;
			$bins{$binmethod}{$binst}{coverage}{$samps}=$coverage;
			$bins{$binmethod}{$binst}{rpkm}{$samps}=$rpkm;
			printf outfile1 "$binst\t$binmethod\t%.2f\t%.2f\t%.2f\t$samps\n",$coverage,$rpkm,$tpm;
			# print "$binst**$samps**$mapped{$binmethod}{$binst}{$samps}{reads}**$bins{$binmethod}{$binst}{size}**$totalreadcount{$binmethod}{$samps}**\n";
			}
		}

	}
		
	#-- Create table of results	
					   
	my $outputfile=$bintable;
	print "  Creating table in $outputfile\n";
	print syslogfile "  Creating table in $outputfile\n";
	open(outfile3,">$outputfile") || die "Can't open $outputfile for writing\n";
	
	#-- Headers
	
	print outfile3 "# Created by $0, ",scalar localtime,"\n";
	print outfile3 "Bin ID\tMethod\tTax\tTax 16S\tLength\tGC perc\tNum contigs\tDisparity\tCompleteness\tContamination\tStrain het";
	foreach my $countfile(keys %allsamples) { print outfile3 "\tCoverage $countfile\tRPKM $countfile\tTPM $countfile"; }
	print outfile3 "\n";
	
	#-- Data
	
	foreach my $method(sort keys %bins) {
		foreach my $thisbin(sort { $bins{$method}{$b}{complete}<=>$bins{$method}{$a}{complete} } keys %{ $bins{$method} }) { 
			my $taxrna;
			foreach my $cd(sort keys %{ $bins{$method}{$thisbin}{rna} }) { 
				if($cd) { $taxrna.="$cd|"; }
				}
			chop $taxrna;
			printf outfile3 "$thisbin\t$method\t$bins{$method}{$thisbin}{consensus}\t$taxrna\t$bins{$method}{$thisbin}{size}\t%.2f\t$bins{$method}{$thisbin}{contignum}\t$bins{$method}{$thisbin}{chimerism}\t$bins{$method}{$thisbin}{complete}\t$bins{$method}{$thisbin}{contamination}\t$bins{$method}{$thisbin}{strain}",$bins{$method}{$thisbin}{gc};			   
			foreach my $countfile(keys %allsamples) { printf outfile3 "\t%.3f\t%.3f\t%.3f",$bins{$method}{$thisbin}{coverage}{$countfile},$bins{$method}{$thisbin}{rpkm}{$countfile},$bins{$method}{$thisbin}{tpm}{$countfile}; }
			print outfile3 "\n";
		}
	}
close outfile3;
close outfile1;
close outfile2;
close syslogfile;

print "  Done!\n";
print "=============\nBIN TABLE CREATED: $outputfile\n=============\n\n";

 
