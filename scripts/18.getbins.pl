#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Makes the Bins table

use strict;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$bincov,$contigcov,%bindirs,%dasdir,$contigsinbins,$resultpath,$bintable);
my(%bins,%contigs,%allsamples,%mapped,%totalreadcount,%taxrna);

	#-- Read 16S in contigs

my $rnafile="$resultpath/02.$project.16S.txt";
open(infile1,$rnafile) || warn "Cannot open $rnafile\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @t=split(/\t/,$_);
	my $contigid=$t[0]; 
	$contigid=~s/\_RNA\d+$//;
	$taxrna{$contigid}{$t[4]}++;
	}
close infile1;
	
	#-- Create the bin coverage table and the contigsinbins file

if(-e $bincov) { system("rm $bincov"); }
open(outfile1,">>$bincov") || die;
print outfile1 "#--Created by $0,",scalar localtime,"\n";
print outfile1 "# Bin ID\tMethod\tCoverage\tRPKM\tSample\n";

if(-e $contigsinbins) { system("rm $contigsinbins"); }
open(outfile2,">>$contigsinbins") || die;
print outfile2 "#--Created by $0,",scalar localtime,"\n";
print outfile2 "# Contig\tMethod\tBin ID\n";

foreach my $binmethod(sort keys %dasdir) {

	#-- For all the binning methods

	print "Method:$binmethod\n";
	my $bindir=$dasdir{$binmethod};
	my $checkmfile="$resultpath/17.$project.$binmethod.checkM";
	print "Reading checkM results in $checkmfile\n";
	
		#-- Read checkM results for each bin
	
	open(infile2,$checkmfile) || warn "Cannot find checkM results in $checkmfile\n";
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

	#-- Read data for each bin (tax, size, chimerism)

	opendir(indir1,$bindir) || die;
	my @files=grep(/tax$/,readdir indir1);
	closedir indir1;

	foreach my $tfil(@files) {
	my $bin=$tfil;
	$bin=~s/\.fa.tax|\.fasta.tax//g;
	print "Reading data for bin $bin           \r";
	open(infile3,"$bindir/$tfil") || die;
	while(<infile3>) {
 		chomp;
		if($_=~/^Consensus/) {
			my($consensus,$size,$chimerism)=split(/\t/,$_);
			$consensus=~s/Consensus\: //g;
			$size=~s/Total size\: //g;
			$chimerism=~s/Chimerism\: //g;
			$bins{$binmethod}{$bin}{consensus}=$consensus;
			$bins{$binmethod}{$bin}{size}=$size;
			$bins{$binmethod}{$bin}{chimerism}=$chimerism;
			}
		else { 
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
	open(infile4,"$bindir/$fasta");
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

	print "\nCalculating coverages\n";
	open(infile5,$contigcov) || die "Cannot open contig coverage file $contigcov\n";
	while(<infile5>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_); 
		my $bincorr=$contigs{$binmethod}{$k[0]}; 
		next if(!$bincorr);
		my $mappedbases=$k[1]*$k[3];
		my $sample=$k[$#k];
		$allsamples{$sample}=1;
		$mapped{$binmethod}{$bincorr}{$sample}{bases}+=$mappedbases;
		$mapped{$binmethod}{$bincorr}{$sample}{reads}+=$k[4];
		$totalreadcount{$binmethod}{$sample}+=$k[4];
		}
	close infile5;

	foreach my $binst(sort keys %{ $mapped{$binmethod} }) {
		foreach my $samps(sort keys %{ $mapped{$binmethod}{$binst} }) {
			my $coverage=$mapped{$binmethod}{$binst}{$samps}{bases}/$bins{$binmethod}{$binst}{size};
			my $rpkm=($mapped{$binmethod}{$binst}{$samps}{reads}*1000000000)/($bins{$binmethod}{$binst}{size}*$totalreadcount{$binmethod}{$samps});
			$bins{$binmethod}{$binst}{coverage}{$samps}=$coverage;
			$bins{$binmethod}{$binst}{rpkm}{$samps}=$rpkm;
			printf outfile1 "$binst\t$binmethod\t%.3f\t%.3f\t$samps\n",$coverage,$rpkm;
			}
		}

	}
		
	#-- Create table of results	
					   
	my $outputfile=$bintable;
	print "Creating table in $outputfile\n";
	open(outfile3,">$outputfile") || die;
	
	#-- Headers
	
	print outfile3 "# Created by $0, ",scalar localtime,"\n";
	print outfile3 "Bin ID\tMethod\tTax\tTax 16S\tSize\tGC perc\tNum contigs\tChimerism\tCompleteness\tContamination\tStrain Het";
	foreach my $countfile(sort keys %allsamples) { print outfile3 "\tCoverage $countfile\tRPKM $countfile"; }
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
			foreach my $countfile(sort keys %allsamples) { printf outfile3 "\t%.3f\t%.3f",$bins{$method}{$thisbin}{coverage}{$countfile},$bins{$method}{$thisbin}{rpkm}{$countfile}; }
			print outfile3 "\n";
		}
	}
close outfile3;
close outfile1;
close outfile2;

print "Done!\nTable created in $outputfile\n";

 
