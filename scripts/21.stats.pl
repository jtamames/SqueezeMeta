#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades) for several metagenomes that will be merged in the next step (merged mode).
#-- Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$tempdir,$basedir,$prinseq_soft,$mincontiglen,$resultpath,$contigsfna,$contigtable,$mergedfile,$bintable,$evalue,$miniden,$mincontiglen,$assembler,$mode);

my %sampledata;

my $resultfile="$resultpath/21.$project.stats";
open(outfile1,">$resultfile") || die "Cannot open $resultfile\n";

	#-- Read bases and reads

my @ranks=('superkingdom','phylum','class','order','family','genus','species');
my($totalbases,$totalreads);
my $mapfile="$resultpath/09.$project.mappingstat";
open(infile1,"$resultpath/09.$project.mappingstat") || die "Cannot open $mapfile\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	$sampledata{$k[0]}{reads}=$k[1];
	$sampledata{$k[0]}{bases}=$k[4];
	$totalbases+=$k[4];
	$totalreads+=$k[1];
	}
close infile1;

	#-- Statistics on contigs with prinseq
	
my(%contigs,%contax);
my $command="$prinseq_soft -stats_info -stats_assembly -stats_len -fasta $contigsfna > $tempdir/stats.txt";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
open(infile2,"$tempdir/stats.txt");
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	if($k[1] eq "N50") { $contigs{N50}= $k[2]; }
	if($k[1] eq "N90") { $contigs{N90}= $k[2]; }
	elsif($k[1] eq "bases") { $contigs{bases}= $k[2]; }
	elsif($k[1] eq "reads") { $contigs{num}= $k[2]; }
	elsif($k[1] eq "max") { $contigs{max}= $k[2]; }
	elsif($k[1] eq "min") { $contigs{min}= $k[2]; }
	}
close infile2;	

	#-- Statistics on contigs (disparity, assignment..)

open(infile3,$contigtable) || die "Cannot open $contigtable\n";
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my @mtax=split(/\;/,$k[1]);
	foreach my $p(@mtax) {
		my($trank,$ttax)=split(/\:/,$p);
		$contigs{$trank}++;
		$contax{$trank}{$ttax}++;
		}
	if($k[2]==0) { $contigs{chimerism}{0}++; }
	else { $contigs{chimerism}{more0}++; }
	if($k[2]>=0.1) { $contigs{chimerism}{0.1}++; }
	if($k[2]>=0.25) { $contigs{chimerism}{0.25}++; }
	}
close infile3;

	#-- Statistics on genes (disparity, assignment..)
		
my $header;
my @head;
my %genes;
open(infile4,$mergedfile) || die "Cannot open $mergedfile\n";
while(<infile4>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	if(!$header) { $header=$_; @head=split(/\t/,$header); next; }
	$genes{totgenes}++;
	my @k=split(/\t/,$_);
	if($k[6]) { $genes{kegg}++; }
	if($k[9]) { $genes{cog}++; }
	if($k[12]) { $genes{pfam}++; }
	if($k[0]=~/\_RNA/) { $genes{rnas}++; }
	foreach(my $pos=0; $pos<=$#k; $pos++) {
		my $f=$head[$pos];
		if($f=~/RAW COUNTS (.*)/) {
			my $tsam=$1;
			if($k[$pos]>0) {
				$genes{$tsam}{totgenes}++;
				if($k[6]) { $genes{$tsam}{kegg}++; }
				if($k[9]) { $genes{$tsam}{cog}++; }
				if($k[12]) { $genes{$tsam}{pfam}++; }
				if($k[0]=~/\_RNA/) { $genes{$tsam}{rnas}++; }
				}
			}
		}
	}
close infile4;

	#-- Statistics on bins (disparity, assignment..)

my %bins;
if($mode ne "sequential") {
	my $header;
	if(-e $bintable) {
		open(infile5,$bintable) || die "Cannot open $bintable\n";
		while(<infile5>) {
			chomp;
			next if(!$_ || ($_=~/^\#/));
			if(!$header) { $header=$_; next; }
			my @k=split(/\t/,$_);
			my $method=$k[1];
			$bins{$method}{num}++;
			my @mtax=split(/\;/,$k[1]);
			foreach my $p(@mtax) {
				my($trank,$ttax)=split(/\:/,$p);
				$bins{$method}{$trank}++;
				}
			if($k[8]>=50) { $bins{$method}{complete}{50}++; }
			if($k[8]>=75) { $bins{$method}{complete}{75}++; }
			if($k[8]>=90) { $bins{$method}{complete}{90}++; }
			if($k[9]<10) { $bins{$method}{contamination}{10}++; }
			if($k[9]>=50) { $bins{$method}{contamination}{50}++; }
			if($k[7]==0) { $bins{$method}{chimerism}{0}++; }
			else { $bins{$method}{chimerism}{more0}++; }
			if($k[7]>=0.1) { $bins{$method}{chimerism}{0.1}++; }
			if($k[7]>=0.25) { $bins{$method}{chimerism}{0.25}++; }
			if(($k[8]>=90) && ($k[9]<10)) { $bins{$method}{hiqual}++; }
			if(($k[8]>=75) && ($k[9]<10)) { $bins{$method}{goodqual}++; }
			}
		close infile5;
		}
	}

	#-- Date of the start of the run

open(infile6,"$basedir/$project/syslog");
my $startdate=<infile6>;
chomp $startdate;
close infile6; 
		

	#-- Write output

	#-- Data for the run
	
print outfile1 "# Created by $0 (SqueezeMeta), ",scalar localtime,"\n";
print outfile1 "# $startdate\n";
print outfile1 "\n";
print outfile1 "# Parameters\n";
print outfile1 "Assembler: $assembler\n";
print outfile1 "Min length of contigs: $mincontiglen\n";
print outfile1 "Max evalue for Diamond search: $evalue\n";
print outfile1 "Min identity for Diamond COG/KEGG search: $miniden\n";
print outfile1 "\n";

	#-- Reads
	
print outfile1 "# Statistics on reads\n";
print outfile1 "#\tAssembly";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$sample"; }
print outfile1 "\n";
print outfile1 "Number of reads\t$totalreads";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$sampledata{$sample}{reads}"; }
print outfile1 "\n";
print outfile1 "Number of bases\t$totalbases";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$sampledata{$sample}{bases}"; }
print outfile1 "\n";

	#-- Contigs

print outfile1 "\n# Statistics on contigs\n";
print outfile1 "#\tAssembly\n";
print outfile1 "Number of contigs\t$contigs{num}\nTotal length\t$contigs{bases}\nLongest contig\t$contigs{max}\nShortest contig\t$contigs{min}\n";
print outfile1 "N50\t$contigs{N50}\nN90\t$contigs{N90}\n";
foreach my $rk(@ranks) { 
	my @ctk=keys %{ $contax{$rk} };
	my $nmtax=$#ctk+1;
	print outfile1 "Contigs at $rk rank\t$contigs{$rk}, in $nmtax taxa\n"; 	
	}
print outfile1 "Congruent\t$contigs{chimerism}{0}\nDisparity >0\t$contigs{chimerism}{more0}\nDisparity >= 0.25\t$contigs{chimerism}{0.25}\n";
print outfile1 "\n";

	#-- Genes

print outfile1 "# Statistics on ORFs\n";
print outfile1 "#\tAssembly";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$sample"; }
print outfile1 "\n";
print outfile1 "ORF Number\t$genes{totgenes}";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$genes{$sample}{totgenes}"; }
print outfile1 "\n";
print outfile1 "RNAs\t$genes{rnas}";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$genes{$sample}{rnas}"; }
print outfile1 "\n";
print outfile1 "KEGG annotations\t$genes{kegg}";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$genes{$sample}{kegg}"; }
print outfile1 "\n";
print outfile1 "COG annotations\t$genes{cog}";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$genes{$sample}{cog}"; }
print outfile1 "\n";
print outfile1 "Pfam annotations\t$genes{pfam}";
foreach my $sample(sort keys %sampledata) { print outfile1 "\t$genes{$sample}{pfam}"; }
print outfile1 "\n";

	#-- Bins

if($mode ne "sequential") {
	print outfile1 "\n# Statistics on bins\n#";
	foreach my $method(sort keys %bins) { print outfile1 "\t$method"; }
	print outfile1 "\n";
	print outfile1 "Number of bins";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{num}"; }
	print outfile1 "\n";
	print outfile1 "Complete >= 50%";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{complete}{50}"; }
	print outfile1 "\n";
	print outfile1 "Complete >= 75%";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{complete}{75}"; }
	print outfile1 "\n";
	print outfile1 "Complete >= 90%";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{complete}{90}"; }
	print outfile1 "\n";
	print outfile1 "Contamination < 10%";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{contamination}{10}"; }
	print outfile1 "\n";
	print outfile1 "Contamination >= 50%";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{contamination}{50}"; }
	print outfile1 "\n";
	print outfile1 "Congruent bins";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{chimerism}{0}"; }
	print outfile1 "\n";
	print outfile1 "Disparity >0";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{chimerism}{more0}"; }
	print outfile1 "\n";
	print outfile1 "Disparity >= 0.25";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{chimerism}{0.25}"; }
	print outfile1 "\n";
	print outfile1 "Hi-qual bins (>90% complete,<10% contam)";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{hiqual}"; }
	print outfile1 "\n";
	print outfile1 "Good-qual bins (>75% complete,<10% contam)";
	foreach my $method(sort keys %bins) { print outfile1 "\t$bins{$method}{goodqual}"; }
        print outfile1 "\n";
	}

close outfile1;
print "Output in $resultfile\n";
