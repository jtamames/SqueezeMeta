#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Counts the coverage of all functions

use strict;
use Tie::IxHash;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
my $taxreq=$ARGV[1];	#-- Invoke it with a name of taxon to get just functions for that taxon

do "$project/squeezeM_conf.pl";


#-- Configuration variables from conf file

our($datapath,$resultpath,$kegglist,$coglist,$ntfile,$fun3tax,$fun3kegg,$fun3cog,$coveragefile,$rpkmfile);

print "Calculating coverage for functions\n";
my(%funs,%taxf,%validid,%tfun,%totalbases,%totalreads,%allsamples,%funstat,%longorfs,%taxcount);

	#-- Reading KEGG functions and pathways

open(infile1,$kegglist) || warn "Missing KEGG equivalence file $kegglist\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$funs{kegg}{$t[0]}{name}=$t[1];
	$funs{kegg}{$t[0]}{fun}=$t[2];
	$funs{kegg}{$t[0]}{path}=$t[3];
	}
close infile1;

	#-- Reading COG functions and pathways

open(infile2,$coglist) || warn "Missing COG equivalence file $coglist\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$funs{cog}{$t[0]}{name}=$t[1];
	$funs{cog}{$t[0]}{fun}=$t[2];
	$funs{cog}{$t[0]}{path}=$t[3];
	}
close infile2;

my %equival;
tie %equival,"Tie::IxHash";
%equival=('superkingdom','k','phylum','p','class','c','order','o','family','f','genus','g','species','s');

	#-- Read the taxonomic assignment for each gene

open(infile3,"$fun3tax.wranks") || warn "Cannot open wranks file $fun3tax.wranks\n";
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my @t=split(/\;/,$k[1]);
	
	#-- Loop for all ranks for the gene
	
	foreach my $cm(@t) {
		my ($rank,$ttax)=split(/\:/,$cm);
		my $cortax=$equival{$rank};
		next if(!$cortax);
		$taxf{$k[0]}{$cortax}=$ttax;
		if($taxreq && ($ttax eq $taxreq)) { $validid{$k[0]}=1; }
		# print "$k[0] $cortax $ttax\n";
		}
	}
close infile3;

	#-- Read sequence length from fasta file (old version, now it can be obtained from the coverage file)

#open(in,$ntfile) || die "I need the DNA sequences from the prediction\n";
#while(<in>) {	
#	chomp;
#	if($_=~/^\>([^ ]+)/) { 
#		if($seq) {  $longorfs{$thisorf}=(length $seq)+1; }
#		$thisorf=$1;
#		$seq="";
#		}
#	else { $seq.=$_; }		     
#	}
#close in;
#if($seq) { $allorfs{$thisorf}=$seq; }

	#-- Reading KEGG assignments

open(infile4,$fun3kegg)  || die;
while(<infile4>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	$tfun{$k[0]}{kegg}=$k[1];
	}
close infile4;

	#-- Reading COG assignments

open(infile5,$fun3cog)  || die;
while(<infile5>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	$tfun{$k[0]}{cog}=$k[1];
	}
close infile5;

	#-- Reading coverages for all genes

print "Reading coverage in $coveragefile\n";
open(infile6,$coveragefile) || die;
while(<infile6>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $cfun_kegg=$tfun{$k[0]}{kegg};	#-- Corresponding KEGG for this gene
	my $cfun_cog=$tfun{$k[0]}{cog}; 	#-- Corresponding COG for this gene
	my $sample=$k[3];
	$longorfs{$k[0]}=$k[4];		#-- Length of the gene
	# print "$k[0] $cfun_kegg $cfun_cog $sample $longorfs{$k[0]}\n";
	$totalbases{$sample}+=$k[1];
	next if((!$cfun_kegg) && (!$cfun_cog));
	next if($taxreq && (!$validid{$k[0]}));
	$allsamples{$sample}++;
	if($k[1]) {
	
		#-- Counting KEGGs
	
		if($cfun_kegg) {
		$funstat{kegg}{$cfun_kegg}{$sample}{copies}++;
		# $funstat{$cfun_kegg}{$sample}{length}+=$k[4];
		$funstat{kegg}{$cfun_kegg}{$sample}{length}+=$longorfs{$k[0]}; 
		$funstat{kegg}{$cfun_kegg}{$sample}{bases}+=$k[1];
		foreach my $tk(keys %equival) {
			my $krank=$equival{$tk};
			my $itax=$taxf{$k[0]}{$krank};
			if($itax) { $taxcount{kegg}{$cfun_kegg}{$sample}{$krank}{$itax}++; }
			}
		}
	
		#-- Counting COGs
	
	if($cfun_cog) { 
		$funstat{cog}{$cfun_cog}{$sample}{copies}++;
		# $funstat{$cfun_cog}{$sample}{length}+=$k[4];
		$funstat{cog}{$cfun_cog}{$sample}{length}+=$longorfs{$k[0]}; #-- Para leer las longitudes directamente de las secuencias (para ficheros coverage antiguos que no la tienen)
		$funstat{cog}{$cfun_cog}{$sample}{bases}+=$k[1];
		# print "$k[0]*$cfun_cog*$sample*$longorfs{$k[0]}*$funstat{cog}{$cfun_cog}{$sample}{length}\n";
		foreach my $tk(keys %equival) {
			my $krank=$equival{$tk};
			my $itax=$taxf{$k[0]}{$krank};
			if($itax) { $taxcount{cog}{$cfun_cog}{$sample}{$krank}{$itax}++; }
 			}
		}
	}
}
close infile6;


	#-- Reading RPKMs for all genes

print "Reading rpkm in $rpkmfile\n";
open(infile7,$rpkmfile) || die;
while(<infile7>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $cfun_kegg=$tfun{$k[0]}{kegg}; 
	my $cfun_cog=$tfun{$k[0]}{cog}; 
	my $sample=$k[3];
	$totalreads{$sample}+=$k[2];
	next if((!$cfun_kegg) && (!$cfun_cog));
	next if($taxreq && (!$validid{$k[0]}));
	if($k[2]) {
		if($cfun_kegg) { $funstat{kegg}{$cfun_kegg}{$sample}{reads}+=$k[2]; }
		if($cfun_cog) { $funstat{cog}{$cfun_cog}{$sample}{reads}+=$k[2]; }  
		}
	}
close infile7;

	#-- Creating output files

foreach my $classfun(sort keys %funstat) {
	print "Now creating $classfun coverage output\n";
	open(outfile1,">$resultpath/11.$project.$classfun.funcover") || die;
	print outfile1 "#-- Created by $0 from $coveragefile, ",scalar localtime;
	if($taxreq) { print outfile1 ", for taxon $taxreq"; }
	print outfile1 "\n";
	print outfile1 "# $classfun ID\tSample\tCopy number\tTotal length\tTotal bases\tCoverage\tNorm Coverage\tDistribution\tName\tFunction\n";
	
	#-- Calculation of coverage, norm coverage, and RPKM

	foreach my $samp(sort keys %allsamples) {
		foreach my $kid(sort keys %{ $funstat{$classfun} }) {
			next if(!$funstat{$classfun}{$kid}{$samp}{length}); 
			my $cover=$funstat{$classfun}{$kid}{$samp}{bases}/$funstat{$classfun}{$kid}{$samp}{length};
			my $normcover=(($funstat{$classfun}{$kid}{$samp}{bases}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalbases{$samp}));
			my $rpkm=(($funstat{$classfun}{$kid}{$samp}{reads}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalreads{$samp}));
 
			my $stringtax=""; 
			foreach my $tk(keys %equival) {
				my $krank=$equival{$tk};
				my $countt=0;
				foreach my $tt(keys %{ $taxcount{$classfun}{$kid}{$samp}{$krank} }) { $countt++; }
				$stringtax.="$krank:$countt;";
				}
			chop $stringtax;
			printf outfile1 "$kid\t$samp\t$funstat{$classfun}{$kid}{$samp}{copies}\t$funstat{$classfun}{$kid}{$samp}{length}\t$funstat{$classfun}{$kid}{$samp}{bases}\t%.3f\t%.3f\t$stringtax\t$funs{$classfun}{$kid}{name}\t$funs{$classfun}{$kid}{fun}\n",$cover,$normcover; 
 			}
		}
close outfile1;
	}
