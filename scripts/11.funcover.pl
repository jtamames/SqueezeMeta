#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Counts the coverage of all functions

use strict;
use Tie::IxHash;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
my $taxreq=$ARGV[1];	#-- Invoke it with a name of taxon to get just functions for that taxon

$project=~s/\/$//;
do "$project/SqueezeMeta_conf.pl";


#-- Configuration variables from conf file

our($datapath,$resultpath,$kegglist,$coglist,$ntfile,$fun3tax,$fun3kegg,$fun3cog,$nokegg,$nocog,$coveragefile,$rpkmfile);

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
	$funs{cog}{$t[0]}{fun}=$t[1];
	$funs{cog}{$t[0]}{path}=$t[2];
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

if(!$nokegg) {
	open(infile4,$fun3kegg) || die "Cannot open KEGG assignments in $fun3kegg\n";;
	while(<infile4>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		$tfun{$k[0]}{kegg}=$k[1];
		}
	close infile4;
	}
	
	#-- Reading COG assignments

if(!$nocog) {
	open(infile5,$fun3cog) || die "Cannot open COG assignments in $fun3cog\n";
	while(<infile5>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		$tfun{$k[0]}{cog}=$k[1];
		}
	close infile5;
	}
	
	#-- Reading coverages for all genes

print "Reading coverage in $coveragefile\n";
open(infile6,$coveragefile) || warn "Cannot open coverages in $coveragefile\n";
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
			my @kegglist=split(/\;/,$cfun_kegg);	#-- Support for multiple COGS (in annotations such as COG0001;COG0002, all COGs get the counts)
			foreach my $tlist_kegg(@kegglist) {
				$funstat{kegg}{$tlist_kegg}{$sample}{copies}++;
				$funstat{kegg}{$tlist_kegg}{$sample}{length}+=$longorfs{$k[0]}; 
				$funstat{kegg}{$tlist_kegg}{$sample}{bases}+=$k[1];
				foreach my $tk(keys %equival) {
					my $krank=$equival{$tk};
					my $itax=$taxf{$k[0]}{$krank};
					if($itax) { $taxcount{kegg}{$tlist_kegg}{$sample}{$krank}{$itax}++; }
					}
				}
			}
	
		#-- Counting COGs
	
	if($cfun_cog) { 
		my @coglist=split(/\;/,$cfun_cog);	#-- Support for multiple COGS (in annotations such as COG0001;COG0002, all COGs get the counts)
		foreach my $tlist_cog(@coglist) {
			$funstat{cog}{$tlist_cog}{$sample}{copies}++;
			$funstat{cog}{$tlist_cog}{$sample}{length}+=$longorfs{$k[0]}; #-- Para leer las longitudes directamente de las secuencias (para ficheros coverage antiguos que no la tienen)
			$funstat{cog}{$tlist_cog}{$sample}{bases}+=$k[1];
			# print "$k[0]*$cfun_cog*$sample*$longorfs{$k[0]}*$funstat{cog}{$cfun_cog}{$sample}{length}\n";
			foreach my $tk(keys %equival) {
				my $krank=$equival{$tk};
				my $itax=$taxf{$k[0]}{$krank};
				if($itax) { $taxcount{cog}{$tlist_cog}{$sample}{$krank}{$itax}++; }
 				}
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
		my @kegglist=split(/\;/,$cfun_kegg);
		my @coglist=split(/\;/,$cfun_cog);
		foreach my $tlist_kegg(@kegglist) { 
			if($tlist_kegg) { $funstat{kegg}{$tlist_kegg}{$sample}{reads}+=$k[2]; }
			}
		foreach my $tlist_cog(@coglist) { 
			if($tlist_cog) { $funstat{cog}{$tlist_cog}{$sample}{reads}+=$k[2]; }  
			}
		}
	}
close infile7;

	#-- Creating output files

my $rawf;
foreach my $classfun(sort keys %funstat) {
	$rawf="$resultpath/11.$project.$classfun.funcover";
	print "Now creating $classfun coverage output in $rawf\n";
	open(outfile1,">$rawf") || die;
	print outfile1 "#-- Created by $0 from $coveragefile, ",scalar localtime;
	if($taxreq) { print outfile1 ", for taxon $taxreq"; }
	print outfile1 "\n";
	print outfile1 "# $classfun ID\tSample\tCopy number\tTotal length\tTotal bases\tCoverage\tNorm Coverage\tNorm Coverage per copy\tDistribution\tName\tFunction\n";
	
	#-- Calculation of coverage, norm coverage, and RPKM

	foreach my $samp(sort keys %allsamples) {
		foreach my $kid(sort keys %{ $funstat{$classfun} }) {
			next if(!$funstat{$classfun}{$kid}{$samp}{length}); 
			my $cover=$funstat{$classfun}{$kid}{$samp}{bases}/$funstat{$classfun}{$kid}{$samp}{length};
			my $normcover=(($funstat{$classfun}{$kid}{$samp}{bases}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalbases{$samp}));
			my $normcoverpercopy=$normcover*$funstat{$classfun}{$kid}{$samp}{copies};
			# print "$classfun*$kid*$samp*$funstat{$classfun}{$kid}{$samp}{length}*$totalreads{$samp}\n";
			my $rpkm=(($funstat{$classfun}{$kid}{$samp}{reads}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalreads{$samp}));
 
			my $stringtax=""; 
			foreach my $tk(keys %equival) {
				my $krank=$equival{$tk};
				my $countt=0;
				foreach my $tt(keys %{ $taxcount{$classfun}{$kid}{$samp}{$krank} }) { $countt++; }
				$stringtax.="$krank:$countt;";
				}
			chop $stringtax;
			printf outfile1 "$kid\t$samp\t$funstat{$classfun}{$kid}{$samp}{copies}\t$funstat{$classfun}{$kid}{$samp}{length}\t$funstat{$classfun}{$kid}{$samp}{bases}\t%.3f\t%.3f\t%.3f\t$stringtax\t$funs{$classfun}{$kid}{name}\t$funs{$classfun}{$kid}{fun}\n",$cover,$normcover,$normcoverpercopy; 
 			}
		}
close outfile1;
	}

my $minraw=200;
foreach my $classfun(sort keys %funstat) {
	$rawf="$resultpath/11.$project.$classfun.rawcounts";
	print "Now creating $classfun raw reads output in $rawf\n";
	open(outfile2,">$rawf") || die;
	if($classfun eq "cog") { print outfile2 "$classfun class\t$classfun ID"; }
        else { print outfile2 "$classfun ID"; }
foreach my $samp(sort keys %allsamples) { print outfile2 "\t$samp"; }
	print outfile2 "\n";
	foreach my $kid(sort keys %{ $funstat{$classfun} }) {
                my($cstring,$accum,$funid,$funclass);
                $funclass=$funs{$classfun}{$kid}{path};
                if(!$funclass || ($funclass=~/\|/)) { $funclass="Unclassified"; } 
                if($classfun eq "cog") { $funid="$funclass\t"; }
                $funid.="$kid:$funs{$classfun}{$kid}{fun}";
		# print outfile2 "$funid";
		$cstring.="$funid";		
		foreach my $samp(sort keys %allsamples) { 
			my $rnum=$funstat{$classfun}{$kid}{$samp}{reads} || "0";
			$accum+=$rnum;
			$cstring.="\t$rnum"; 
			#print outfile2 "\t$funstat{$classfun}{$kid}{$samp}{reads}"; 
			}
		if($accum>=$minraw) { print outfile2 "$cstring\n"; }
		# print outfile2 "\n";
		}
	close outfile2;
	}


	


	
