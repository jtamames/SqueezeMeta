#!/usr/bin/env perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Needed to create LCA database. Creates a list of species and their corresponding taxa,
#-- from taxonomy (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy)


use strict;

$|=1;

my $lca_dir=$ARGV[0];
die if(!$lca_dir);

	#-- Output files

my $outfile="$lca_dir/taxatree.txt";
my $parentfile="$lca_dir/parents.txt";

open(outfile1,">$outfile");

        #-- Download NCBI files.

my $command="wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz -nc -P $lca_dir";
system $command;
system("tar -xvzf $lca_dir/new_taxdump.tar.gz -C $lca_dir");

	#-- Read taxonomy structure

my(%ranks,%parents,%virus);
open(infile1,"$lca_dir/nodes.dmp") || die;
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t\|\t/,$_);
	my $tid=$t[0]; 
	
	#-- Don't consider some taxonomy divisions (Synthetic and Chimeric,Unassigned,Environmental samples)  
	
	next if($t[4]=~/7|8|11/);
	$ranks{$tid}=$t[2]; 
	$parents{$tid}=$t[1]; 
	
	#-- Store viruses separately (for being able to take the full name of them later)
	
	if($t[4]=~/3|9/) { $virus{$tid}=1; }
	}
close infile1;

	#-- Read taxonomy names

my %names;
open(infile2,"$lca_dir/names.dmp") || die;
while(<infile2>) {
	chomp;
	next if !$_;
	my @t=split(/\t\|\t/,$_);
	next if($t[3]!~/scientific/);	#-- Only store scientific names, not synonyms
	$names{$t[0]}=$t[1];
	}
close infile2;

my %seen;
foreach my $print(sort keys %ranks) {
	my $specname;
        if(!$names{$print}){print "MISSING NAME FOR $print!!\n";}
	#print "*****$print\t$names{$print}\t$ranks{$print}\n";
	next if($ranks{$print} ne "species");
	my @k=split(/\s+/,$names{$print});
	if($virus{$print}) { $specname=$names{$print}; }
	elsif($k[0] eq "Candidatus") { $specname="$k[0] $k[1] $k[2]"; }
	elsif( ($k[1] eq "symbiotica") or ($k[1] eq "symbioticus") or ($k[1] eq "asymbiotica") or ($k[1] eq "asymbioticus") ) { $specname="$k[0] $k[1]"; }
	elsif( $names{$print}=~/symbio/i )  { $specname=$names{$print}; }
	else { $specname="$k[0] $k[1]"; }
	next if($seen{$specname});
	$seen{$specname}=1;
	my $found=1;
	my $string="$ranks{$print}:$print:$specname";
        #print "$string\n";
	my $current=$print;
	
	#-- Travel the taxonomy upwards while we find parents
	
	while($found) {
		my $parent=$parents{$current};
		if($parent) {
		$string.=";$ranks{$parent}:$parent:$names{$parent}";
		$current=$parent;
		}
		
	#-- Cases in which we don't find parents
		
		else { $found=0; }
		if($ranks{$parent} eq "superkingdom") { $found=0; }
		if($parent==1) { $found=0; }
	}
 
	print outfile1 "$string\n";
}

close outfile1;
print "Output created in $outfile\n";
   
	#-- Write a file of parents for traveling taxonomy more easily

my %yseen;
open(outfile2,">$parentfile") || die;
open(infile3,$outfile) || die "Cannot open file $outfile";
while(<infile3>) {
	chomp;
	next if !$_;
	my @k=split(/\;/,$_);
	my $string;
	for(my $pos1=0; $pos1<=$#k; $pos1++) {
		my($root,$string,$idroot)="";
		for(my $pos2=$pos1; $pos2<=$#k; $pos2++) { 
			my @n=split(/\:/,$k[$pos2]);
			if(!$n[0]) { $n[0]="no rank"; }
			$string="$n[0]:$n[2];".$string; 
			if(!$root) { $root="$n[2]"; $idroot=$n[1]; }
			}
		chop $string;	
		$string=~s/ \<prokaryotes\>//;					   
		if(!$yseen{$root}) { print outfile2 "$root\t$string\t$idroot\n"; }
		$yseen{$root}=1;					   
		}
	}
	
close infile3;
close outfile2;
print "Output created in $parentfile\n";
