#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Last Common Ancestor (LCA) taxonomic assignment from a Diamond file. 
#-- 20/05/2018 Includes rank-directed identity level checks

$|=1;

use strict;
use DBI;
use Tie::IxHash;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$databasepath,$taxdiamond,$lca_db,$fun3tax,$evalue);
my $infile=$taxdiamond;

#-- Some parameters for the algorithm

my @ranks=('species','genus','family','order','class','phylum','superkingdom');
my %idenrank=('species',85,'genus',60,'family',55,'order',50,'class',46,'phylum',42,'superkingdom',40);
my $scoreratio=0.8;   #-- Ratio first score/currsent score for the hit to be considered
my $diffiden=10;       #-- Maximim identity difference with the first
my $flex=0.2;           #-- Allows this PERCENTAGE (if less than one) or NUMBER (if greater than one) of hits from different taxa than LCA
my $minhits=2;        #-- Minimum number of hits for the taxa (if there is only one valid hit, this value sets to one automatically
my $verbose=0;
my $thereareresults=0;

#-- Prepare the LCA database (containing the acc -> tax correspondence)

my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1}) or die $DBI::errstr;

#-- Reads the taxonomic tree (parsed from NCBI's taxonomy in the parents.txt file)

my %parents;
open(infile1,"$databasepath/LCA_tax/parents.txt") || die;
while(<infile1>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	$parents{$tax}{wranks}=$par;
	my @m=split(/\;/,$par);
	foreach my $y(@m) {
		my($rt,$gtax)=split(/\:/,$y);
		$parents{$tax}{noranks}.="$gtax;"; 
		}
	chop $parents{$tax}{noranks};
	}
close infile1;

#-- Preparing the output files

open(outfile1,">$fun3tax") || die;
print outfile1 "# Created by $0 from $infile, ",scalar localtime,", evalue=$evalue, scoreratio=$scoreratio, diffiden=$diffiden, flex=$flex, minhits=$minhits\n";
open(outfile2,">$fun3tax.wranks") || die;
print outfile2 "# Created by $0 from $infile, ",scalar localtime,", evalue=$evalue, scoreratio=$scoreratio, diffiden=$diffiden, flex=$flex, minhits=$minhits\n";

#-- Parsing of the diamond file

my(%provhits,%accum,%providen);
my($thisorf,$lastorf,$validhits,$tothits,$refscore,$refiden,%giden);
tie %provhits,"Tie::IxHash";
tie %accum,"Tie::IxHash";

if($infile=~/gz/) { open(infile2,"zcat $infile|") || die; }			#-- If file is gzipped
else { open(infile2,$infile) || die "Cannot open Diamond file $infile\n"; }	#-- or if it is not

while(<infile2>) { 
	chomp;
	my $string;
	next if(!$_ || ($_=~/^\#/));	
	$_=~s/\;\_//g;
	my @fields=split(/\t/,$_);
	$thisorf=$fields[0];
	if(!$lastorf) { $lastorf=$thisorf; }
	if($lastorf && ($thisorf ne $lastorf)) {	

		#-- We finished reading the hits for an ORF,then we will query the database for the taxa of the hits
		#-- and then we clean all the data corresponding to the past ORF 
		
		query();
   
		(%accum,%provhits,%providen,%giden)=();
 		($validhits,$tothits)=0;
		$string="";
		$lastorf=$thisorf;	
		($refscore,$refiden)=0;	
		}

	#-- If we are reading a hit, we store its bitscore and identity value
	#-- If it is the first (and therefore best) hit, we take that values as reference	

	if(!$refscore) { $refscore=$fields[$#fields]; }
	if(!$refiden) { $refiden=$fields[2]; }  			   
	$provhits{$fields[1]}=$fields[$#fields];
	$providen{$fields[1]}=$fields[2];
        my @hitfields=split(/\|/,$fields[1]);
        $giden{$fields[1]}=$fields[2];
	# print "PROVHIT: $fields[1] $provhits{$fields[1]} $providen{$fields[1]}\n";
	$tothits++;			   
       }
close infile2;

#-- For the last ORF, we have to call again query() because the file ended before we could make the call

$lastorf=$thisorf;
query();    

close outfile1;
close outfile2;

if(!$thereareresults) { die "Tax assignment done in $fun3tax but no results found. Aborting\n"; }

print "Tax assignment done! Result stored in file $fun3tax\n";


sub query {
	my($refcc,$genocc,$unicc,$ratioscore,$idendiff)=0; 
	# if($lastorf=~/k141_440_IC1025022_4/) { $verbose=1; } else { $verbose=0; }
	print "refscore: $refscore refiden: $refiden\n" if $verbose; 
	
	#-- We start building the query to the database, using the acc numbers
	
	my $query="select * from taxid where (";
	foreach my $lhits(keys %provhits) {

		#-- We evaluate if the hit is valid (bitscore and identity are within the limits)

		if($refscore) { $ratioscore=$provhits{$lhits}/$refscore; }
    		next if($ratioscore<=$scoreratio);
		if($refiden) { $idendiff=$refiden-$providen{$lhits}; }
    		next if($idendiff>$diffiden);    

		#-- If it is, we add its acc to the query

    		my @hitfields=split(/\|/,$lhits);
      		if($refcc) { $query.=" or "; }
      		else { $refcc=1; } 
      		$query.="id=\"$hitfields[0]\"";
		}

	$query.=");";	
	print "*$query*\n" if $verbose;			     

	#-- If the query is not empty, we make the call to the database	

	if($refcc) {
		my $sth = $dbh->prepare($query);  
		$sth->execute();

		#-- And start reading the result

		my @list;
		while(@list=$sth->fetchrow()) {
			print "$lastorf\t@list\n" if $verbose;
			for(my $pos=1; $pos<=7; $pos++) {	#-- We loop on the taxonomic ranks for the hit  #-- Fixed bug when updating to new nr db. Missing field makes parsing incorrectly the taxonomy. JT 14/11/18
				my $rank=$ranks[$pos-1];	#-- Retrieve the rank
				my $tax=$list[$pos];		#-- and the taxon
				print "$lastorf $rank $giden{$list[0]} $idenrank{$rank}\n" if $verbose;
				if($giden{$list[0]}>=$idenrank{$rank}) { $accum{$rank}{$tax}++; }		#-- and add a count for that taxon in that rank
			}
		if(($list[7]) && ($giden{$list[0]}>=$idenrank{'superkingdom'})) { $validhits++;  }			#-- Count the number of valid hits
		}
	}

	#-- Now, if there are some results, we will find the LCA

	my($minreqhits,$required);
	if($validhits==1) { $minreqhits=1; } else { $minreqhits=$minhits; }
	if($flex<1) { $required=$validhits-($flex*$validhits); } else { $required=$validhits-$flex; }  
	print "$lastorf Hits: $tothits; Valid: $validhits; Min: $minreqhits; Required: $required\n" if $verbose;
	my $lasttax="";
		
	#-- Looping on the ranks from top to bottom, looking if it fulfills the conditions for LCA
	
	foreach my $k(@ranks) {	
		print "   $k\n" if $verbose;
		foreach my $t(keys %{ $accum{$k} }) {
			print "      $t $accum{$k}{$t}\n" if $verbose;
			if(($accum{$k}{$t}>=$required) && ($accum{$k}{$t}>=$minreqhits) && ($required>0)) { 	#-- REQUIREMENTS FOR VALID LCA
				print "$k -> $t\n" if $verbose;
				$lasttax=$t; 
				#  if($t) { $string="$t;$string"; }
				}
			}

		last if($lasttax);		
 		}
	print outfile1 "$lastorf\t$parents{$lasttax}{noranks}\n";
	print outfile2 "$lastorf\t$parents{$lasttax}{wranks}\n";		
	print "$lastorf\t$parents{$lasttax}{noranks}\n" if $verbose;
	if($parents{$lasttax}{noranks}) { $thereareresults=1; }	
}

