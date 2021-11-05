#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Last Common Ancestor (LCA) taxonomic assignment from a Diamond file. 
#-- 20/05/2018 Includes rank-directed identity level checks

$|=1;

use strict;
use DBI;
use DBD::SQLite::Constants qw/:file_open/;
use Tie::IxHash;
use Linux::MemInfo;
use Cwd;
use lib ".";
use threads;

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$databasepath,$interdir,$tempdir,$taxdiamond,$lca_db,$fun3tax,$numthreads,$evalue,$scoreratio6,$diffiden6,$flex6,$minhits6,$noidfilter6,$syslogfile);
my $infile=$taxdiamond;

#-- Set maximum number of threads (we aim for 1 thread per Gb of free RAM)
my %mem=get_mem_info;
my $ram=$mem{"MemAvailable"}/(1024*1024);
my $ramstr=sprintf('%.2f',$ram);
my $ramthreads=int($ram)-2;
if($ramthreads < $numthreads) {
	print "  AVAILABLE (free) RAM memory: $ramstr Gb\n  We will set the number of threads to $ramthreads for this step\n";
	$numthreads = $ramthreads;
	}



#-- Some parameters for the algorithm

my @ranks=('species','genus','family','order','class','phylum','superkingdom');
my %idenrank=('species',85,'genus',60,'family',55,'order',50,'class',46,'phylum',42,'superkingdom',40);
my $verbose=0;
my $thereareresults=0;

my(%provhits,%accum,%accumnofilter,%providen);
my($thisorf,$lastorf,$validhits,$validhitsnofilter,$tothits,$refscore,$refiden,$assigned,%giden);
tie %provhits,"Tie::IxHash";
tie %accum,"Tie::IxHash";
tie %accumnofilter,"Tie::IxHash";

open(syslogfile,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Reads the taxonomic tree (parsed from NCBI's taxonomy in the parents.txt file)

my %parents;
open(infile1,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	$tax=~s/\[|\]//g;
	$parents{$tax}{wranks}=$par;
	my @m=split(/\;/,$par);
	foreach my $y(@m) {
		my($rt,$gtax)=split(/\:/,$y);
		$parents{$tax}{noranks}.="$gtax;"; 
		}
	chop $parents{$tax}{noranks};
	}
close infile1;

#-- Split Diamond file

splitfiles();

#-- Launch threads

print "  Starting multithread LCA in $numthreads threads\n";
print syslogfile "  Starting multithread LCA in $numthreads threads\n";
my $threadnum;
for($threadnum=1; $threadnum<=$numthreads; $threadnum++) {
#print "$threadnum ";
my $thr=threads->create(\&current_thread,$threadnum);
}
#print "\n";
$_->join() for threads->list();

my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.wranks "; }
$catcommand.="> $fun3tax.wranks";
print "  Creating $fun3tax.wranks file\n";
print syslogfile "  Creating $fun3tax.wranks file: $catcommand\n";
my $ecode = system($catcommand);

my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.noidfilter.wranks "; }
$catcommand.=" > $fun3tax.noidfilter.wranks";
print "  Creating $fun3tax.noidfilter.wranks file\n";
print syslogfile "  Creating $fun3tax.noidfilter.wranks file: $catcommand\n";
$ecode = system $catcommand;
print syslogfile "  Removing temporaty diamond files in $tempdir\n";
system("rm $tempdir/diamond_lca.*.m8");
close syslogfile;
 

sub splitfiles {
	print "  Splitting Diamond file\n";
	print syslogfile "  Splitting Diamond file\n";
	system("wc -l $infile > $tempdir/wc");
	open(intemp,"$tempdir/wc");
	my $wc=<intemp>;
	close intemp;
	chomp $wc;
	$wc=~s/\s+.*//;    #-- Number of lines in the diamond result
	my $splitlines=int($wc/$numthreads);
	print syslogfile "  Total lines in Diamond: $wc; Allocating $splitlines in $numthreads threads\n";
	my $nextp=$splitlines;
	my ($filelines,$splitorf);
	my $numfile=1;
	print syslogfile "  Opening file $numfile in line $filelines (estimated in $nextp)\n";
	open(outfiletemp,">$tempdir/diamond_lca.$numfile.m8");
	open(infile2,$infile) || die "Can't open Diamond file $infile\n"; 
	while(<infile2>) {
		$filelines++;
		my @f=split(/\t/,$_);
		if($filelines==$nextp) { $splitorf=$f[0];  }
		if($filelines<=$nextp) { print outfiletemp $_; next; }                      
		elsif($f[0] ne $splitorf) { 
			close outfiletemp;
			$numfile++;
			print syslogfile "  Opening file $numfile in line $filelines (estimated in $nextp)\n"; 
			open(outfiletemp,">$tempdir/diamond_lca.$numfile.m8");
			print outfiletemp $_;
			$nextp+=$splitlines;
			}
                else { print outfiletemp $_; }
		}
	close infile2;
	}


sub current_thread {

#-- Preparing the output files

my $threadnum=shift;
print syslogfile "  Starting thread $threadnum\n";

#-- Prepare the LCA database (containing the acc -> tax correspondence)

my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1, sqlite_open_flags => SQLITE_OPEN_READONLY }) or die $DBI::errstr;
$dbh->sqlite_busy_timeout( 120 * 1000 );

my $currentfile="$tempdir/diamond_lca.$threadnum.m8";
open(outfile2,">$tempdir/fun3tax\_$threadnum.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.wranks for writing\n";
if($noidfilter6) {
        open(outfile4,">$tempdir/fun3tax\_$threadnum.noidfilter.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.noidfilter.wranks for writing\n";
        }
open(infile2,$currentfile) || die "Cannot open $currentfile\n"; 
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
		
		query($dbh);
   
		(%accum,%accumnofilter,%provhits,%providen,%giden)=();
 		($validhits,$validhitsnofilter,$tothits)=0;
		$string="";
		$lastorf=$thisorf;	
		($refscore,$refiden)=0;	
		$assigned++;
	#	if(!($assigned%1000)) { print "  $assigned ORFs processed     \r"; }
		}

	#-- If we are reading a hit, we store its bitscore and identity value
	#-- If it is the first (and therefore best) hit, we take that values as reference	

	if(!$refscore) { $refscore=$fields[$#fields]; }
	if(!$refiden) { $refiden=$fields[2]; }  			   
	if($fields[$#fields]>$provhits{$fields[1]}) { $provhits{$fields[1]}=$fields[$#fields]; }
  	if($fields[2]>$providen{$fields[1]}) { $providen{$fields[1]}=$fields[2]; }
        my @hitfields=split(/\|/,$fields[1]);
        $giden{$fields[1]}=$fields[2];
	# print "PROVHIT: $fields[1] $provhits{$fields[1]} $providen{$fields[1]}\n";
	$tothits++;			   
       }
close infile2;

#-- For the last ORF, we have to call again query() because the file ended before we could make the call

$lastorf=$thisorf;
query($dbh);    

# close outfile1;
close outfile2;
#close outfile3;
close outfile4;
# print "  Total: $assigned ORFs processed\n";

if(!$thereareresults) { die "Tax assignment done in $fun3tax\_$threadnum.wranks but no results found. Aborting\n"; }

print "  Tax assignment done! Result stored in file $fun3tax\_$threadnum.wranks\n" if $verbose;
}

sub query {
	my $dbh=shift;
	my($refcc,$genocc,$unicc,$ratioscore,$idendiff)=0; 
	# if($lastorf=~/k141_440_IC1025022_4/) { $verbose=1; } else { $verbose=0; }
	print "refscore: $refscore refiden: $refiden\n" if $verbose; 
	
	#-- We start building the query to the database, using the acc numbers
	
	my $query="select * from taxid where (";
	foreach my $lhits(keys %provhits) {

		#-- We evaluate if the hit is valid (bitscore and identity are within the limits)

		if($refscore) { $ratioscore=$provhits{$lhits}/$refscore; }
    		next if($ratioscore<=$scoreratio6);
		if($refiden) { $idendiff=$refiden-$providen{$lhits}; }
    		next if($idendiff>$diffiden6);    

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
				if($giden{$list[0]}>=$idenrank{$rank}) { $accum{$rank}{$tax}++;  }		#-- and add a count for that taxon in that rank
                                $accumnofilter{$rank}{$tax}++; 		#-- Not considering identity filters for ranks
			}
		if(($list[7]) && ($giden{$list[0]}>=$idenrank{'superkingdom'}))  { $validhits++; }			#-- Count the number of valid hits
                if(($list[7])) { $validhitsnofilter++; }
		}
	}

	#-- Now, if there are some results, we will find the LCA

	my($minreqhits,$required);
	# $validhitsnofilter=$validhits;
	if($validhits==1) { $minreqhits=1; } else { $minreqhits=$minhits6; }
	if($flex6<1) { $required=$validhits-($flex6*$validhits); } else { $required=$validhits-$flex6; }  
	print "$lastorf Hits: $tothits; Valid: $validhits; Min: $minreqhits; Required: $required\n" if $verbose;
	my($lasttax,$lasttaxnofilter)="";
		
	#-- Looping on the ranks from top to bottom, looking if it fulfills the conditions for LCA
	
	foreach my $k(@ranks) {	
		print "   $k\n" if $verbose;
		my $maxp=0;
		foreach my $t(sort { $accum{$k}{$a}<=>$accum{$k}{$b}; } keys %{ $accum{$k} }) {
			if($t && ($accum{$k}{$t}==$maxp)) { $lasttax=""; next; }    #-- Equality of hits, don´t choose any
			if($t) { $maxp=$accum{$k}{$t}; }
			print "      $t $accum{$k}{$t}\n" if $verbose;
			if(($accum{$k}{$t}>=$required) && ($accum{$k}{$t}>=$minreqhits) && ($required>0) && ($parents{$t}{wranks})) { 	#-- REQUIREMENTS FOR VALID LCA
				print "$k -> $t\n" if $verbose;
				$lasttax=$t; 
				#  if($t) { $string="$t;$string"; }
				}
			}

		last if($lasttax);		
 		}
		
	if($noidfilter6) {
		if($validhitsnofilter==1) { $minreqhits=1; } else { $minreqhits=$minhits6; }
		if($flex6<1) { $required=$validhitsnofilter-($flex6*$validhitsnofilter); } else { $required=$validhitsnofilter-$flex6; }
		print "NOFILTER $lastorf Hits: $tothits; Valid: $validhits; Min: $minreqhits; Required: $required\n" if $verbose;
		$lasttaxnofilter="";			
		foreach my $k(@ranks) {
			print "   NOFILTER $k\n" if $verbose;
			foreach my $t(keys %{ $accumnofilter{$k} }) {
			print "      NOFILTER $t $accumnofilter{$k}{$t}\n" if $verbose;
				if(($accumnofilter{$k}{$t}>=$required) && ($accumnofilter{$k}{$t}>=$minreqhits) && ($required>0) && ($parents{$t}{wranks})) { $lasttaxnofilter=$t; }
				print "NOFILTER $k -> $t\n" if $verbose;
				}

			last if($lasttaxnofilter);
			}
		}		
		
	# print outfile1 "$lastorf\t$parents{$lasttax}{noranks}\n";
	my $abb=$parents{$lasttax}{wranks};
	
	#-- Changing nomenclature to abbreviations
	
	$abb=~s/sub\w+\:/n_/g;
	$abb=~s/sub\w+\_/n_/g;
	$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
	# print outfile2 "$lastorf\t$parents{$lasttax}{wranks}\n";		
	print outfile2 "$lastorf\t$abb\n";		
	if($noidfilter6) {
		# print outfile3 "$lastorf\t$parents{$lasttaxnofilter}{noranks}\n";
		my $abb=$parents{$lasttaxnofilter}{wranks};
		$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g; 
		# print outfile4 "$lastorf\t$parents{$lasttaxnofilter}{wranks}\n";	
		print outfile4 "$lastorf\t$abb\n";	
		}	
	 print "$lastorf\t$parents{$lasttax}{noranks}\n" if $verbose;
	#print "$lastorf\t$abb\n" if $verbose;
	if($parents{$lasttax}{noranks}) { $thereareresults=1; }	
}

