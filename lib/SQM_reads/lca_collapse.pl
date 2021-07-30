#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 29/01/2019 Original version, (c) Javier Tamames, CNB-CSIC
#-- Last Common Ancestor (LCA) taxonomic assignment from a Diamond file. For blastx collapsed format 

$|=1;

use strict;
use DBI;
use DBD::SQLite::Constants qw/:file_open/;
use Tie::IxHash;
use Cwd;
use lib ".";
use threads;


my $pwd=cwd();
my($infile,$outdir,$scriptdir,$thisfile,$numthreads)=@ARGV;

do "$scriptdir/SqueezeMeta_conf.pl";
do "$scriptdir/parameters.pl";


our($datapath,$resultpath,$databasepath,$tempdir,$syslogfile,$taxdiamond,$lca_db,$fun3tax,$fun3tax_blastx,$evalue,$scoreratio6,$diffiden6,$flex6,$minhits6,$noidfilter6);

$resultpath=$outdir;

my @ranks=('species','genus','family','order','class','phylum','superkingdom');
my %idenrank=('species',85,'genus',60,'family',55,'order',50,'class',46,'phylum',42,'superkingdom',40);
my $noidentical=0;  #-- Drops the first 100% identical hit (for mock)
my $verbose=0;
my $bhitforced=0;	#-- Forces that assignment cannot differ from best hit

my(%provhits,%providen,%giden,%accum,%accumnofilter);
my($validhits,$validhitsnofilter,$tothits,$skipidentical,$refscore,$refiden,$string,$posinit,$posend,$syslogfile);
tie %provhits,"Tie::IxHash";
tie %accum,"Tie::IxHash";
tie %accumnofilter,"Tie::IxHash";

$syslogfile="$outdir/syslog";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Reads the taxonomic tree (parsed from NCBI's taxonomy in the parents.txt file)

my(%parents);
open(infile1,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my($tax,$par)=split(/\t/,$_);
	$parents{$tax}{wranks}=$par;
	my @m=split(/\;/,$par);
	foreach my $y(@m) {
		my($rt,$gtax)=split(/\:/,$y);
		$parents{$tax}{noranks}.="$gtax;"; 
		}
	chop $parents{$tax}{noranks};
	}
close infile1;
$tempdir="$outdir/temp";
system("mkdir $tempdir");

#-- Split Diamond file

splitfiles();

#-- Launch threads

print "  Starting multithread LCA in $numthreads threads: ";
print syslogfile "  Starting multithread LCA in $numthreads threads\n";
my $threadnum;
for($threadnum=1; $threadnum<=$numthreads; $threadnum++) {
	print "$threadnum ";
	my $thr=threads->create(\&current_thread,$threadnum);
}
print "\n";
$_->join() for threads->list();

my $wrankfile="$resultpath/$thisfile.fun3.blastx.tax.wranks";
my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.wranks "; }
$catcommand.=" > $wrankfile";
print "  Creating $wrankfile file\n";
print syslogfile "  Creating $wrankfile file: $catcommand\n";
system $catcommand;

my $wrankfile="$resultpath/$thisfile.fun3.blastx.tax_nofilter.wranks";
my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.nofilter.wranks "; }
$catcommand.=" > $wrankfile";
print "  Creating $wrankfile file\n";
print syslogfile "  Creating $wrankfile file: $catcommand\n";
system $catcommand;

print syslogfile "  Removing temporaty diamond files in $tempdir\n";
system("rm $tempdir/diamond_lca.*.m8");

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
        # print syslogfile "Total lines in Diamond: $wc; Allocating $splitlines in $numthreads threads\n";

        my $nextp=$splitlines;
        my ($filelines,$splitorf);
        my $numfile=1;
        # print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
        open(outfiletemp,">$tempdir/diamond_lca.$numfile.m8");
        open(infile2,$infile) || die "Can't open Diamond file $infile\n";
        while(<infile2>) {
                $filelines++;
                my @f=split(/\t/,$_);
                if($filelines==$nextp) { $splitorf=$f[0]; }
                if($filelines<=$nextp) { print outfiletemp $_; next; }
                elsif($f[0] ne $splitorf) {
                        close outfiletemp;
                        $numfile++;
                        # print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
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
	print syslogfile "Starting thread $threadnum\n";
	open(outc,">$tempdir/fun3tax\_$threadnum.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.wranks for writing\n";
	open(outcnof,">$tempdir/fun3tax\_$threadnum.nofilter.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.nofilter.wranks for writing\n";

	#-- Prepare the LCA database (containing the acc -> tax correspondence)

	my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1, sqlite_open_flags => SQLITE_OPEN_READONLY }) or die $DBI::errstr;
	$dbh->sqlite_busy_timeout( 120 * 1000 );

	my $currentfile="$tempdir/diamond_lca.$threadnum.m8";
	open(infile2,$currentfile) || die "Cannot open $currentfile\n";
	my $lastorf;
	while(<infile2>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));	
		$_=~s/\;\_//g;
		my @fields=split(/\t/,$_);
		my $thisorf=$fields[0];
		if($noidentical && (!$skipidentical) && ($fields[2] eq "100.0")) { $skipidentical=1; next; }			   
		if(!$refscore) { $refscore=$fields[11]; }
		if(!$refiden) { $refiden=$fields[2]; }  
		$posinit=$fields[6];			   
		$posend=$fields[7];
 		$provhits{$fields[1]}=$fields[11];
 		$providen{$fields[1]}=$fields[2];
		$tothits++;			   
		if($thisorf) { 
			# print "!!! $thisorf $lastorf\n";
			$lastorf=$thisorf;	
			query($dbh,$lastorf);
			(%provhits,%providen,%giden)=();
			($validhits,$validhitsnofilter,$tothits,$skipidentical)=0;
			$string="";
			($refscore,$refiden)=0;	
			}
		}
	close infile2;
	close outc;
	close outcnof;
	close outsyslog;
	print "  Tax assignment done! Result stored in file $fun3tax\_$threadnum.wranks\n" if $verbose;
}


sub query {
	my $dbh=shift;
	my $lastorf=shift;
	my($refcc,$genocc,$unicc)=0; 
	# if($lastorf=~/NODE_1_length_433318_cov_12.8415_1/) { $verbose=1; } else { $verbose=0; }
	print "refscore: $refscore refiden: $refiden\n" if $verbose; 
	my (%giden,%bhit)=();
	my($besthit,$ratoscore,$idendiff,$lasttax,$nuquery);
	my $query="select * from taxid where (";
	foreach my $lhits(keys %provhits) {
		print ">*>$lhits $provhits{$lhits}\n" if $verbose;
		my @gh=split(/\;/,$lhits);
		foreach my $fhit(@gh) {
			my @e=split(/\|/,$fhit);
			my $thishit=$e[0];
			my $thisscore=$e[1];
			my $thisiden=$e[2];
			my $ratioscore;
			$giden{$thishit}=$thisiden;
			# print "  ----- $thishit $thisscore $refscore\n";
			if($refscore) { $ratioscore=$thisscore/$refscore; }
			next if($ratioscore<=$scoreratio6);
                        $nuquery++;
                        last if($nuquery>=100);
			if($refiden) { $idendiff=$refiden-$thisiden; }
			next if($idendiff>$diffiden6);    
			if($refcc) { $query.=" or "; }
 			else { $refcc=1; }
			$query.="id=\"$thishit\"";
			if(!$besthit) { $besthit=$thishit; }
			}
		}
				     
	$query.=");";	
	print "*$query*\n" if $verbose;			     
	my (%accum,%accumnofilter)=();		     
	if($refcc) {
		my $sth = $dbh->prepare($query);  
		$sth->execute();
		while(my @list=$sth->fetchrow()) {
			print "$lastorf\t@list\n" if $verbose;
			for(my $pos=2; $pos<=8; $pos++) {
				my $rank=$ranks[$pos-1];
				my $tax=$list[$pos];
				print " $rank $tax $giden{$list[0]} $idenrank{$rank}\n" if $verbose;
				if($list[0] eq $besthit) { $bhit{$rank}=$tax; }
				if($giden{$list[0]}>=$idenrank{$rank}) { $accum{$rank}{$tax}++; }		#-- and add a count for that taxon in that rank
				$accumnofilter{$rank}{$tax}++; 		#-- Not considering identity filters for ranks
				}
			#if(($list[8]) && ($giden{$list[0]}>=$idenrank{'superkingdom'})) { $validhits++;  }			#-- Count the number of valid hits
			if(($list[0])) { 
				$validhitsnofilter++; 
				if($giden{$list[0]}>=$idenrank{'superkingdom'}) { $validhits++; }		#-- Count the number of valid hits 
				}			
			
			}
		}
	my($required,$minreqhits);
	if($validhits==1) { $minreqhits=1; } else { $minreqhits=$minhits6; }
	if($flex6<1) { $required=$validhits-($flex6*$validhits); } else { $required=$validhits-$flex6; }
	print "$lastorf Hits: $tothits; Valid: $validhits; Min: $minreqhits; Required: $required\n" if $verbose;
	$lasttax="";			
	foreach my $k(@ranks) {
		print "   $k\n" if $verbose;
                my $maxp=0;
                foreach my $t(sort { $accum{$k}{$a}<=>$accum{$k}{$b}; } keys %{ $accum{$k} }) {
                        if($t && ($accum{$k}{$t}==$maxp)) { $lasttax=""; next; }    #-- Equality of hits, don´t choose any
                        if($t) { $maxp=$accum{$k}{$t}; }
                        print "      $t $accum{$k}{$t}\n" if $verbose;
			if(($accum{$k}{$t}>=$required) && ($accum{$k}{$t}>=$minreqhits)) { 
				next if(($t ne $bhit{$k}) && ($bhitforced));
				print "$k -> $t\n" if $verbose;
				$lasttax=$t; 
				#  if($t) { $string="$t;$string"; }
				}
			}
		last if($lasttax);		
		}
	
	if($validhitsnofilter==1) { $minreqhits=1; } else { $minreqhits=$minhits6; }
	if($flex6<1) { $required=$validhitsnofilter-($flex6*$validhitsnofilter); } else { $required=$validhitsnofilter-$flex6; }
	my $lasttaxnofilter="";			
	foreach my $k(@ranks) {
		foreach my $t(keys %{ $accumnofilter{$k} }) {
			if(($accumnofilter{$k}{$t}>=$required) && ($accumnofilter{$k}{$t}>=$minreqhits)) { $lasttaxnofilter=$t; }
			}

		last if($lasttaxnofilter);		
		}

 
	my $abb=$parents{$lasttax}{wranks};
	
	#-- Changing nomenclature to abbreviations
	
	$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
	#print outc "$lastorf\t$parents{$lasttax}{wranks}\n";		
	print outc "$lastorf\t$abb\n";		
	my $abb=$parents{$lasttaxnofilter}{wranks};
	$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
	# print outcnof "$lastorf\t$parents{$lasttaxnofilter}{wranks}\n";		
	print outcnof "$lastorf\t$abb\n";		
	print "$lastorf\t$parents{$lasttax}{noranks}\n" if $verbose;
	(%accum,%accumnofilter)=();	
       }
