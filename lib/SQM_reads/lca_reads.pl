#!/usr/bin/perl

# LCA tax assignment from a m8 file (c) Javier Tamames Dec 2017

$|=1;

use DBI;
use DBD::SQLite::Constants qw/:file_open/;
use Tie::IxHash;
use Cwd;
use threads;
use strict;
my $pwd=cwd();

use File::Basename;
our $auxdir = dirname(__FILE__);

do "$auxdir/../../scripts/SqueezeMeta_conf.pl";
our($databasepath,$numthreads);

my $lca_db="$databasepath/LCA_tax/taxid.db";

my @ranks=('species','genus','family','order','class','phylum','superkingdom');
my %idenrank=('species',85,'genus',60,'family',55,'order',50,'class',46,'phylum',42,'superkingdom',40);
my $scoreratio=0.8;   #-- Ratio first score/currsent score for the hit to be considered
my $diffiden=10;       #-- Maximim identity difference with the first
my $flex=0.2;           #-- Allows this PERCENTAGE (if less than one) or NUMBER (if grater than one) of hits from different taxa than LCA
my $minhits=1;        #-- Minimum number of hits for the taxa (if there is only one valid hit, this value sets to one automatically
my $noidentical=0;  #-- Drops the first 100% identical hit (for mock)
my $miniden=30;
my $verbose=0;
my $bhitforced=0;	#-- Forces that assignment cannot differ from best hit

my $outname;
my $infile=$ARGV[0];
my $resultpath = dirname($infile);
my $outname = basename($infile);
$outname=~s/\.m8//;
my $tempdir = $resultpath;
if($ARGV[1]) { $numthreads=$ARGV[1]; } else { $numthreads=12; }

my(%provhits,%providen,%giden,%accum,%accumnofilter,%parents);
my($validhits,$validhitsnofilter,$tothits,$skipidentical,$refscore,$refiden,$string,$posinit,$posend,$syslogfile);
tie %provhits,"Tie::IxHash";
tie %accum,"Tie::IxHash";
tie %accumnofilter,"Tie::IxHash";

if(!$infile) { die "Usage: lca.pl <infile> <outfile>\n"; }

open(in,"$databasepath/LCA_tax/parents.txt") || die;
while(<in>) {
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
close in;

#-- Split Diamond file

splitfiles();

#-- Launch threads

print "  Starting multithread LCA in $numthreads threads: ";
my $threadnum;
for($threadnum=1; $threadnum<=$numthreads; $threadnum++) {
        print "$threadnum ";
        my $thr=threads->create(\&current_thread,$threadnum);
}
print "\n";
$_->join() for threads->list();

my $wrankfile="$resultpath/$outname.wranks";
my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.wranks "; }
$catcommand.=" > $wrankfile";
print "  Creating $wrankfile file\n";
system $catcommand;

my $wrankfile="$resultpath/$outname\_nofilter.wranks";
my $catcommand="cat ";
for(my $h=1; $h<=$numthreads; $h++) { $catcommand.="$tempdir/fun3tax\_$h.nofilter.wranks "; }
$catcommand.=" > $wrankfile";
print "  Creating $wrankfile file\n";
system $catcommand;

system("rm $tempdir/diamond_lca.*.m8");
system("rm $tempdir/fun3tax*");
system("rm $tempdir/wc");


sub splitfiles {
        print "  Splitting Diamond file\n";
        system("wc -l $infile > $tempdir/wc");
        open(intemp,"$tempdir/wc");
        my $wc=<intemp>;
        close intemp;
        chomp $wc;
        $wc=~s/\s+.*//;    #-- Number of lines in the diamond result
        my $splitlines=int($wc/$numthreads);

        my $nextp=$splitlines;
        my ($filelines,$splitorf);
        my $numfile=1;
        print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
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
                        print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
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
       open(outc,">$tempdir/fun3tax\_$threadnum.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.wranks for writing\n";
       open(outcnof,">$tempdir/fun3tax\_$threadnum.nofilter.wranks") || die "Can't open $tempdir/fun3tax\_$threadnum.nofilter.wranks for writing\n";

        #-- Prepare the LCA database (containing the acc -> tax correspondence)

        my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1, sqlite_open_flags => SQLITE_OPEN_READONLY }) or die $DBI::errstr;
        $dbh->sqlite_busy_timeout( 120 * 1000 );
        my $currentfile="$tempdir/diamond_lca.$threadnum.m8";
        open(inc,$currentfile) || die "Cannot open $currentfile\n";
	my($lastorf,$thisorf);
	while(<inc>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));	
		$_=~s/\;\_//g;
		my @fields=split(/\t/,$_);
  		$thisorf=$fields[0];
  		if(!$lastorf) { $lastorf=$thisorf; }
  		if($lastorf && ($thisorf ne $lastorf)) { 
   			# print "!!! $thisorf $lastorf\n";
   			query($dbh,$lastorf);
   			(%accum,%accumnofilter,%provhits,%providen,%giden)=();
   			($validhits,$validhitsnofilter,$tothits,$skipidentical)=0;
   			$string="";
   			$lastorf=$thisorf;	
   			($refscore,$refiden)=0;	
			}
		if($noidentical && (!$skipidentical) && ($fields[2] eq "100.0")) { $skipidentical=1; next; }			   
		if(!$refscore) { $refscore=$fields[$#fields]; }
		if(!$refiden) { $refiden=$fields[2]; } 
		if($fields[$#fields]>$provhits{$fields[1]}) { $provhits{$fields[1]}=$fields[$#fields]; }
		if($fields[2]>$providen{$fields[1]}) { $providen{$fields[1]}=$fields[2]; }
		# print "PROVHIT: $fields[1] $provhits{$fields[1]} $providen{$fields[1]}\n";
		$tothits++;			   
		}
	close inc;
	$lastorf=$thisorf;
	query($dbh,$lastorf);    
	close outc;
	close outcnof;
        print "  Tax assignment done! Result stored in file fun3tax\_$threadnum.wranks\n" if $verbose;
}


sub query {
	my $dbh=shift;
	my $lastorf=shift;
   	my($refcc,$genocc,$unicc)=0; 
	  # if($lastorf=~/NODE_1_length_433318_cov_12.8415_1/) { $verbose=1; } else { $verbose=0; }
	print "refscore: $refscore refiden: $refiden\n" if $verbose; 
	my(%giden,%bhit)=();
	my($besthit,$ratioscore,$idendiff,$lasttax,$nuquery,$minreqhits,$required,$lasttaxnofilter,$thereareresults);
	my $query="select * from taxid where (";
	foreach my $lhits(keys %provhits) {
    		if($refscore) { $ratioscore=$provhits{$lhits}/$refscore; }
    		print ">*>$lhits $provhits{$lhits} $ratioscore $providen{$lhits}\n" if $verbose;
    		next if($ratioscore<=$scoreratio);
    		next if($providen{$lhits}<$miniden);
    		if($refiden) { $idendiff=$refiden-$providen{$lhits}; }
    		next if($idendiff>$diffiden);    
    		my @hitfields=split(/\|/,$lhits);
     		$giden{$hitfields[0]}=$providen{$lhits}; 
     		if($refcc) { $query.=" or "; }
      		else { $refcc=1; } 
      		$query.="id=\"$hitfields[0]\"";
      		if(!$besthit) { $besthit=$hitfields[0]; }
		}
    	$query.=");";	
    	print "*$query*\n" if $verbose;			     
    	if($refcc) {
		my $sth = $dbh->prepare($query);  
		$sth->execute();
		while(my @list=$sth->fetchrow()) {
		print "$lastorf\t@list\n" if $verbose;
		for(my $pos=1; $pos<=7; $pos++) {
			my $rank=$ranks[$pos-1];
			my $tax=$list[$pos];
			if($list[0] eq $besthit) { $bhit{$rank}=$tax; }
			if($giden{$list[0]}>=$idenrank{$rank}) { 
				$accum{$rank}{$tax}++; #-- and add a count for that taxon in that rank
				print ">>>> $pos $rank $tax $accum{$rank}{$tax}\n" if $verbose;
				}		
 			$accumnofilter{$rank}{$tax}++; 
     			#  $accum{$rank}{$tax}++;
			}
		if(($list[7]) && ($giden{$list[0]}>=$idenrank{'superkingdom'})) { $validhits++;  }			#-- Count the number of valid hits
		if(($list[7])) { $validhitsnofilter++;  }			#-- Count the number of valid hits
		}
	}
	if($validhits==1) { $minreqhits=1; } else { $minreqhits=$minhits; }
   	if($flex<1) { $required=$validhits-($flex*$validhits); } else { $required=$validhits-$flex; }
   	print "$lastorf Hits: $tothits; Valid: $validhits; Min: $minreqhits; Required: $required\n" if $verbose;
   	$lasttax="";			
   	foreach my $k(@ranks) {
    		print "   $k\n" if $verbose;
		my $maxp=0;
                foreach my $t(sort { $accum{$k}{$a}<=>$accum{$k}{$b}; } keys %{ $accum{$k} }) {
                        if($t && ($accum{$k}{$t}==$maxp)) { $lasttax=""; next; }    #-- Equality of hits, don´t choose any
                        if($t) { $maxp=$accum{$k}{$t}; }
     			print "      $t $accum{$k}{$t}\n" if $verbose;
     			if(($accum{$k}{$t}>=$required) && ($accum{$k}{$t}>=$minreqhits) && ($parents{$t}{wranks})) { 
      			next if(($t ne $bhit{$k}) && ($bhitforced));
      			print "$k -> $t\n" if $verbose;
      			$lasttax=$t; 
    			#  if($t) { $string="$t;$string"; }
			}
		}
		last if($lasttax);		
	}
    	if($validhitsnofilter==1) { $minreqhits=1; } else { $minreqhits=$minhits; }
	if($flex<1) { $required=$validhitsnofilter-($flex*$validhitsnofilter); } else { $required=$validhitsnofilter-$flex; }
	$lasttaxnofilter="";			
	foreach my $k(@ranks) {
		foreach my $t(keys %{ $accumnofilter{$k} }) {
			if(($accumnofilter{$k}{$t}>=$required) && ($accumnofilter{$k}{$t}>=$minreqhits) && ($parents{$t}{wranks})) { 
				$lasttaxnofilter=$t; 
				}
			}
		last if($lasttaxnofilter);		
		}
	my $abb=$parents{$lasttax}{wranks};
	
	#-- Changing nomenclature to abbreviations
	
	$abb=~s/sub\w+\:/n_/g;
	$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g;
	# print outfile2 "$lastorf\t$parents{$lasttax}{wranks}\n";		
	print outc "$lastorf\t$abb\n";		
		# print outfile3 "$lastorf\t$parents{$lasttaxnofilter}{noranks}\n";
		my $abb=$parents{$lasttaxnofilter}{wranks};
		$abb=~s/superkingdom\:/k_/; $abb=~s/phylum\:/p_/; $abb=~s/order\:/o_/; $abb=~s/class\:/c_/; $abb=~s/family\:/f_/; $abb=~s/genus\:/g_/; $abb=~s/species\:/s_/; $abb=~s/no rank\:/n_/g; $abb=~s/\w+\:/n_/g; 
		# print outfile4 "$lastorf\t$parents{$lasttaxnofilter}{wranks}\n";	
		print outcnof "$lastorf\t$abb\n";	
	 print "$lastorf\t$parents{$lasttax}{noranks}\n" if $verbose;
	#print "$lastorf\t$abb\n" if $verbose;
	if($parents{$lasttax}{noranks}) { $thereareresults=1; }	
       }
