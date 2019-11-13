#!/usr/bin/perl

# LCA tax assignment from a m8 file (c) Javier Tamames Dec 2017

$|=1;

use DBI;
use Tie::IxHash;
use Cwd;
$pwd=cwd();

use File::Basename;
our $auxdir = dirname(__FILE__);

do "$auxdir/../../scripts/SqueezeMeta_conf.pl";
our($databasepath);

my $lca_db="$databasepath/LCA_tax/taxid.db";

my $dbh = DBI->connect("dbi:SQLite:dbname=$lca_db","","",{ RaiseError => 1}) or die $DBI::errstr;

@ranks=('species','genus','family','order','class','phylum','superkingdom');
my %idenrank=('species',85,'genus',60,'family',55,'order',50,'class',46,'phylum',42,'superkingdom',40);
$scoreratio=0.8;   #-- Ratio first score/currsent score for the hit to be considered
$diffiden=10;       #-- Maximim identity difference with the first
$flex=0.2;           #-- Allows this PERCENTAGE (if less than one) or NUMBER (if grater than one) of hits from different taxa than LCA
$minhits=1;        #-- Minimum number of hits for the taxa (if there is only one valid hit, this value sets to one automatically
$noidentical=0;  #-- Drops the first 100% identical hit (for mock)
$miniden=30;
$verbose=1;
$bhitforced=0;	#-- Forces that assignment cannot differ from best hit

$infile=$ARGV[0];
$outname=$infile;
$outname=~s/\.m8//;

if(!$outname) { die "Usage: lca.pl <infile> <outfile>\n"; }

open(in,"$databasepath/LCA_tax/parents.txt") || die;
while(<in>) {
 chomp;
 next if !$_;
 ($tax,$par)=split(/\t/,$_);
 $parents{$tax}{wranks}=$par;
 @m=split(/\;/,$par);
 foreach my $y(@m) {
  ($rt,$gtax)=split(/\:/,$y);
  $parents{$tax}{noranks}.="$gtax;"; 
                   }
 chop $parents{$tax}{noranks};
             }
close in;

 tie %provhits,"Tie::IxHash";
 tie %accum,"Tie::IxHash";
 if($infile=~/\.gz$/) { open(inc,"zcat $infile|") || die; }
 else { open(inc,$infile) || die "Cannot open Diamond file $infile\n"; }
 my $outfilter="$outname.wranks";
 my $outnofilter="$outname\_nofilter.wranks";
 
# open(out,">$outname") || die;
 open(outc,">$outfilter") || die;
# open(outnof,">$outname\_nofilter") || die;
 open(outcnof,">$outnofilter") || die;
 while(<inc>) { 
  chomp;
  next if(!$_ || ($_=~/^\#/));	
  $_=~s/\;\_//g;
  @fields=split(/\t/,$_);
  $thisorf=$fields[0];
  if(!$lastorf) { $lastorf=$thisorf; }
  if($lastorf && ($thisorf ne $lastorf)) { 
   # print "!!! $thisorf $lastorf\n";
   query();
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
query();    
# close out;
close outc;
#close outnof;
close outcnof;
# print "Tax assignment done! Result stored in file $outfilter\n";


sub query {
   ($refcc,$genocc,$unicc)=0; 
  # if($lastorf=~/NODE_1_length_433318_cov_12.8415_1/) { $verbose=1; } else { $verbose=0; }
  print "refscore: $refscore refiden: $refiden\n" if $verbose; 
   (%giden,%bhit)=();
   my $besthit;
  $query="select * from taxid where (";
   foreach my $lhits(keys %provhits) {
    if($refscore) { $ratioscore=$provhits{$lhits}/$refscore; }
    print ">*>$lhits $provhits{$lhits} $ratioscore $providen{$lhits}\n" if $verbose;
    next if($ratioscore<=$scoreratio);
    next if($providen{$lhits}<$miniden);
    if($refiden) { $idendiff=$refiden-$providen{$lhits}; }
    next if($idendiff>$diffiden);    
    @hitfields=split(/\|/,$lhits);
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
     while(@list=$sth->fetchrow()) {
     print "$lastorf\t@list\n" if $verbose;
      for(my $pos=1; $pos<=7; $pos++) {
       $rank=$ranks[$pos-1];
       $tax=$list[$pos];
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
    foreach my $t(keys %{ $accum{$k} }) {
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

 # print out "$lastorf\t$parents{$lasttax}{noranks}\n";
 # print outc "$lastorf\t$parents{$lasttax}{wranks}\n";		
 # print outnof "$lastorf\t$parents{$lasttaxnofilter}{noranks}\n";
 # print outcnof "$lastorf\t$parents{$lasttaxnofilter}{wranks}\n";		
 # print "$lastorf\t$parents{$lasttax}{noranks}\n" if $verbose;	
       }
