#!/usr/bin/perl

$|=1;

#-- Creates a list of species and their corresponding taxa, from taxonomy (ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy)
#-- Javier Tamames, Dec 2017

$outfile="taxatree.txt";
$parentfile="parents.txt";

open(out,">$outfile");
open(in,"/media/mcm/jtamames/databases/taxonomy/nodes.dmp") || die;
while(<in>) {
 chomp;
 next if !$_;
 @t=split(/\t\|\t/,$_);
 $tid=$t[0]; 
 # next if($t[4]=~/3|7|8|9|11/);
 next if($t[4]=~/7|8|11/);
 $ranks{$tid}=$t[2]; 
 $parents{$tid}=$t[1]; 
 if($t[4]=~/3|9/) { $virus{$tid}=1; }
            }
close in;

open(in,"/media/mcm/jtamames/databases/taxonomy/names.dmp") || die;
while(<in>) {
 chomp;
 next if !$_;
 @t=split(/\t\|\t/,$_);
 $corr=$store{$t[0]};
 next if($t[3]!~/scientific/);
 $names{$t[0]}=$t[2];
	    }
close in;

foreach my $print(sort keys %ranks) {
 # print "*****$print\t$names{$print}\t$ranks{$print}\n";
 next if($ranks{$print} ne "species");
 @k=split(/\s+/,$names{$print});
 if($virus{$print}) { $specname=$names{$print}; }
 elsif($k[0] eq "Candidatus") { $specname="$k[0] $k[1] $k[2]"; } 
 else { $specname="$k[0] $k[1]"; }
 next if($seen{$specname});
 $seen{$specname}=1;
 $string="$ranks{$print}:$print:$specname";
  $found=1;
 $current=$print;
 while($found) {
  $parent=$parents{$current};
  # print "   >>$parent $ranks{$parent}\n";
  if($parent) {
  # if($ranks{$parent} ne "no rank") { 
  #  $string.=";$ranks{$parent}:$parent:$names{$parent}";
   #                                  }
   $string.=";$ranks{$parent}:$parent:$names{$parent}";
   $current=$parent;
              }
  else { $found=0; }
  if($ranks{$parent} eq "superkingdom") { $found=0; }
  if($parent==1) { $found=0; }
               }
 print out "$string\n";
	                            }
close out;
print "Output created in $outfile\n";
   

open(out,">$parentfile") || die;
open(in,"$outfile") || die;
while(<in>) {
 chomp;
 next if !$_;
 @k=split(/\;/,$_);
 for(my $pos1=0; $pos1<=$#k; $pos1++) {
  ($root,$string)="";
  for(my $pos2=$pos1; $pos2<=$#k; $pos2++) { 
   @n=split(/\:/,$k[$pos2]);
  # $string.="$n[0]:$n[2];"; 
   $string="$n[0]:$n[2];".$string; 
   if(!$root) { $root="$n[2]"; }
                                           }
  chop $string;	
  $string=~s/ \<prokaryotes\>//;					   
  if(!$yseen{$root}) { print out "$root\t$string\n"; }
  $yseen{$root}=1;					   
				      }
	   }
close in;
close out;
print "Output created in $parentfile\n";

