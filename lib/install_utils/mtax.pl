#!/usr/bin/env perl

#-- Creates alltaxlist.txt, needed for mcount

open(in,"/media/mcm/jtamames/databases/taxonomy/names.dmp") || die;
while(<in>) {
 chomp;
 next if !$_;
 $_=~s/\t\|$//g;
 @t=split(/\t\|\t/,$_);
 next if($t[3]!~/scientific/);
 # print "$t[0]\t$t[1]\n";
 $names{$t[0]}=$t[2];
            }
close in;

open(in,"/media/mcm/jtamames/databases/taxonomy/nodes.dmp") || die;
while(<in>) {
 chomp;
 next if !$_;
 $_=~s/\t\|$//g;
 @t=split(/\t\|\t/,$_);
 if($names{$t[0]}) { print "$t[0]\t$names{$t[0]}\t$t[2]\n"; }
             }
close in;
