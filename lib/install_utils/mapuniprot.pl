#!/usr/bin/env perl

$|=1;

$command="wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz";
system $command;

$file="/home/tamames/fun/fun3/uniprot_trembl.dat.gz";

print "# Created by $0 with data from $file, ",scalar localtime,"\n"; 

open(in,"zcat $file|") || die "Cannot open $file\n";
while(<in>) {
 chomp;
 next if !$_;
 if($_=~/^DR\s+EMBL/) {
  @f=split(/\; /,$_);
  $acc=$f[2];
  next if($acc!~/\w/);
  push(@gid,$acc);

                     } 
 if($_=~/^DR\s+RefSeq/) {
  @f=split(/\; /,$_);
  $acc=$f[1];
  next if($acc!~/\w/);
  push(@rid,$acc);

                     } 
 if($_=~/^DR\s+eggNOG/) {
  @f=split(/\; /,$_);
  $acc=$f[1];  
  next if($acc!~/\w/);
  push(@cid,$acc);
  $m=1;
                     }
 if($_=~/^DR\s+KO/) {
  @f=split(/\; /,$_);
  $acc=$f[1];
  next if($acc!~/\w/);
  push(@kid,$acc);
  $m=1;
                     }
if($_=~/\/\//) {
 if($m) {
  $cog=join(";",sort @cid);
  $kegg=join(";",sort @kid);
  foreach my $e(@rid) {
 #  print "$e\t$cog\t$kegg\n";
   print "$e\t$kegg\n";
                      }
	}
 $m=0;
 (@gid,@rid,@cid,@kid)=();
              }
         }   		     
close in;
