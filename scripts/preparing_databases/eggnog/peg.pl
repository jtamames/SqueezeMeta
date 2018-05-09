#!/usr/bin/perl

# Genera eggnog4.gb

open(in,"/home/bcamara/database/NOG.members.tsv") || die;  # De http://eggnogdb.embl.de/download/eggnog_4.5/data/NOG/
while(<in>) {
 chomp;
 @f=split(/\t/,$_);
 @e=split(/\,/,$f[5]);
 map { $cogs{$_}=$f[1]; } @e;
            }
close in;

open(in,"eggnog4.clustered_proteins.blastp");  # Resultado de correr Diamond de eggnog4.proteins.all.fa (http://eggnogdb.embl.de/download/eggnog_4.5) sobre nr
while(<in>) {
 chomp;
 @f=split(/\t/,$_);
 if(($f[2]=~/100/) && ($f[4] eq "0")  && ($f[5] eq "0")) { 
  print "$f[0]\t$f[1]\t$cogs{$f[0]}\n"; 
                                                         }
            }
close in;
