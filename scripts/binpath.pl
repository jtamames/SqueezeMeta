#!/usr/bin/perl

$|=1;

$mincompletion=50;
$pathreq=$ARGV[0];

open(in,"/media/mcm/jtamames/metagenomic_data/kegg2ec.txt") || die;
while(<in>) {
 chomp;
 next if !$_;
 @k=split(/\t/,$_);
 $k[1]=~s/EC\://g;
 $ec{$k[0]}{$k[1]}=1;
            }
close in;

open(in,"/media/mcm/jtamames/metagenomic_data/ontopath.txt") || die;
while(<in>) {
 chomp;
 next if !$_;
 ($id,$rest)=split(/ -> /,$_);
 @onto=split(/\t/,$rest);
 $pathname=$onto[0];
 $biocyc{$id}=$pathname;
 for(my $pos=1; $pos<=$#onto; $pos++) { $ontocyc{$pathname}{$pos}{$onto[$pos]}=1; }
 $fullpath=join(";",@onto);
 $upper{$pathname}=$fullpath;
            }
close in;

open(in,"/media/mcm/jtamames/artico022/results/09.artico022.bintable") || die;
while(<in>) {
 chomp;
 next if(!$_ || ($_=~/^\#/));
 @e=split(/\t/,$_);
 next if($e[1] ne "metabat2"); 
 $bintax{$e[0]}=$e[2];
 $bincomplete{$e[0]}=$e[8];
            }
close in;

open(in,"/media/mcm/jtamames/artico022/results/09.artico022.contigsinbins") || die;
while(<in>) {
 chomp;
 next if(!$_ || ($_=~/^\#/));
 @e=split(/\t/,$_);
 next if($e[1] ne "metabat2"); 
 $contig{$e[2]}{$e[0]}=$e[1];
            }
close in;

open(in,"/media/mcm/jtamames/artico022/results/04.artico022.fun3.kegg");
while(<in>) {
 chomp;
 next if(!$_ || ($_=~/^\#/));
 @u=split(/\t/,$_);
 $contig=$u[0];
 $contig=~s/\_\d+$//;
 $keggid=$u[2];
 next if !$keggid;
 foreach my $thisec(keys %{ $ec{$keggid} }) { $accum{$contig}{$thisec}=$u[0]; }
            }
close in;

foreach my $thisbin(sort keys %contig) {
 $minfile="/media/mcm/jtamames/artico022/temp/$thisbin.ec";
 $report="/media/mcm/jtamames/artico022/temp/$thisbin.report";
 open(outb,">$minfile") || die "Cannot open $minfile\n";
 # print "Bin $thisbin\n";
 foreach my $thiscontig(sort keys %{ $contig{$thisbin} }) { 
  foreach my $thisec(sort keys %{ $accum{$thiscontig} }) {
   print outb "$accum{$thiscontig}{$thisec} $thisec\n";
                                                         }
                                                          }
 close outb;
 $command="/home/jtamames/software/MinPath/MinPath1.4.py -any $minfile -map ec2path -report $report > /dev/null";
# system $command; 
 open(inr,$report) || warn "Cannot open $report\n";
 while(<inr>) {
  chomp;
  next if !$_;
  next if($_!~/minpath 1/);
  @k=split(/\s+/,$_);
  $pathid=$k[$#k];
  $pathname=$biocyc{$pathid};
  $paths{$pathname}{$thisbin}=1;
              }
                                        }

foreach my $apath(sort keys %paths) { 
 $ont=$upper{$apath};
 if((!$pathreq) || ($ont=~/\Q$pathreq\E/i)) { $goodpath{$apath}=1; }
                                    }					

$outfile="/media/mcm/jtamames/artico022/temp/pathway.table";
open(outt,">$outfile") || die;
print outt "\tTax\tCompleteness";
foreach my $apath(sort keys %goodpath) { print outt "\t$apath"; }
print outt "\n";
#foreach my $bin(sort { $bintax{$a} cmp $bintax{$b} } keys %bintax) {
foreach my $bin(sort { $bincomplete{$b}<=>$bincomplete{$a} } keys %bincomplete) {
 next if($bincomplete{$bin}<$mincompletion);
 print outt "$bin\t$bintax{$bin}\t$bincomplete{$bin}";
 foreach my $apath(sort keys %goodpath) { 
  $value=$paths{$apath}{$bin} || "0";
  print outt "\t$value";	
                                     }
 print outt "\n";
                                   }			   
close outt;   		
print "Output in $outfile\n"; 
 
 
