#!/usr/bin/perl

$|=1;

$infile=$ARGV[0];
$outfile=$ARGV[1];

$mindif=0.1;  # percentage of difference in scores for assigning bestaver
$maxhits=1;   #  Maximum number of hits for averaging

open(out,">$outfile") || die "Cannot open outfile $outfile\n";
open(in,$infile) || die "Cannot open m8 file $infile\n";
print out "# Created by $0, ",scalar localtime,", from $infile\n";
print out "#ORF\tBESTHIT\tBESTAVER\n";
while(<in>) {
 chomp;
 next if(!$_ || ($_=~/^\#/));
 @f=split(/\t/,$_);
 if($currorf && ($f[0] ne $currorf)) {
  foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }
  @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
  $score1=$accum{$list[0]}; 
  $score2=$accum{$list[1]};
  $dif=(($score1-$score2)/$score1); 
  if($dif>$mindif) { $cog{$currorf}{bestaver}=$list[0]; }
  print out "$currorf\t$cog{$currorf}{besthit}\t$cog{$currorf}{bestaver}\n";
  $currorf=$f[0];
  (%accum,%count)=();
                       }
 if($f[1]<$f[3]) { $minali=$f[1]; } else { $minali=$f[3]; }
 $olap=($f[5]*3)/$minali;
 next if($olap<$minolap);
 $currorf=$f[0];
 @c=split(/\|/,$f[2]);
 $khit=$c[1];
 if(!$cog{$currorf}{besthit}) {
  $cog{$currorf}{besthit}=$khit;
 # print out "$currorf\t$kegg{$currorf}\n";
                               }
 next if($count{$khit}>=$maxhits);
 $count{$khit}++;
 $accum{$khit}+=$f[7];
	     }
close in;
  foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }
  @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
  $score1=$accum{$list[0]}; 
  $score2=$accum{$list[1]};
  $dif=(($score1-$score2)/$score1); 
  if($dif>$mindif) { $cog{$currorf}{bestaver}=$list[0]; }
  print out "$currorf\t$cog{$currorf}{besthit}\t$cog{$currorf}{bestaver}\n";
  $currorf=$f[0];
close out;

print "   Output created in $outfile\n";

