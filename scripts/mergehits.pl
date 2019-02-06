#!/usr/bin/perl


# -- Une hits contiguos y no necesariamente solapantes, como los procedentes de frameshifts

$|=1;

$file=$ARGV[0];
open(in,$file) || die "Cannot open $file\n";
while(<in>) {
 chomp;
 @fields=split(/\t/,$_);
 @culist=split(/\|/,$fields[0]);
 $contig=$culist[0];
 $contig=~s/\_[^_]+$//;
 $current=$contig;
 if($current eq $last) {
  @f=split(/\t/,$_);
  $long=$f[7]-$f[6];
  $store{$_}=$long;
  $frames{"$current $f[6]-$f[7]"}=$f[$#f];
  $nh++;
                        }
 else {
  if($nh>1) { mergehits(); }
  foreach my $pk(keys %store) { print "$pk\n"; }
  %store=();
  %frames=();
  $nh=1;
  @f=split(/\t/,$_);
  $long=$f[7]-$f[6];
  $store{$_}=$long;
  $frames{"$current $f[6]-$f[7]"}=$f[$#f];
      }
 $last=$current;
            }
close in;
if($nh>1) { mergehits(); }


sub mergehits {
 $changes=1;
 $merged=0;
 while($changes) {
  $changes=0;
  @ak=sort { $store{$b}<=>$store{$a}; } keys %store;
 #  return if(!$#ak);
  for(my $pos1=0; $pos1<=$#ak; $pos1++) {
   for(my $pos2=$pos1+1; $pos2<=$#ak; $pos2++) {
    %accum=();	  
    ($in1,$in2)=0;  
    @fields1=split(/\t/,$ak[$pos1]);
    @fields2=split(/\t/,$ak[$pos2]);
    @ids1=split(/\;/,$fields1[1]);    
    @ids2=split(/\;/,$fields2[1]);
    foreach my $nh(@ids1) {
     @sf=split(/\|/,$nh);
     if($sf[0] eq "gi") { $tomerge=$sf[1]; } else { $tomerge=$sf[0]; }
     if($accum{$tomerge}!~/1/) { $in1++; }
     $accum{$tomerge}="1";
                          }         		       
    foreach my $nh(@ids2) {
     @sf=split(/\|/,$nh);
     if($sf[0] eq "gi") { $tomerge=$sf[1]; } else { $tomerge=$sf[0]; }
     if($accum{$tomerge}!~/2/) { $in2++; }
     $accum{$tomerge}.="2";
                          }
    my $mli=min($in1,$in2);		    
    my($single,$multiple)=0;	
    my $muls;
    foreach my $print(keys %accum) { 
     if($accum{$print}=~/12/) { $multiple++; $muls.="$print;"; } else { $single++; }	       		       
                                   }
    $mperc=$multiple/$mli;
    if($mperc>=0.5) { 
  #  if($multiple) { 
  #   print "Unir:\n**$ak[$pos1]\n**$ak[$pos2]\n$multiple\n\n\n"; 
     $long1=$fields1[3];
     $long2=$fields2[3];
  #   print "Unir por $muls:\n**$ak[$pos1]\n**$ak[$pos2]\n\n\n"; 
     if($long1>$long2) { @newl=@fields1; } else { @newl=@fields2; }
  #   $frames{"$fields1[6]-$fields1[7]"}=$fields1[$#fields1];
  #   $frames{"$fields2[6]-$fields2[7]"}=$fields2[$#fields2];
     if(($fields1[6]<=$fields1[7]) && ($fields2[6]<=$fields2[7])) {
      splice(@newl,6,1,min($fields1[6],$fields2[6]));
      splice(@newl,7,1,max($fields1[7],$fields2[7]));
                                 }
     elsif(($fields1[6]>$fields1[7]) && ($fields2[6]<=$fields2[7])) {
      splice(@newl,6,1,min($fields1[7],$fields2[6]));
      splice(@newl,7,1,max($fields1[6],$fields2[7]));
                                 }
     elsif(($fields1[6]<=$fields1[7]) && ($fields2[6]>$fields2[7])) {
      splice(@newl,6,1,min($fields1[6],$fields2[7]));
      splice(@newl,7,1,max($fields1[7],$fields2[6]));
                                 }
     elsif(($fields1[6]>$fields1[7]) && ($fields2[6]>$fields2[7])) {
      splice(@newl,6,1,min($fields1[7],$fields2[7]));
      splice(@newl,7,1,max($fields1[6],$fields2[6]));
                                 }
      			 
     splice(@newl,8,1,min($fields1[8],$fields2[8]));
     splice(@newl,9,1,max($fields1[9],$fields2[9]));
     if($newl[6]<$newl[7]) { splice(@newl,0,1,"$current\_$newl[6]-$newl[7]"); }
     else { splice(@newl,0,1,"$current\_$newl[7]-$newl[6]"); }
     $newlong=$newl[9]-$newl[8];
     splice(@newl,3,1,$newlong);
     $newline=join("\t",@newl);
  #   print "Genera:\n**$newline\n\n\n";
  # print "$newline\n";     
     delete $store{$ak[$pos1]};
     delete $store{$ak[$pos2]};     
     $store{$newline}=1;
     $changes=1;
     $merged=1;
                  }
      last if $changes;	  
                                          }
      last if $changes;	  
					 }
                   }
 if($merged) { 
 # foreach my $tf(sort keys %frames) { print "Frame: $tf\t$frames{$tf}\n"; }
             }
                }

sub min {
 ($a1,$a2)=@_;
 if($a1<$a2) { return $a1; } else { return $a2; } 
        }
 
sub max {
 ($a1,$a2)=@_;
 if($a1<$a2) { return $a2; } else { return $a1; } 
        }
 
 
