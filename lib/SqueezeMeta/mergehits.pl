#!/usr/bin/env perl

use strict;

# -- This program joins ORFs that are contiguous, partial, and match the same entry, for instance these product of frameshifts

$|=1;

my $file=$ARGV[0];
my(%store,%frames);
my($nh,$current,$last);
my $verbose=0;

open(in,$file) || die "Can't open $file\n";
while(<in>) {
	chomp;
	my @fields=split(/\t/,$_);
	#my @culist=split(/\|/,$fields[0]);
	#my $contig=$culist[0];
	#$contig=~s/\_[^_]+$//;
	my $contig=$fields[0];
	$contig=~s/\_[^_]+$//;
	my $current=$contig;
print "$current\n";
	if($current eq $last) {
		my @f=split(/\t/,$_);
		my $long=$f[7]-$f[6];
		$store{$_}=$long;
		$frames{"$current $f[6]-$f[7]"}=$f[$#f];
		$nh++;
		}
	else {
		if($nh>1) { mergehits($last); }
		foreach my $pk(keys %store) { print "$pk\n"; }
		%store=();
		%frames=();
		$nh=1;
		my @f=split(/\t/,$_);
		my $long=$f[7]-$f[6];
		$store{$_}=$long;
		$frames{"$current $f[6]-$f[7]"}=$f[$#f];
		}
	$last=$current;
	}
close in;
if($nh>1) { mergehits($last); }


sub mergehits {
	my $currcontig=shift;
	my $changes=1;
	my $merged=0;
	while($changes) {
		$changes=0;
		my @ak=sort { $store{$b}<=>$store{$a}; } keys %store;
		#  return if(!$#ak);
		for(my $pos1=0; $pos1<=$#ak; $pos1++) {
			for(my $pos2=$pos1+1; $pos2<=$#ak; $pos2++) {
			my %accum=();	  
			my($in1,$in2,$tomerge)=0;  
			my @fields1=split(/\t/,$ak[$pos1]);
			my @fields2=split(/\t/,$ak[$pos2]);
			my @ids1=split(/\;/,$fields1[1]);    
			my @ids2=split(/\;/,$fields2[1]);
			foreach my $nh(@ids1) {
				my @sf=split(/\|/,$nh);
				if($sf[0] eq "gi") { $tomerge=$sf[1]; } else { $tomerge=$sf[0]; }
				if($accum{$tomerge}!~/1/) { $in1++; }
				$accum{$tomerge}="1";
				}         		       
			foreach my $nh(@ids2) {
				my @sf=split(/\|/,$nh);
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
			my $mperc=$multiple/$mli;
			my @newl;
			if($mperc>=0.5) { 
				#  if($multiple) { 
				#   print "Unir:\n**$ak[$pos1]\n**$ak[$pos2]\n$multiple\n\n\n"; 
				my $long1=$fields1[3];
				my $long2=$fields2[3];
				if($verbose) { print " >>Joining because $muls:\n**$ak[$pos1]\n**$ak[$pos2]\n\n\n"; }
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
				if($newl[6]<$newl[7]) { splice(@newl,0,1,"$currcontig\_$newl[6]-$newl[7]"); }
				else { splice(@newl,0,1,"$currcontig\_$newl[7]-$newl[6]"); }
				my $newlong=$newl[9]-$newl[8];
				splice(@newl,3,1,$newlong);
				my $newline=join("\t",@newl);
				if($verbose) { print " --Creating:\n**$newline\n\n\n"; }
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
}

sub min {
	my $a1=shift;
	my $a2=shift;
	if($a1<$a2) { return $a1; } else { return $a2; } 
        }
 
sub max {
	my $a1=shift;
	my $a2=shift;
	if($a1<$a2) { return $a2; } else { return $a1; } 
        }
 
 
