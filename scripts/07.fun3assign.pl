#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Performs functional assignment (COG, KEGG, Pfam) using previous Diamond/hmmer runs

use strict;
use Cwd;

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 

do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$nocog,$nokegg,$nopfam,$cogdiamond,$fun3cog,$evalue,$miniden,$keggdiamond,$fun3kegg,$pfamlist,$fun3pfam,$pfamhmmer);

#-- Some local conf variables

my $mindif=0.1;		#-- Percentage of difference in scores for assigning bestaver
my $maxhits=3;		#-- Maximum number of hits for averaging
my $minolap=30;		#-- Minimum alignment percentage between query and hit

print "Assigning";
if(!$nocog) { print " COGS"; }
if(!$nokegg) { print " KEGG"; }
if(!$nopfam) { print " PFAM"; }
print "\n";

#----------------------------------- COG assignment -------------------------------------

if(!$nocog) {
	open(infile1,$cogdiamond) || die "Cannot open cog file $cogdiamond\n";
	open(outfile1,">$fun3cog") || die "Cannot open $fun3cog\n";
	print outfile1 "# Created by $0, ",scalar localtime,", evalue=$evalue, miniden=$miniden, minolap=$minolap\n";
	print outfile1 "#ORF\tBESTHIT\tBESTAVER\n";

	#-- We start reading the Diamond against COG datafile

	my($currorf,$minali,$score1,$score2,$dif);
	my(%accum,%count,%cog);
	my @f;
	while(<infile1>) { 
 		chomp;
 		next if(!$_ || ($_=~/^\#/));
 		@f=split(/\t/,$_);
		
		#-- If we finished reading the hits for last ORF, we output its best hit,
		#-- and also calculate and output the best average, if any
		
 		if($currorf && ($f[0] ne $currorf)) {	
		
			#-- Calculate the best average
		
			foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }	#-- Average score for the COG
  			my @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
  			$score1=$accum{$list[0]};		#-- Score of the best COG
			if($score1) {
  				if($list[1]) { $score2=$accum{$list[1]} } else { $score2=0; } 	#-- Score of the 2nd best
  				$dif=(($score1-$score2)/$score1); 
  				if($dif>$mindif) { $cog{$currorf}{bestaver}=$list[0]; }	#-- If difference is enough, then the first is best average

				#-- Output the result (Best hit and best aver)

  				print outfile1 "$currorf\t$cog{$currorf}{besthit}\t$cog{$currorf}{bestaver}\n";
				    }
  			$currorf=$f[0];		#-- And current ORF is the new one just read
  			(%accum,%count)=();
			}
			
		#-- If we are still reading the hits for current ORF, just store them if they pass the filters

		$currorf=$f[0];
		if($f[1]<$f[3]) { $minali=$f[1]; } else { $minali=$f[3]; }
		my $olap=$f[5]*100/$minali;		#-- Percentage of the query covered by the hit
		next if($olap<$minolap);	#-- Partial hits are not allowed
		my @c=split(/\|/,$f[2]);
		my $khit=$c[1];			#-- This is the COG for the hit
		if(!$cog{$currorf}{besthit}) {
			$cog{$currorf}{besthit}=$khit;	#-- If it is the first hit, then it is best hit
			# print out "$currorf\t$kegg{$currorf}\n";
			}
		next if($count{$khit} && ($count{$khit}>=$maxhits));	#-- If we have already $maxhits hits for that COG, skip this hit
		$count{$khit}++;
		$accum{$khit}+=$f[$#f];
		}
		
	close infile1;

	#-- We need to proccess also the last ORF in the file

	foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }
	my @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
	$score1=$accum{$list[0]}; 
	if($score1) {
 		if($list[1]) { $score2=$accum{$list[1]} } else { $score2=0; } 	#-- Score of the 2nd best
		$dif=(($score1-$score2)/$score1); 
		if($dif>$mindif) { $cog{$currorf}{bestaver}=$list[0]; }
		print outfile1 "$currorf\t$cog{$currorf}{besthit}\t$cog{$currorf}{bestaver}\n";
		}
	close outfile1;
	
	}		#-- END of COG assignment


#----------------------------------- KEGG assignment -------------------------------------

if(!$nokegg) {
	open(infile2,$keggdiamond) || die "Cannot open $keggdiamond\n";
	open(outfile2,">$fun3kegg") || die "Cannot open $fun3kegg\n";
	print outfile2 "# Created by $0, ",scalar localtime,", evalue=$evalue, miniden=$miniden, minolap=$minolap\n";
	print outfile2 "#ORF\tBESTHIT\tBESTAVER\n";

	#-- We start reading the Diamond against KEGG datafile

	my($currorf,$minali,$score1,$score2,$dif);
	my(%accum,%count,%kegg);
	my @f;
	
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		@f=split(/\t/,$_);
		
		#-- If we finished reading the hits for last ORF, we output its best hit,
		#-- and also calculate and output the best average, if any

		if($currorf && ($f[0] ne $currorf)) {
		
			#-- Calculate the best average
		
			foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }	#-- Average score for the KEGG
			my @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
			$score1=$accum{$list[0]};		#-- Score of the best KEGG
			if($score1) { 
 				if($list[1]) { $score2=$accum{$list[1]} } else { $score2=0; } 	#-- Score of the 2nd best
				$dif=(($score1-$score2)/$score1); 
				if($dif>$mindif) { $kegg{$currorf}{bestaver}=$list[0]; }

				#-- Output the result (Best hit and best aver)

 				print outfile2 "$currorf\t$kegg{$currorf}{besthit}\t$kegg{$currorf}{bestaver}\n";
				   }
			$currorf=$f[0];			#-- And current ORF is the new one just read
			(%accum,%count)=();
			}
			
		#-- If we are still reading the hits for current ORF, just store them if they pass the filters
			
			
		if($f[1]<$f[3]) { $minali=$f[1]; } else { $minali=$f[3]; }
		my $olap=$f[5]*100/$minali;		#-- Percentage of the query covered by the hit
		next if($olap<$minolap);	#-- Partial hits are not allowed
		$currorf=$f[0];
		my @c=split(/\|/,$f[2]);
		my $khit=$c[1];			#-- This is the KEGG for the hit
		if(!$kegg{$currorf}{besthit}) {
			$kegg{$currorf}{besthit}=$khit;	#-- If it is the first hit, then it is best hit
			# print outfile2 "$currorf\t$kegg{$currorf}\n";
			}
		next if($count{$khit}>=$maxhits);	#-- If we have already $maxhits hits for that KEGG, skip this hit
		$count{$khit}++;
		$accum{$khit}+=$f[$#f];
	     }
	close infile2;

	#-- We need to proccess also the last ORF in the file

	foreach my $k(keys %accum) { $accum{$k}/=$count{$k}; }
	my @list=sort { $accum{$b}<=>$accum{$a}; } keys %accum;
	$score1=$accum{$list[0]}; 
	if($score1) {
 		if($list[1]) { $score2=$accum{$list[1]} } else { $score2=0; } 	#-- Score of the 2nd best
		$dif=(($score1-$score2)/$score1); 
		if($dif>$mindif) { $kegg{$currorf}{bestaver}=$list[0]; }
		print outfile2 "$currorf\t$kegg{$currorf}{besthit}\t$kegg{$currorf}{bestaver}\n";
		}
	close outfile2;
	}		#-- END of COG assignment
	    

#----------------------------------- PFAM assignment -------------------------------------

if(!$nopfam) {

	#-- Read the Pfam data for the pfam.dat file

	my(%pfamname,%hits);
	open(infile3,$pfamlist) || die;
	while(<infile3>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		$pfamname{$k[0]}=$k[2];
		}
	close infile3;	   

	#-- We start reading the hmmer results

	open(outfile3,">$fun3pfam") || warn "Cannot open Pfam output file $fun3pfam\n";
	print outfile3 "# Created by $0, ",scalar localtime,"\n";
	open(infile4,$pfamhmmer);
	while(<infile4>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\s+/,$_);
		my $pfam=$k[4];
		$pfam=~s/\.\d+//;
		$hits{$k[0]}{$pfam}=1;		#-- Simply store the pfam(s) for that ORF
		}
	close infile4;

	foreach my $y(sort keys %hits) {
		print outfile3 "$y";
		my $string="\t";
		foreach my $n(sort keys %{ $hits{$y} }) { $string.="$n [$pfamname{$n}];"; }
		chop $string;
		print outfile3 "$string\n";
		}
			      
	close outfile3;
	}		#-- END of Pfam assignment
