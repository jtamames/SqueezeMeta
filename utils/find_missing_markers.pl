#!/usr/bin/perl

#-- Looks for contigs contaiing missing markers for the specified bin, or all of them if none given

local $; = '#';
use strict;
use Cwd;
use lib ".";
use Getopt::Long;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(!$project) { die "Usage: $0 <project> [bin name]\n"; }
my $binreq=$ARGV[1];
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project name ok?\n"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

my $verbose=0;		# For debugging purposes

my $mode="strict";	# "strict" (taxonomy of the contig must match exactly the one for the bin) or "relaxed" (taxonomy of the contig must not contradict the one for the bin)

our($alllog,%bindir,%dasdir,$installpath,$checkm_soft,$numthreads,$datapath,$tempdir,$taxlist,$mergedfile,$contigtable,$bindir,$databasepath,$contigsinbins,$contigtable);
my %branks=('k','domain','p','phylum','c','class','o','order','f','family','g','genus','s','species');
my $markerdir="$datapath/checkm_markers";
my $checktemp="$tempdir/checkm_batch";
my $tempc="$tempdir/checkm_prov.txt";


my $bindir=$dasdir{DASTool};
my @binlist;
my(%centroid,%tigrfam);

my $outfile="$bindir/missingmarkers.txt";
open(outfile1,">$outfile") || die;
print outfile1 "# Created by $0 for project $project, ",scalar localtime,"\n";

if($binreq) { push(@binlist,$binreq); }
else {
	opendir(indir1,$bindir) || die "Cannot open binning directory $bindir\n";
	@binlist=(grep/\.fa$|\.fasta$/,readdir indir1);
	closedir indir1;
	}

my $tigr2pfam="$databasepath/pfam/tigrfam2pfam.tsv";
open(infile0,$tigr2pfam) || die "Cannot open tigr2pfam file $tigr2pfam\n";
while(<infile0>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my($pf,$tg)=split(/\t/,$_);
	$pf=~s/pfam/PF/;
	$tigrfam{$tg}=$pf;
	}
close infile0;

	
open(infile0,$mergedfile) || die;
my($pfamf,$header);
while(<infile0>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	if(!$header) {
	$header=$_;
	my @k=split(/\t/,$header);
	for(my $pos=0; $pos<=$#k; $pos++) {
		if($k[$pos] eq "PFAM") { $pfamf=1; last; }
		}
	}
	last if($header);
	}
close infile0;
if(!$pfamf) { die "Your gene table $mergedfile does not have PFAM annotations, this is required for finding missing markers\n"; }
	
 	#-- Read NCBI's taxonomy 

my %tax;
open(infile1,$taxlist) || die "Can't find taxonomy list $taxlist\n";
print "Reading $taxlist\n";
while(<infile1>) {
 	chomp;
 	next if !$_;
 	 my @t=split(/\t/,$_);
 	 $tax{$t[1]}=$t[2];
 	 }
 close infile1;
		

my $skip;
my($taxbin,$size,$chim,$chimlev,$thres)="";
foreach my $bin(sort @binlist) {
	my $fastaname="$bindir/$bin";
	my $binname="$bindir/$bin\.tax";
	next if (!(-e $binname));
	my(%consensus,%alltaxa);
	my @misslist;
	$thres=0;
	open(infile1,$binname) || die "Can't open $binname\n";
	print "\nFound bin $bin: ";
	while(<infile1>) { 
		chomp;
		if($_=~/Consensus/) {
 			($taxbin,$size,$chim,$chimlev)=split(/\t/,$_);
 			$taxbin=~s/Consensus\: //;
 			$size=~s/Total size\: //g;
 			$consensus{$bin}=$taxbin;
 		 	my @k=split(/\;/,$taxbin);
		 	if($taxbin!~/k\_/) { print "No consensus tax. Skipping\n"; $skip=1; last; }
 		 	print "$taxbin\n";
			
			 		 #-- We store the full taxonomy for the bin because not all taxa have checkm markers
 	 
 			 foreach my $ftax(reverse @k) { 
 				 my($ntax,$rank);
 				 if($ftax!~/\_/) { $ntax=$ftax; } else { ($rank,$ntax)=split(/\_/,$ftax); }
 				 $ntax=~s/unclassified //gi;
 				 $ntax=~s/ \<.*\>//gi; 
 				 if($tax{$ntax} && ($rank ne "n") && ($rank ne "s")) { 
 				 push( @{ $alltaxa{$binname} },"$branks{$rank}\_$ntax");
 			 	#   print "$m\t$ntax\t$tax{$ntax}\n";
 					 }
 				 }
			}
		}
	close infile1;
	next if($skip);
	my $genetable=$mergedfile;
	if(!$genetable) { die "Missing gene table\n"; }
	if(!$contigtable) { die "Missing contig table\n"; }

 # my $intdist=centroid($bin);
 # print "\nCentroid: $intdist";
 my $inloop=1;
 
 #-- We will find the deepest taxa with checkm markers
 while($inloop) {  
 	 my $taxf=shift(@{ $alltaxa{$binname}  }); 
 	 if(!$taxf) { last; $inloop=0; }
 	 my($rank,$tax)=split(/\_/,$taxf);
 	 $tax=~s/ \<.*//g;
 	 $tax=~s/\s+/\_/g;
 	 # my $rank=$equival{$grank};
 	 if($rank eq "superkingdom") { $rank="domain"; }
 	 print " Using profile for $rank rank : $tax\n";   
 	 my $marker="$markerdir/$tax.ms"; 

 	 #-- Use already existing tax profile or create it

 	 if(-e $marker) {} else { 
 		 my $command="export PATH=\"$installpath/bin/pplacer\":\$PATH; $checkm_soft taxon_set $rank $tax $marker > /dev/null 2>&1"; #Override $PATH for external dependencies of checkm. (FPS).
 		 my $ecode = system $command;
 		 if($ecode!=0) { die "Error running command:	$command"; }
 		 }

 	 #-- If it was not possible to create the profile, go for the upper rank
 	 
 	 if(-e $marker) {} else { next; }

 	 my $fastafile=$binname;
 	 $fastafile=~s/\.tax//;
 	 $fastafile=~s/.*\///;
 	 # print ">>> $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null\n";
 	 # system("$checkm_soft analyze -t $numthreads -x $bins{$binname} $marker $bindir $checktemp > /dev/null");
 	 my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null 2>&1";
 	 # print "$command\n";
 	 my $ecode = system $command;
 	 if($ecode!=0) { die "Error running command:	$command"; }

 	 my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft qa -t $numthreads $marker $checktemp -f $tempc > /dev/null 2>&1"; #Override $PATH for external dependencies of checkm. (FPS).
 	 # print "$command\n";
	 my $ecode = system $command;
 	 if($ecode!=0) { die "Error running command:	$command"; }
 	 $inloop=0;
 	 }


	my $ablist;
	my(%ml,%incontig);
	open(infile2,"$checktemp/storage/bin_stats_ext.tsv") || die;
	while(<infile2>) {
		chomp;
		next if !$_;
		if($_=~/GCN0\'\: \[(.*?)\]/) { $ablist=$1; }
		}
	close infile2;

	$ablist=~s/\s+//g;
	$ablist=~s/\'//g;
	@misslist=split(/\,/,$ablist);
	my $mk;
	foreach my $tmark(@misslist) {
		if($tmark=~/^TIGR/) { $mk=$tigrfam{$tmark}; } else { $mk=$tmark; }
		$ml{$mk}=$tmark;
		}
	print "Missing markers: @misslist\n";

	open(infile3,$genetable) || die;
	my($pfampos,$header);
	while(<infile3>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if(!$header) {
		$header=$_;
		my @k=split(/\t/,$header);
		for(my $pos=0; $pos<=$#k; $pos++) {
			if($k[$pos] eq "PFAM") { $pfampos=$pos; last; }
			}
		next;
		}	
      		my @k=split(/\t/,$_);
		# next if($k[5]!~/$rank\:$lookingfor/);
		my $pfam=$k[$pfampos];
		next if(!$pfam);
		my @pf=($pfam=~/PF\d+/g);
		my $idex;
		foreach my $pg(@pf) {
			if($ml{$pg}) { 
			# print "$k[0]\t$k[1]\t$pg\t$k[5]\n";
			if($pg eq $ml{$pg}) { $idex=$pg; } else { $idex=$ml{$pg}; }
			$incontig{$k[1]}{$idex}=$k[0]; 
			$thres=1;
			}
		}
	}
close infile3;

next if(!$thres);
my(%strings,%strnum);
print outfile1 "$bin; Consensus: $taxbin\n"; 
print outfile1 "Missing markers for $bin: @misslist\n";
print outfile1 "Bin\tContig\tMarker\tGene\tContig tax\n";
open(infile4,$contigtable);
while(<infile4>) {
        chomp;
        next if(!$_ || ($_=~/^\#/));
        my @k=split(/\t/,$_);
	my $contigid=$k[0];
	my $contigtax=$k[1];
	foreach my $igen(keys %{ $incontig{$contigid} }) {
		next if(($mode eq "strict") && ($contigtax!~/$taxbin/));
		next if(($mode eq "relaxed") && ($taxbin!~/$contigtax/));
		#my $contigdist=dist_to_centroid($contigid);
		#my $normdist=$contigdist/$intdist;
		# print outfile1 "$bin\t$contigid\t$igen\t$incontig{$contigid}{$igen}\t$contigtax\n";
		$strings{$contigid}.="$bin\t$contigid\t$igen\t$incontig{$contigid}{$igen}\t$contigtax\n";
		$strnum{$contigid}++;
		}
	}
foreach my $p(sort { $strnum{$b}<=>$strnum{$a}; } keys %strnum) { print outfile1 $strings{$p}; }
close infile4;
print outfile1 "\n";

	}
close outfile1;



sub centroid {
	my $bin=shift;
	$bin=~s/\.fa$|\.fasta$//g;
	my(%cinbin,%tpm,%gc,%tpmsets,%intrad,@headerp,$headerp,$dimensions,$intradist,@e);
	%centroid=();
	open(infile5,$contigtable) || die;
	while(<infile5>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if(!$headerp) { 
			$headerp=$_; 
			@headerp=split(/\t/,$headerp); 
			for(my $pos=0; $pos<=$#headerp; $pos++) { 
				if($headerp[$pos]=~/^TPM/) {  $tpmsets{$headerp[$pos]}=1; $dimensions++; }
				}
			next; 
			}
		@e=split(/\t/,$_);
		if($e[6]=~/$bin/) { 
			$cinbin{$e[0]}=1; 			
			for(my $pos=0; $pos<=$#e; $pos++) {
				if($headerp[$pos]=~/^TPM/) { $tpm{$e[0]}{$headerp[$pos]}=$e[$pos]; }
				if($headerp[$pos]=~/^GC perc/) { $gc{$e[0]}=$e[$pos]; }
				}
			}
		}
	close infile5;
	foreach my $dp(sort keys %tpmsets) {
		foreach my $points(sort keys %tpm) { $centroid{$dp}+=$tpm{$points}{$dp}; }
		$centroid{$dp}/=$dimensions;	
		# print "$dp: $centroid{$dp} ";		
		}
	# print "\n";

	my $dpp;
	foreach my $dp(sort keys %tpmsets) {
		foreach my $points(sort keys %tpm) { 
			$intrad{$dp}+=($tpm{$points}{$dp}-$centroid{$dp})**2;
			$dpp++; 
			}
		}
	foreach my $r(keys %intrad) { $intradist+=$intrad{$r}; }
	$intradist=(sqrt($intradist))/$dpp;
	return $intradist;	
	}

		
sub dist_to_centroid {
	my $thiscontig=shift;
	my(%tpm,$gc,%tpmsets,%intrad,%distcent,@headerp,$headerp,$dimensions,$intradist,@e,$distaccum);
	open(infile5,$contigtable) || die;
	while(<infile5>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if(!$headerp) { 
			$headerp=$_; 
			@headerp=split(/\t/,$headerp); 
			for(my $pos=0; $pos<=$#headerp; $pos++) { 
				if($headerp[$pos]=~/^TPM/) {  $tpmsets{$headerp[$pos]}=1; $dimensions++; }
				}
			next; 
			}
		@e=split(/\t/,$_);
		if($e[0] eq $thiscontig) { 
			for(my $pos=0; $pos<=$#e; $pos++) {
				if($headerp[$pos]=~/^TPM/) { $tpm{$headerp[$pos]}=$e[$pos]; }
				if($headerp[$pos]=~/^GC perc/) { $gc=$e[$pos]; }
				}
			}
		}
	close infile5;
	foreach my $dp(sort keys %tpmsets) {
		foreach my $points(sort keys %tpm) { $distaccum+=($tpm{$dp}-$centroid{$dp})**2; }
		}
	$distaccum=sqrt $distaccum;
	return $distaccum;
	}
	
	



