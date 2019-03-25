#!/usr/bin/perl

use Cwd;
use Getopt::Long;
use strict;

$|=1;

my($complete_cutoff,$contamination_cutoff,$funclass,$reqfunctions,$hel);

#-- Define help text

my $helptext = <<END_MESSAGE;
Usage: makeitolplot.pl <options> project name

Options:

   -completion [number]: Select only bins with completion ABOVE that threshold (Default: 30)
   -contamination [number]: Select only bins with contamination BELOW that threshold (Default: 100)
   -classification [metacyc|kegg]: Functional classification to use (Default:metacyc)
   -functions [file]: File containing the name of the fucntions to be considered
   -h: This help
     
END_MESSAGE

my $result = GetOptions ("completion=i" => \$complete_cutoff,
                     "contamination=i" => \$contamination_cutoff,
                     "classification=s" => \$funclass,
		     "functions=s" => \$reqfunctions,
		     "h" => \$hel
 		    );

if($hel) { die "$helptext\n"; } 

my $pwd=cwd();
my $project=pop @ARGV;
$project=~s/\/$//;
if(!$project) { die "Please provide project name\nUsage: makeitolplot.pl <options> project name\n"; }

do "$project/SqueezeMeta_conf.pl";

our($extdatapath,$contigsinbins,$mergedfile,$aafile,$tempdir,$resultpath,$minpath_soft,$bintable,$extpath,%bindirs,%dasdir);

if(!$complete_cutoff) { $complete_cutoff=30; }		#-- Do not consider bins below this level of completion
if(!$contamination_cutoff) { $contamination_cutoff=100; }		#-- Do not consider bins above this level of contamination
if(!$funclass) { $funclass="metacyc"; }
my $numtaxalabels=4;

my $dirbin=$dasdir{DASTool};

my @colors=("#fc05ea","#0000ff","#ff0000","#00ff00","#7a4304","#f77300","#00adf7","#04721c","#7a30db","#c4c107");	#-- Magenta, Blue, Red, Green, brown, Orange, ligth blue, dark green, yellow
my @colors2=("#f77300","#00adf7","#04721c","#7a30db","#c4c107");	#-- Orange, ligth blue, dark green, yellow

print "Reading bins in $bintable...\n";

my($header,$totalbins,$numbins,);
my(%newnames,%complete,%tax,%store,%samcode,%contamination,%phylo,%countphylo,%requested,%fun,%allfun,%color,%contigs,%bins,%seqs,%distance);
my @header;
open(infile1,$bintable) || die "Cannot open $bintable";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	if(!$header) { 
		$header=$_; 
		@header=split(/\t/,$_);
		next; 
		}
	my @k=split(/\t/,$_);
	next if($k[1] ne "DASTool");
	$totalbins++;
	next if($k[8]<$complete_cutoff);
	next if($k[9]>$contamination_cutoff);
	$numbins++;
	chop $k[2]; 
	my @mtx=split(/\;/,$k[2]);
	my ($lrank,$ltax)=split(/\:/,$mtx[$#mtx]);
	my $binname="$k[0]_$ltax";
	$binname=~s/\s+/\_/g;
	$binname=~s/sp\./sp/g;
	$newnames{$k[0]}=$binname;
	$complete{$binname}=$k[8];
	$tax{$binname}=$k[2];
	for(my $pos=0; $pos<=$#k; $pos++) {
		if($header[$pos]=~/RPKM (.*)/) {
			my $corrdata=$1;
			$store{$binname}{$corrdata}=$k[$pos];
			$samcode{$corrdata}=1;
			}
		}
	$contamination{$binname}=$k[9];

	if($k[2]=~/phylum\:([^;]+)/) { $phylo{$binname}{phylum}=$1; $countphylo{phylum}{$1}++; }
	if($k[2]=~/genus\:([^;]+)/) { $phylo{$binname}{genus}=$1; $countphylo{genus}{$1}++;}
	}
close infile1;
print "Found $totalbins bins\nWorking with $numbins bins with more than $complete_cutoff % completion and less than $contamination_cutoff % contamination\n";


if($reqfunctions) {
	open(infile1,$reqfunctions) || die "Cannot open requested functions file $reqfunctions\n";
	print "Reading requested functions from $reqfunctions\n";
	while(<infile1>) {
		chomp;
		$requested{$_}=1;
		}
	close infile1;
	}

my $funfile="$resultpath/21.$project.$funclass.pathways";
print "Reading functions in $funfile\n";
open(infile1,$funfile) || die;
$_=<infile1>;
my $headerf=$_; 
chomp $headerf;
my @headerf=split(/\t/,$_);
if(!$reqfunctions) {
	for(my $pos=3; $pos<=$#headerf; $pos++) { $requested{$headerf[$pos]}=1; }
	}
my $pathnum=<infile1>;
chomp $pathnum;
my @pathnum=split(/\t/,$pathnum); 
while(<infile1>) {
	chomp;
	my @k=split(/\t/,$_);
	my $biname=$newnames{$k[0]};
	next if(!$tax{$biname});
	for(my $pos=3; $pos<=$#k; $pos++) {
		my $refun=$headerf[$pos];			#-- Full label
		my @lifun=split(/\;/,$headerf[$pos]);
		my $shortfun=$lifun[$#lifun];				#-- Just pathway as a label
		foreach my $u(keys %requested) {
			if($refun=~/$u/i) {
				my $tratio;
				if($k[$pos] eq "NF") { $tratio=0; } else { $tratio=$k[$pos]/$pathnum[$pos]; }
				$fun{$biname}{$shortfun}=$tratio; 
				$allfun{$shortfun}=1;
				}
			}
		}
	}
close infile1;
				


my $cphy=0;
foreach my $phyl(sort { $countphylo{phylum}{$b}<=>$countphylo{phylum}{$a}; } keys %{ $countphylo{phylum} }) {
	$cphy++;
	if($cphy>$numtaxalabels) { delete $countphylo{phylum}{$phyl}; } else { $color{phylum}{$phyl}=$colors[$cphy-1]; }
	}
	
my $cphy=0;
foreach my $phyl(sort { $countphylo{genus}{$b}<=>$countphylo{genus}{$a}; } keys %{ $countphylo{genus} }) {
	$cphy++;
	if($cphy>$numtaxalabels) { delete $countphylo{genus}{$phyl}; } else { $color{genus}{$phyl}=$colors2[$cphy-1]; }
	}

my $distfile="$resultpath/maxbin/aai/dist.out";
my $treefile="$extpath/$project.nw";
my ($contig);

if(-e $treefile) {}
else {						#-- Run compareM only if we have not done it before
	open(infile2,$contigsinbins) || die "Cannot open contigsinbins file $contigsinbins\n";
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		my $bin=$k[2];
		my $binname=$newnames{$bin};
		if($complete{$binname}) { 
			$contigs{$k[0]}{$binname}=1; 
			$bins{$binname}{$k[0]}=1;
			}
		}
	close infile2;

	print "Reading sequences in $aafile\n";
	open(infile3,$aafile) || die "Cannot open amino acid fasta file $aafile\n";
	while(<infile3>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if($_=~/^\>([^ ]+)/) { 
			my $orf=$1; 
			my @sf=split(/\_/,$orf);
			my $ipos=pop @sf;
			$contig=join("_",@sf);
			}
		$_=~s/\*//g;
		if($contigs{$contig}) { $seqs{$contig}.="$_\n"; }
		}
	close infile3;

	my $tempoutdir="$tempdir/comparem";
	if(-d $tempoutdir) { system("rm -r $tempoutdir"); }
	system("mkdir $tempoutdir"); 
	foreach my $bin(sort keys %bins) {
		next if($complete{$bin}<$complete_cutoff);
		print "Getting sequences for bin $bin                 \r";
		my $outfile="$tempoutdir/$bin.faa";
		open(outfile1,">$outfile") || die "Cannot open $outfile\n";
		foreach my $ct(sort keys %{ $bins{$bin} }) { 
			$seqs{$ct}=~s/\>/\>$bin\_/g;
	 		print outfile1 "$seqs{$ct}\n"; 
			}
		close outfile1;
		}

	my $command="comparem aai_wf --proteins $tempoutdir -x faa $resultpath/maxbin";	
	print "\nRunning CompareM: $command\n";
	system($command);

	open(infile4,"$resultpath/maxbin/aai/aai_summary.tsv") || die;
	$_=<infile4>;
	while(<infile4>) {
		chomp;
		next if !$_;
		my @k=split(/\t/,$_);
		next if(($k[1]<1500) || ($k[3]<1500));
		my $dist=100-$k[5];
		$distance{$k[0]}{$k[2]}=$dist;
		$distance{$k[2]}{$k[0]}=$dist;
		}
	close infile4;

	open(outfile2,">$distfile") || die;
	print outfile2 "Genome A\tGenome B\tDissimilarity\n";

	foreach my $p1(sort keys %distance) { 
		foreach my $p2(sort keys %distance) { 
			next if($p1 ge $p2);
			my $dista=$distance{$p1}{$p2}; 
			print outfile2 "$p1\t$p2\t$dista\n";
			}
		}
	close outfile2;

	$command="comparem hclust $distfile $treefile";
	system $command;
	print "Tree created in $treefile\n";
	}

my $out2="$extpath/heatmap_abund.txt";
open(out1,">$out2") || die "Cannot open $out2\n";
print out1 "DATASET_HEATMAP\nSEPARATOR TAB\nDATASET_LABEL\tabundances\nCOLOR\t#ff0000\n";
print out1 "LEGEND_TITLE\tPhyla (Tree branches)\nLEGEND_SHAPES";
for(my $num=1; $num<=$numtaxalabels; $num++) { print out1 "\t2"; }
print out1 "\nLEGEND_COLORS";
for(my $num=1; $num<=$numtaxalabels; $num++) {
	my $tcolor=$colors[$num-1];
	print out1 "\t$tcolor"; 
	}
#foreach my $cl(@colors) { print out1 "\t$cl"; } 
print out1 "\nLEGEND_LABELS";
foreach my $phyl(sort { $countphylo{phylum}{$b}<=>$countphylo{phylum}{$a}; } keys %{ $countphylo{phylum} }) { print out1 "\t$phyl"; }
print out1 "\n";
print out1 "FIELD_LABELS";
foreach my $sam(sort keys %samcode) { print out1 "\t$sam"; }
print out1 "\n";
print out1 "DATA\n";
foreach my $meta(keys %fun) { 
	print out1 "$meta";
	foreach my $sam(sort keys %samcode) {
		my $value= $store{$meta}{$sam};
		if($value>20) { $value=20; }
		print out1 "\t$value"; 
		}
	print out1 "\n";
	}
close out1;	
print "Dataset 1 (abundances) created in $out2\n";

$out2="$extpath/heatmap_pathways.txt";
open(out1,">$out2") || die;
print out1 "DATASET_HEATMAP\nSEPARATOR TAB\nDATASET_LABEL\tpathways\nCOLOR\t#ff00ff\n";
print out1 "LEGEND_TITLE\tPhyla (Tree branches)\nLEGEND_SHAPES";
for(my $num=1; $num<=$numtaxalabels; $num++) { print out1 "\t2"; }
print out1 "\nLEGEND_COLORS";
for(my $num=1; $num<=$numtaxalabels; $num++) {
	my $tcolor=$colors[$num-1];
	print out1 "\t$tcolor"; 
	}
#foreach my $cl(@colors) { print out1 "\t$cl"; } 
print out1 "\nLEGEND_LABELS";
foreach my $phyl(sort { $countphylo{phylum}{$b}<=>$countphylo{phylum}{$a}; } keys %{ $countphylo{phylum} }) { print out1 "\t$phyl"; }
print out1 "\n";
print out1 "FIELD_LABELS";
foreach my $fu(sort keys %allfun) { print out1 "\t$fu"; }
print out1 "\n";
print out1 "DATA\n";
foreach my $meta(keys %fun) { 
	print out1 "$meta";
	foreach my $fu(sort keys %allfun) {
		my $value= $fun{$meta}{$fu};
		printf out1 "\t%.3f",$value; 
		}
	print out1 "\n";
	}
close out1;	
print "Dataset 2 (pathways) created in $out2\n";

	
my $out3="$extpath/heatmap_completion.txt";
open(out2,">$out3") || die;
print out2 "DATASET_HEATMAP\nSEPARATOR TAB\nDATASET_LABEL\tcompletion\nCOLOR\t#ffff00\n";
#print out2 "LEGEND_TITLE\tGenera (Tree labels)\nLEGEND_SHAPES\t2\t2\t2\t2\nLEGEND_COLORS";
#foreach my $cl(@colors2) { print out2 "\t$cl"; } 
#print out2 "\nLEGEND_LABELS";
#foreach my $phyl(sort { $countphylo{genus}{$b}<=>$countphylo{genus}{$a}; } keys %{ $countphylo{genus} }) { print out2 "\t$phyl"; }
#print out2 "\n";
print out2 "FIELD_LABELS\tCompletion\n";
print out2 "DATA\n";
foreach my $meta(keys %complete) { 
	my $value=$complete{$meta};
	if(!$value) { $value="0"; }
	next if($value<$complete_cutoff);
	print out2 "$meta\t$value\n"; }
close out2;
print "Dataset 3 (Completion) created in $out3\n";

my $out4="$extpath/labels.txt";
open(out3,">$out4") || die;
print out3 "TREE_COLORS\n";
print out3 "SEPARATOR TAB\n";
print out3 "DATA\n";
foreach my $binp(keys %phylo) {  
	my $phylum=$phylo{$binp}{phylum};
	if($color{phylum}{$phylum}) { print out3 "$binp\tbranch\t$color{phylum}{$phylum}\tnormal\n"; }
	my $genus=$phylo{$binp}{genus};
	if($color{genus}{$genus}) { print out3 "$binp\tlabel\t$color{genus}{$genus}\tnormal\n"; }
	}
close out3;
print "Dataset 4 (Labels) created in $out4\n";
my $tarfile="$extpath/plot.tar";
my $command="cd $extpath; tar cvf $tarfile $project.nw heatmap_abund.txt heatmap_completion.txt labels.txt > /dev/null";
system $command;
print "tar file created in $tarfile\n";

print "Now upload tree and datasets in https://itol.embl.de\n";
