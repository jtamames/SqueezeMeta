#!/usr/bin/env perl
use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

my $pwd=cwd();

if($ARGV[0] eq "-h") {
 print "\nUsage: $0 <project> <KEGG pathway> <mode>\n"; 
 print "  Mode can be \"single\" (default) or \"multiple\".\n  Single will create a map for each sample.\n  Multiple will plot all samples in one map\n";
 print "\n  Information of Pathway codes can be found in data/kegg_pathlist.txt\n\n";
 die;
 }
my $projectpath=$ARGV[0];
my $path=$ARGV[1];
my $mode=$ARGV[2];
if(!$projectpath) { die "Please provide a valid project name or project path (use -h for help)\n"; }
if(!$path) { die "Please provide a KEGG pathway code (use -h for help)\n"; }
if($mode && (($mode!~/single/i) && ($mode!~/multiple/i))) { die "Invalid mode, please use single or multiple (use -h for help)\n"; };
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

our($datapath,$resultpath);

if($mode=~/multiple/i) { $mode="multiple"; } else { $mode="single"; }

my(%samples,%insample);

open(in,"$resultpath/12.$project.kegg.funcover") || die;
while(<in>) {
	chomp;
	next if($_!~/^K/);
	my @t=split(/\t/,$_);
	$insample{$t[0]}{$t[1]}=$t[6];
	$samples{$t[1]}=1;
	}
close in;

my @l=keys %samples;
my $numsamples=$#l+1;
open(out,">gtable.txt") || die;
print out "Gene";
foreach my $p1(sort keys %samples) { print out " $p1"; }
print out "\n";
foreach my $k(sort keys %insample) {
	print out "$k";
	foreach my $p1(sort keys %samples) { 
		my $inum=$insample{$k}{$p1} || "0";
		print out " $inum"; 
		}
	print out "\n";
	}
close out;

if($mode eq "multiple") {
	open(outR,">script.R") || die;
	print outR "library('pathview')\n";
	print outR "myd<-read.table(\"gtable.txt\",row.names=1,header=TRUE)\n";
	#print outR "pv.out <- pathview(gene.data = as.matrix(myd)[, 1], gene.idtype=\"kegg\", pathway.id = \"$path\",species = \"ko\", out.suffix = \"$project\", kegg.native = T, both.dirs=list(gene=FALSE,cpd=FALSE),limit=list(gene=100, cpd=100))\n";

	print outR "pv.out <- pathview(gene.data = as.matrix(myd)[, 1:$numsamples], gene.idtype=\"kegg\", pathway.id = \"$path\",species = \"ko\", limit =list(gene=10), discrete=list(gene=F), bins=list(gene=16), na.col=\"white\", out.suffix = \"$project\", kegg.native = T, same.layer=FALSE,both.dirs=list(gene=FALSE),mid=\"white\")\n";

	close outR;
	system("R CMD BATCH script.R");
	}
	
elsif($mode eq "single") {
	for(my $n=1; $n<=$numsamples; $n++) {
		open(outR,">script.R") || die;
		print outR "library('pathview')\n";
		print outR "myd<-read.table(\"gtable.txt\",row.names=1,header=TRUE)\n";
		print outR "cname=colnames(myd)[$n]\n";
		print outR "cpath=paste0(\"$project.\",cname)\n";
		print outR "pv.out <- pathview(gene.data = as.matrix(myd)[, $n], gene.idtype=\"kegg\", pathway.id = \"$path\",species = \"ko\", limit =list(gene=10), discrete=list(gene=F), bins=list(gene=16), na.col=\"white\", out.suffix = cpath, kegg.native = T, same.layer=FALSE,both.dirs=list(gene=FALSE),mid=\"white\")\n";
		close outR;
		system("R CMD BATCH script.R");
		}
	}
