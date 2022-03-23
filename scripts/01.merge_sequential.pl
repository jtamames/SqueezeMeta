#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Merges individual assemblies using minimus2, for merged mode. It also uses cd-hit for excluding identical contigs

use strict;
use Cwd;
use lib "."; 

my $pwd=cwd();

$|=1;

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($installpath,$scriptdir,$resultpath,$interdir,$tempdir,$cdhit_soft,$extassembly,$minimus2_soft,$toamos_soft,$prinseq_soft,$numthreads,$methodsfile,$syslogfile,$force_overwrite);

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";
open(out_tr,">$tempdir/merge.order");

#-- Merges the assemblies in a single dataset

my $finalcontigs="$resultpath/01.$project.fasta";
my($ecode,$command,$mergestep,$merged);
opendir(indir0,$interdir);
my @indassemblies=grep(/01.*fasta$/,readdir indir0);
my $numassem=$#indassemblies+1;
closedir indir0;

if($extassembly) { 
	print "  External assembly provided: $extassembly. Overriding merging\n";
	system("cp $extassembly $finalcontigs");
	}
else {
	print "  Starting assembly merge\n";
	if(-e "$tempdir/mergelog") { system("rm $tempdir/mergelog > /dev/null 2>&1"); }
	system("rm $interdir/merged* > /dev/null 2>&1");
	while($numassem>1) {
		$mergestep++;	
		if(-e $finalcontigs) { system("rm $finalcontigs > /dev/null 2>&1"); }
		$merged="$interdir/merged_$mergestep.$project.fasta";
                my $command="$installpath/lib/SqueezeMeta/kmerdist.pl $projectdir $mergestep >> $syslogfile 2>&1";
		print outsyslog "Calculating distances between metagenomes: $command\n";
		my $ecode=system($command);
		if($ecode!=0) { die "Error running command:    $command"; }
		open(infile0,"$tempdir/$project.2merge") || die;
		$_=<infile0>;
		chomp;
		my($sample1,$sample2,$mdist)=split(/\t/,$_);
		$sample1.=".fasta";
		$sample2.=".fasta";
		close infile0;
		my $md=sprintf ('%.2f', $mdist);
		print "  MERGE $mergestep, $sample1 and $sample2 (dist $md) -> merged_$mergestep.$project.fasta\n";
		print out_tr "MERGE $mergestep, $sample1 and $sample2 ($mdist)\n";
		$command="cat $interdir/$sample1  $interdir/$sample2 > $merged";
		print outsyslog "Merging $sample1 and $sample2 ($mdist): $command\n";
		system $command;
		if(-z $merged) { die "$merged is empty\n"; }
		

			#-- Uses cd-hit to identify and remove contigs contained in others

			my $merged_clustered="$tempdir/mergedassemblies.$project.99.fasta";
			$command="$cdhit_soft -i $merged -o $merged_clustered -T $numthreads -M 0 -c 0.99 -d 100 -aS 0.9 >> $syslogfile";
			print "  Running cd-hit-est\n";
			print outsyslog "Running cd-hit-est: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $merged_clustered) { die "$merged_clustered is empty\n"; }

			#-- Uses Amos to chage format to afg (for minimus2)

			my $afg_format="$tempdir/mergedassemblies.$project.afg";
			$command="$toamos_soft -s $merged_clustered -o $afg_format";
			print "  Transforming to afg format\n";
			print outsyslog "Transforming to afg format: $command\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Now separate reference and query
		
			parseafg($afg_format);
		
			#-- Uses minimus2 to assemble overlapping contigs

			$command="$minimus2_soft\_mod $tempdir/mergedassemblies.$project -D OVERLAP=100 -D MINID=95 -D THREADS=$numthreads > /dev/null 2>&1";
			print outsyslog "Merging with minimus2: $command\n";
			print "  Merging with minimus2\n";
			$ecode = system $command;
			if($ecode!=0) { die "Error running command:    $command"; }
			if(-z $afg_format) { die "$afg_format is empty\n"; }

			#-- Create the final result (overlapping contigs plus singletons)

			$command="cat $tempdir/mergedassemblies.$project.fasta $tempdir/mergedassemblies.$project.singletons.seq > $merged.prov";
			print outsyslog "Joining with singletons: $command\n";
			system($command);
			my $contignm;
			open(outfile0,">$merged") || die;
			open(infile0,"$merged.prov") || die;
			while(<infile0>) {
				chomp;
				if($_=~/^\>/) {
					$contignm++;
					my $newname="Merged$mergestep\_$contignm";
					print outfile0 ">$newname\n"; 
					}
				else { print outfile0 "$_\n"; }
				}
			close infile0;
			close outfile0;
			system("rm $merged.prov > /dev/null 2>&1");
			
			$numassem--;
			open(mlog0,">>$tempdir/mergelog") || die;
			print mlog0 "$sample1\n$sample2\n";
			close mlog0;
			#system("mv $interdir/$sample1 $interdir/$sample1.orig");
			#system("mv $interdir/$sample2 $interdir/$sample2.orig");
			#opendir(indir0,$interdir);
			#@indassemblies=grep(/fasta$/,readdir indir0);
			#$numassem=$#indassemblies+1;
			}
				

		#-- Remove files from temp

		 system("rm -r $tempdir/mergedassemblies* > /dev/null 2>&1");
	}

system("mv $merged $finalcontigs");
print outmet "Totally overlapping contigs were removed using cd-hit(Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
print outmet "Contigs were merged using Minimus2(Treangen et al 2011, Curr Protoc Bioinfomatics 11)\n";


#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $finalcontigs -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
print outsyslog "Using prinseq for contig statistics: $command\n";
$ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	

#-- Count length of contigs (needed later)

my $contigslen="$interdir/01.$project.lon";
print "  Counting lengths of contigs\n";
open(outfile1,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile1,$finalcontigs) || die "Can't open $finalcontigs\n";
my($thisname,$contigname,$numc,$seq);
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) {
	$numc++;
	$thisname=$1;
	if($contigname) {
		my $len=length $seq;
		print outfile1 "$contigname\t$len\n"; 
		}
	$seq="";
	$contigname=$thisname;
	}
	else { $seq.=$_; }
}
close infile1;
if($contigname) { my $len=length $seq; print outfile1 "$contigname\t$len\n"; }
close outfile1;
close out_tr;
close outmet;
print "  Number of contigs: $numc\n";

sub parseafg {				
	my $inafg=shift;
	my($inred,$inpos,$accum,$ts);
	my(%order,%samples);
	open(inp1,$inafg) || die;
	while(<inp1>) {
	chomp;
	next if !$_;
	if($_=~/\{RED/) { $inred=1; }
	elsif($_=~/\}/) { $inred=0; }
	next if(!$inred);
	if($_=~/eid\:(.*)/) {
		$inpos++;
		my @m=split(/\_/,$1);
		if($m[0]=~/^Merged/) { $ts=$m[0]; }
		else {
			shift @m; shift @m;
			 $ts=join("_",@m);
			 $ts=$m[$#m];   #-- Last field of contig name contains the sample ID. Previous solution worked for megahit but not for SPAdes. This one can get trouble if sample names contains "_"
			}
		$order{$inpos}=$ts;
		$samples{$ts}=1;
		}
	}

	close inp1;
	foreach my $p(sort keys %samples) { $accum++; $samples{$p}=$accum; }

	open(out1,">$inafg.1") || die;
	open(out2,">$inafg.2") || die;
	my($headsw,$intofrg,$numfrg,$intored,$numred,$samp,$tofile);
	open(inp2,$inafg) || die;
	while(<inp2>) {
		chomp;
		next if !$_;
		if($_=~/\{FRG/) { $headsw=1; $intofrg=1; $numfrg++; }
		if($_=~/\{RED/) { $headsw=1; $intofrg=0; $intored=1; $numred++; }
		if(!$headsw) {
			print out1 "$_\n";
			print out2 "$_\n";
			next;
			}
		if($intofrg) {	
			$samp=$order{$numfrg};
			$tofile=$samples{$samp}; 
			if($tofile==1) { print out1 "$_\n"; } else { print out2 "$_\n"; } 
			}
		elsif($intored) {
			$samp=$order{$numred};
			$tofile=$samples{$samp}; 
			if($tofile==1) { print out1 "$_\n"; } else { print out2 "$_\n"; } 
			}
		
		}
	close inp2;
	}
