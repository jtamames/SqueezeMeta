#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use Linux::MemInfo;
use Term::ANSIColor qw(:constants);
use lib ".";
use strict;

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my $pwd=cwd();
my $utilsdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $utilsdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$utilsdir/..");
my $scriptdir = "$installpath/scripts";
my $auxdir = "$installpath/lib/SQM_reads";

###

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

my $start_run = time();

do "$scriptdir/SqueezeMeta_conf.pl";
do "$scriptdir/parameters.pl";
#-- Configuration variables from conf file
our($databasepath);

my($numthreads,$project,$equivfile,$rawseqs,$miniden,$evalue,$dietext,$blocksize,$currtime,$nocog,$nokegg,$opt_db,$hel,$printversion,$nodiamond,$euknofilter,$methodsfile,$evaluetax4,$minidentax4);

my $helpshort="Usage: SQM_longreads.pl -p <project name> -s <samples file> -f <raw fastq dir> [options]\n";

my $helptext = <<END_MESSAGE;
Usage: SQM_reads.pl -p <project name> -s <samples file> -f <raw fastq dir> [options]

Arguments:

 Mandatory parameters:
   -p: Project name (REQUIRED)
   -s|-samples: Samples file (REQUIRED)
   -f|-seq: Fastq read files' directory (REQUIRED)
   
 Options:
   --nocog: Skip COG assignment (Default: no)
   --nokegg: Skip KEGG assignment (Default: no)
   --nodiamond: Skip Diamond runs, assuming that you already did it (Default: no)
   --euk: Drop identity filters for eukaryotic annotation  (Default: no)
   -extdb <database file>: List of user-provided databases
   -e|-evalue: max evalue for discarding hits for Diamond run  (Default: 1e-03)
   -i|-miniden: minimum identity for the hits (Default: 30)
   -t: Number of threads (Default: 12)
   -b|-block-size: block size for Diamond run against the nr database (Default: 8)
   -v|version: Print version
   -h: this help

END_MESSAGE

my $result = GetOptions ("t=i" => \$numthreads,
                     "p=s" => \$project,
                     "s|samples=s" => \$equivfile,
                     "f|seq=s" => \$rawseqs, 
		     "e|evalue=f" => \$evalue,   
		     "i|miniden=f" => \$miniden,
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nodiamond" => \$nodiamond,   
		     "nodiamond" => \$nodiamond,   
		     "extdb=s" => \$opt_db, 
		     "euk" => \$euknofilter,
                     "b|block_size=i" => \$blocksize,
		     "v|version" => \$printversion,
		     "h" => \$hel
		    );

if(!$numthreads) { $numthreads=12; }
if(!$evalue) { $evalue=0.001; }
if(!$miniden) { $miniden=30; }         #-- Minimum identity for the hit
if(!$euknofilter) { $euknofilter="0"; }
my $querycover=0;	#-- Minimum coverage of hit in query
	
print BOLD "\nSqueezeMeta on Long Reads v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;
if($printversion) { exit; }

if(!$blocksize) {
        print "\nSetting block size for Diamond\n";
        my %mem=get_mem_info;
        my $ram=$mem{"MemAvailable"}/(1024*1024);
        my $ramstr=sprintf('%.2f',$ram);
        my $block_size_set=sprintf('%.1f',$ram/5);
        if($block_size_set>8) { $block_size_set=16; }
        if($block_size_set<1) { $block_size_set=1; }
        print "  AVAILABLE (free) RAM memory: $ramstr Gb\nWe will set Diamond block size to $block_size_set (Gb RAM/5, Max 16). You can override this setting using the -b option when starting the project.\n\n";
        $blocksize=$block_size_set;
        }

if($hel) { die "$helptext\n"; } 
if(!$project) { $dietext.="MISSING ARGUMENT: -p: Project name\n"; }
if(!$rawseqs) { $dietext.="MISSING ARGUMENT: -f|-seq:Read files' directory\n"; }
if(!$equivfile) { $dietext.="MISSING ARGUMENT: -s|-samples: Samples file\n"; }
if($dietext) { print BOLD "$helpshort"; print RESET; print RED; print "$dietext"; print RESET;  die; }

my(%allsamples,%ident,%noassembly,%accum,%totalseqs,%optaccum,%allext);
my($sample,$file,$iden,$mapreq);
tie %allsamples,"Tie::IxHash";

my $verbose=0;
my $nr_db="$databasepath/nr.dmnd";
my $cog_db="$databasepath/eggnog";
my $kegg_db="$databasepath/keggdb";
my $diamond_soft="$installpath/bin/diamond";
my $coglist="$installpath/data/coglist.txt";    #-- COG equivalence file (COGid -> Function -> Functional class)
my $kegglist="$installpath/data/keggfun2.txt";  #-- KEGG equivalence file (KEGGid -> Function -> Functional class)
my %ranks=('k',1,'p',1,'c',1,'o',1,'f',1,'g',1,'s',1);    #-- Only these taxa will be considered for output

my $resultsdir="$pwd/$project";
print "----$resultsdir---\n";
my @fields=split(/\//, $resultsdir);
my $project=$fields[-1];
if (-d $resultsdir) { print RED "WARNING: Project name $resultsdir already exists\n"; print RESET; print outsyslog "WARNING: Project name $resultsdir already exists\n"; } else { system("mkdir $resultsdir"); }
$methodsfile="$resultsdir/methods.txt";
open(outmet,">$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "Analysis done with SqueezeMeta on Reads v$version (Tamames & Puente-Sanchez 2019, Frontiers in Microbiology 9, 3349)\n";
open(outsyslog,">$resultsdir/syslog") || warn "Cannot open syslog file in $resultsdir/syslog\n";
print outsyslog "Created by $0, ",scalar localtime,"\nCommand: $commandline\n";

if(!$nodiamond) { print outmet "Similarity searches for"; }

my $output_all="$project.out.allreads";
open(outall,">$resultsdir/$output_all") || die;

my $output_counts="$project.out.mappingstat";
open(outcount,">$resultsdir/$output_counts") || die;

#-- Reading the sample file 

print "Now reading samples from $equivfile\n";
print outsyslog "Now reading samples from $equivfile\n";
open(infile1,$equivfile) or do { print RED "Cannot open samples file $equivfile\n"; print RESET; print outsyslog "Cannot open samples file $equivfile\n"; die; };
while(<infile1>) {
	chomp;
	$_=~s/\r//g; # Remove DOS line terminators
	next if(!$_ || ($_=~/^\#/));
	my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
	if($_=~/ /) { print RED "Please do not use blank spaces in the samples file\n"; print RESET; print outsyslog "Please do not use blank spaces in the samples file\n"; die; }
	if(($iden ne "pair1") && ($iden ne "pair2")) { print RED "Samples file, line $_: file label must be \"pair1\" or \"pair2\". For single reads, use \"pair1\"\n";  print RESET; print outsyslog "Samples file, line $_: file label must be \"pair1\" or \"pair2\"\n"; die; }
	if((!$sample) || (!$file) || (!$iden)) { print RED "Bad format in samples file $equivfile. Missing fields\n"; print RESET; print outsyslog "Bad format in samples file $equivfile. Missing fields\n"; die; }
	if(-e "$rawseqs/$file") {} else { print RED "Cannot find sample file $rawseqs/$file for sample $sample in the samples file. Please check\n"; print RESET; print outsyslog "Cannot find sample file $rawseqs/$file for sample $sample in the samples file\n"; die; }
	$allsamples{$sample}{$file}=$iden;
	$ident{$sample}{$file}=$iden;
}
close infile1;

my($extdbname,$extdb,$dblist,$optdbsw);
if($opt_db) {
	open(infile0,$opt_db) or do { print RED "WARNING: Can't open EXTDB file $opt_db\n"; }; 
	while(<infile0>) {
	($extdbname,$extdb,$dblist)=split(/\t/,$_);
	$allext{$extdbname}=$dblist;
			}
	close infile0;
	}

my @nmg=keys %allsamples;
my $numsamples=$#nmg+1;
my $sampnum;
print "$numsamples metagenomes found";
print "\n";
print outsyslog "$numsamples metagenomes found\n";
print outall "# Created by $0 from data in $equivfile, ", scalar localtime,"\n";
print outall "# Sample\tFile\tRead\tTax\tConsensus tax\t";
if(!$nocog) { print outall "\tCOG"; }
if(!$nokegg) { print outall "\tKEGG"; }
if($opt_db) {  foreach my $extdb(sort keys %allext) { print outall "\t$extdb"; } }
print outall "\n";
print outcount "# Created by $0 from data in $equivfile, ", scalar localtime,"\n";
print outcount "# Sample\tFile\tTotal Reads\tReads with hits\tTotal number of hits\n";

	
my(%cogaccum,%keggaccum,%rblast,%iblast,%store,%inputfile);
my($thisfile,$numseqs);
foreach my $thissample(keys %allsamples) {
	%store=();
	$sampnum++;
	print BOLD "\nSAMPLE $sampnum/$numsamples: $thissample\n\n"; print RESET;
	print outsyslog "\nSAMPLE $sampnum/$numsamples: $thissample\n\n"; 
	my $thissampledir="$resultsdir/$thissample";
	if(-d $thissampledir) {} else { system("mkdir $thissampledir"); }
	foreach my $thisfile(sort keys %{ $allsamples{$thissample} }) {                
		(%iblast,%rblast)=();
		my $numseqs=0;
		print "   File: $thisfile\n";
		print outsyslog "   File: $thisfile\n";
		my $idenf=$allsamples{$thissample}{$thisfile};

		#-- Transforming the input to a gunzipped fasta file

		my $fastqfile="$rawseqs/$thisfile";
		my $fastafile="$thissampledir/$thisfile";
		system("cp $fastqfile $fastafile");
		if($fastafile=~/gz$/) { system("gunzip $fastafile"); $fastafile=~s/\.gz$//; }
                my $fastaname=$fastafile;
		$fastaname=~s/fastq$/fasta/;
		if($fastafile=~/fastq$|fq$/) {
			$fastaname=~s/fastq$/fasta/;
			$fastaname=~s/fq$/fasta/;
			open(outfasta,">$fastaname") || die;
			my($nline,$header,$seqline);
			open(infastq,$fastafile) || die;
			while(<infastq>) {
				$header=$_;
				$header=~s/^\@//;
				$header=~s/\s+.*//;
				chomp $header;
				$inputfile{$header}=$thisfile;
				$seqline=<infastq>;
				$_=<infastq>;
				$_=<infastq>;
				print outfasta ">$header\n$seqline";
				}
			close infastq;
			close infasta;
			}
		my $numseqs;
		system("wc -l $fastaname > $thissampledir/rc.txt");
		open(inf,"$thissampledir/rc.txt") || die;
		$numseqs=<inf>;
		close inf;
		$numseqs=~s/\s+.*//;
		$numseqs/=2;
		$totalseqs{$thisfile}=$numseqs;
		system("rm $thissampledir/rc.txt");

	#---- Taxonomic annotation

		$currtime=timediff();
		print CYAN "[",$currtime->pretty,"]: Starting taxonomic annotation\n"; print RESET;		
		print outsyslog "[",$currtime->pretty,"]: Starting taxonomic annotation\n";		
		my $blastxout="$thissampledir/$thissample.nr.blastx";
		my $collapsed="$thissampledir/$thissample.nr.blastx.collapsed.m8";
		my $collapsedmerged=$collapsed;
		$collapsedmerged=~s/\.m8/\.merged\.m8/;
		my $ntseqs="$thissampledir/$thissample.nt.fasta";
		my $wrankfile="$thissampledir/$thissample.fun3.blastx.tax.wranks";
		my $wrankfile_nofilter="$thissampledir/$thissample.fun3.blastx.tax_nofilter.wranks";
		if($nodiamond) { print "   (Skipping Diamond search for taxa because of --nodiamond flag)\n"; } 
		else { 
			print CYAN "[",$currtime->pretty,"]: Running Diamond (Buchfink et al 2015, Nat Methods 12, 59-60) for taxa (GenBank nr, Clark et al 2016, Nucleic Acids Res 44, D67-D72)\n"; print RESET;
			print outsyslog "[",$currtime->pretty,"]: Running Diamond for taxa\n";
			print outmet "GenBank (Clark et al 2016, Nucleic Acids Res 44, D67-D72), ";
			if(-s $blastxout>0) { print "  Diamond result found in $blastxout, not running it again\n"; }
			else { run_blastx($fastafile,$blastxout); }
			if(-s $collapsed>0) { print "  Diamond collapsed result found in $collapsed, not running it again\n"; }
			else { collapse($blastxout,$collapsed); }
			if(-s $collapsedmerged>0) { print "  Diamond collapsed and merged result found in $collapsedmerged, not running it again\n"; }
			else { merge($collapsed,$collapsedmerged); }
			if(-s $ntseqs>0) { print "  ORFs sequences found in $ntseqs, not running it again\n"; }
			else { getseqs($collapsedmerged,$fastaname,$ntseqs); }
			}

		if(-s $wrankfile>0) { print "  Tax annotations found in $wrankfile, not running it again\n"; }
		else { lca($collapsedmerged,$thissampledir,$scriptdir,$thissample,$numthreads); }
		open(outsyslog,">>$resultsdir/syslog");
		my $numtotalhits;
		open(inf,$collapsedmerged);
		while(<inf>) {
			chomp;
			next if !$_;
			my @h=split(/\t/,$_);
			my @u=split(/\_/,$h[0]);
			pop @u;
			my $ctg=join("_",@u);
			$iblast{$ctg}=1;
			$rblast{$h[0]}=1;
			}
		close inf;
		open(infiletax,$wrankfile) || die "Cannot open file $wrankfile\n";
		while(<infiletax>) {
			chomp;
			next if(!$_ || ($_=~/^\#/));
			my @f=split(/\t/,$_);
			$store{$f[0]}{tax}=$f[1];
			}
		close infiletax;
		if($euknofilter) {     #-- Drops the filters for eukaryotes
			open(infiletax,$wrankfile_nofilter) || die "Cannot open file $wrankfile_nofilter\n";
			while(<infiletax>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				if($f[1]=~/Eukaryota/) { $store{$f[0]}{tax}=$f[1]; }
				}
			close infiletax;
			}
		
		#-- Run consensus annotation of reads

		print "  Running consensus annotation: Output in $thissampledir/readconsensus.txt\n";
		my $command="$installpath/lib/SQM_reads/readconsensus.pl $thissample $thissampledir $euknofilter $installpath $databasepath";
		print outsyslog "  Running consensus annotation: Output in $thissampledir/readconsensus.txt: $command\n";
		# print "$command\n";
		system($command);

	#---- Functional annotation
		
		$currtime=timediff();
                print CYAN "[",$currtime->pretty,"]: Starting functional annotation\n"; print RESET; 
                print outsyslog"[",$currtime->pretty,"]: Starting functional annotation\n"; 
		if(!$nocog) {
			print outsyslog "Starting COG annotation\n";
			my $outfile="$thissampledir/$thisfile.cogs.m8";
			my $blastx_command="$diamond_soft blastx -q $ntseqs -p $numthreads -d $cog_db -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
			#print "Running BlastX: $blastx_command\n";
			if($nodiamond) { print "   (Skipping Diamond run for COGs because of --nodiamond flag)\n"; print outsyslog"   (Skipping Diamond run for COGs because of --nodiamond flag)\n"; } 
			else { 
				print CYAN "[",$currtime->pretty,"]: Running Diamond for COGs\n"; print RESET;
				print outsyslog "[",$currtime->pretty,"]: Running Diamond for COGs: $blastx_command\n";
				system($blastx_command); 
				print outmet "eggNOG (Huerta-Cepas et al 2016, Nucleic Acids Res 44, D286-93), ";
			}
			my $outfile_cog="$thissampledir/$thisfile.cogs";
			my $func_command="perl $auxdir/func.pl $outfile $outfile_cog";
			$currtime=timediff();
			print "[",$currtime->pretty,"]: Running functional annotations for COGs\n"; print RESET;
			print outsyslog "[",$currtime->pretty,"]: Running functional annotations for COGs: $func_command\n";
			system($func_command);
			open(infilecog,$outfile_cog) || die;
			while(<infilecog>) { 
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				my $orfid="$f[0]";
				$store{$orfid}{cog}=$f[1];
				if($f[1] eq $f[2]) { $store{$orfid}{cog}.="*"; }
				my $ctg=$orfid;
				$ctg=~s/\_\d+\-\d+$//;
	                        $iblast{$ctg}=1;
	                        $rblast{$orfid}=1;

				}
			close infilecog;
			}
			
		if(!$nokegg) {
			print outsyslog "Starting KEGG annotation\n";
			$currtime=timediff();
			my $outfile="$thissampledir/$thisfile.kegg.m8";
			my $blastx_command="$diamond_soft blastx -q $ntseqs -p $numthreads -d $kegg_db -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
			#print "Running BlastX: $blastx_command\n";
			if($nodiamond) { print "   (Skipping Diamond run for KEGG because of --nodiamond flag)\n"; print outsyslog "   (Skipping Diamond run for KEGG because of --nodiamond flag)\n"; }
			else { 
				print CYAN "[",$currtime->pretty,"]: Running Diamond for KEGG\n"; print RESET;
				print outsyslog "[",$currtime->pretty,"]: Running Diamond for KEGG: $blastx_command\n";
				system($blastx_command); 
				print outmet "KEGG (Kanehisa and Goto 2000, Nucleic Acids Res 28, 27-30), ";
			}
			my $outfile_kegg="$thissampledir/$thisfile.kegg";
			my $func_command="perl $auxdir/func.pl $outfile $outfile_kegg";
			$currtime=timediff();
			print CYAN "[",$currtime->pretty,"]: Running functional annotation for KEGG\n"; print RESET;
			print outsyslog "[",$currtime->pretty,"]: Running functional annotation for KEGG: $func_command\n"; 
			system($func_command);
			open(infilekegg,$outfile_kegg) || die;
			while(<infilekegg>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				my $orfid="$f[0]";
				$store{$orfid}{kegg}=$f[1];
				if($f[1] eq $f[2]) { $store{$orfid}{kegg}.="*"; }
                                my $ctg=$orfid;
                                $ctg=~s/\_\d+\-\d+$//;
                                $iblast{$ctg}=1;
                                $rblast{$orfid}=1;
				}
			close infilekegg;
			}
		if($opt_db) {	
			open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
			while(<infile0>) {
				chomp;
				next if(!$_ || ($_=~/\#/));
				my($extdbname,$extdb,$dblist)=split(/\t/,$_);
				print outsyslog "Starting $extdbname annotation\n";
				$currtime=timediff();
				my $outfile="$thissampledir/$thisfile.$extdbname.m8";
				my $blastx_command="$diamond_soft blastx -q $ntseqs -p $numthreads -d $extdb -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
				#print "Running BlastX: $blastx_command\n";
				if($nodiamond) { print "   (Skipping Diamond run for $extdbname because of --nodiamond flag)\n"; }
				else { 
					print CYAN "[",$currtime->pretty,"]: Running Diamond for $extdbname\n"; print RESET;
					print outsyslog "[",$currtime->pretty,"]: Running Diamond for $extdbname: $blastx_command\n";
					system($blastx_command); 
					print outmet "$extdbname, ";
				}
				my $outfile_opt="$thissampledir/$thisfile.$extdbname";
				my $func_command="perl $auxdir/func.pl $outfile $outfile_opt";
				$currtime=timediff();
				print "[",$currtime->pretty,"]: Running functional annotation for $extdbname\n"; print RESET;
				print outsyslog "[",$currtime->pretty,"]: Running functional annotation for $extdbname: $func_command\n"; 
				system($func_command);
				open(infileopt,$outfile_opt) || die;
				while(<infileopt>) {
					chomp;
					next if(!$_ || ($_=~/^\#/));
					my @f=split(/\t/,$_);
					my $orfid="$f[0]";
					$store{$orfid}{$extdbname}=$f[1];
					if($f[1] eq $f[2]) { $store{$orfid}{$extdbname}.="*"; }
	                                my $ctg=$orfid;
        	                        $ctg=~s/\_\d+\-\d+$//;
                	                $iblast{$ctg}=1;
                        	        $rblast{$orfid}=1;
					}
				close infileopt;
				}
			close infile0;
			}

	#-- Mapping counts

        my @y=keys %iblast;
        my $numhits=($#y)+1;
        my @y=keys %rblast;
        my $numtotalhits=($#y)+1;
	print outsyslog "Counting mapped counts\n";
        print outcount "$thissample\t$thisfile\t$numseqs\t$numhits\t$numtotalhits\n";
#	system("rm $thissampledir/diamond_collapse*; rm $thissampledir/collapsed*m8; rm $thissampledir/rc.txt; rm $thissampledir/wc;");
 	system("rm -r $thissampledir/temp");

		}

	#-- Global statistics

	my %consannot;
	my $consannotation="$thissampledir/readconsensus.txt";
	print outsyslog "Making global statistics: Reading from $consannotation\n";
	
	open(infile5,$consannotation) || die "Cannot open consensus annotation in $consannotation\n";
	while(<infile5>) {
		chomp;
		next if !$_;
		my($readname,$constax)=split(/\t/,$_);
		my @tfields=split(/\;/,$constax);        #-- As this will be a huge file, we do not report the full taxonomy, just the deepest taxon
                my $lastconstax=$tfields[$#tfields];
		$consannot{$readname}=$lastconstax;
		$accum{$thissample}{taxread}{$constax}++;
		}
	close infile5;
		
	if(!$nodiamond) { print outmet " were done using Diamond (Buchfink et al 2015, Nat Methods 12, 59-60)\n"; }		
	my @listorfs;
	foreach my $orf(keys %store) { 
		my @j=split(/\_/,$orf);
		my $pos=pop @j;
		my $contname=join("_",@j);
		my($poinit,$poend)=split("-",$pos); 
		push(@listorfs,{'orf',=>$orf,'contig'=>$contname,'posinit'=>$poinit});
		}
		my @sortedorfs=sort {
		$a->{'contig'} cmp $b->{'contig'} ||
		$a->{'posinit'} <=> $b->{'posinit'}
		} @listorfs;
	#foreach my $k(sort keys %store) {
	foreach my $orf(@sortedorfs) {
		my $k=$orf->{'orf'};
		my @tfields=split(/\;/,$store{$k}{tax});	#-- As this will be a huge file, we do not report the full taxonomy, just the deepest taxon
		my $lasttax=$tfields[$#tfields];
		my $tread=$k;
		$tread=~s/\_\d+\-\d+$//;
		my $ifile=$inputfile{$tread};
		print outall "$thissample\t$ifile\t$k\t$lasttax\t$consannot{$tread}\t";
		if(!$nocog) { print outall "\t$store{$k}{cog}"; }
		if(!$nokegg) { print outall "\t$store{$k}{kegg}"; }
		if($opt_db) { 
			foreach my $extdb(sort keys %allext) { print outall "\t$store{$k}{$extdb}"; }
			}
		print outall "\n";
		$store{$k}{cog}=~s/\*//;
		$store{$k}{kegg}=~s/\*//;
		if($lasttax) { $accum{$thissample}{tax}{$store{$k}{tax}}++; }
		if($store{$k}{cog}) { 
			$accum{$thissample}{cog}{$store{$k}{cog}}++; 
			$cogaccum{$store{$k}{cog}}++;
			}
		if($store{$k}{kegg}) { 
			$accum{$thissample}{kegg}{$store{$k}{kegg}}++;		
			$keggaccum{$store{$k}{kegg}}++;	
			}	
		foreach my $topt(sort keys %allext) {
			if($store{$k}{$topt}) { 		
				$store{$k}{$topt}=~s/\*//;
				$accum{$thissample}{$topt}{$store{$k}{$topt}}++;		
				$optaccum{$topt}{$store{$k}{$topt}}++;	
				}
			}

		}

	}
		
close outall;	
close outcount;


#------------ Global tables --------------#

my(%cog,%kegg,%opt);

	#-- Reading data for KEGGs (names, pathways)

print outsyslog "Making global statistics\n";
open(infile2,$kegglist) or do { print RED "WARNING: Missing KEGG equivalence file\n"; print RESET; print outsyslog "WARNING: Missing KEGG equivalence file\n"; };
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$kegg{$t[0]}{name}=$t[1];
	$kegg{$t[0]}{fun}=$t[2];
	$kegg{$t[0]}{path}=$t[3];
	}
close infile2;

foreach my $idb(keys %allext) {
	open(infile2,$allext{$idb}) or do { print RED "WARNING: Missing $idb equivalence file\n"; print RESET; print outsyslog "WARNING: Missing $idb equivalence file\n"; };
	while(<infile2>) {
		chomp;
		$_=~s/\r//g; # Remove windows line terminators
		next if(!$_ || ($_=~/\#/));
		my @t=split(/\t/,$_);
		$opt{$idb}{$t[0]}{fun}=$t[1];
		#if($t[2]) { $opt{$t[0]}{fun}=$t[2]; }
		#if($t[3]) { $opt{$t[0]}{path}=$t[3]; }
		}
	close infile2;
	}


$currtime=timediff();
print CYAN "\n[",$currtime->pretty,"]: Creating global tables\n"; print RESET;
print outsyslog "\n[",$currtime->pretty,"]: Creating global tables\n"; 
print "   Tax table: $resultsdir/$output_all.mcount\n";
print outsyslog "   Tax table: $resultsdir/$output_all.mcount\n";		
open(outtax,">$resultsdir/$output_all.mcount");
print outtax "# Created by $0 from data in $equivfile", scalar localtime,"\n";
print outtax "Rank\tTax\tTotal ORFs\tTotal reads";
foreach my $sprint(sort keys %accum) { print outtax "\t$sprint ORFs\t$sprint reads"; }
print outtax "\n";
my(%taxaccum,%taxaccumread);
foreach my $isam(sort keys %accum) {
	foreach my $itax(keys %{ $accum{$isam}{tax} }) {
		$itax=~s/\;$//;
		my @stx=split(/\;/,$itax);
		my $thisrank;
		foreach my $tf(@stx) {
			$thisrank.="$tf;";
			$taxaccum{$isam}{$thisrank}+=$accum{$isam}{tax}{$itax};
			$taxaccum{total}{$thisrank}+=$accum{$isam}{tax}{$itax};
			}
		}
        foreach my $itax(keys %{ $accum{$isam}{taxread} }) {
                $itax=~s/\;$//;
                my @stx=split(/\;/,$itax);
                my $thisrank;
                foreach my $tf(@stx) {
                        $thisrank.="$tf;";
                        $taxaccumread{$isam}{$thisrank}+=$accum{$isam}{taxread}{$itax};
                        $taxaccumread{total}{$thisrank}+=$accum{$isam}{taxread}{$itax};
                        }
                }

	}
foreach my $ntax(sort { $taxaccum{total}{$b}<=>$taxaccum{total}{$a}; } keys %{ $taxaccum{total} }) {
	my @stx=split(/\;/,$ntax);
	my($lastrank,$lasttax)=split(/\_/,$stx[$#stx]);
	next if(!$ranks{$lastrank});
	print outtax "$lastrank\t$ntax\t$taxaccum{total}{$ntax}\t$taxaccumread{total}{$ntax}";
	foreach my $isam(sort keys %accum) {
		my $dato=$taxaccum{$isam}{$ntax} || "0";
		my $datoread=$taxaccumread{$isam}{$ntax} || "0";
		print outtax "\t$dato\t$datoread";
		}
	print outtax "\n";
	}

close outtax;	 
	
if(!$nocog) {
	open(infile1,$coglist) || warn "Missing COG equivalence file\n";
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my @t=split(/\t/,$_);
		$cog{$t[0]}{fun}=$t[1];
		$cog{$t[0]}{path}=$t[2]; 
           	 }
	close infile1;

	print "   COG table: $resultsdir/$output_all.funcog\n";		
	print outsyslog "   COG table: $resultsdir/$output_all.funcog\n";		
	open(outcog,">$resultsdir/$output_all.funcog");
	print outcog "# Created by $0 from data in $equivfile", scalar localtime,"\n";
	print outcog "COG\tTotal";
	foreach my $sprint(sort keys %accum) { print outcog "\t$sprint"; }
	print outcog "\tFunction\tClass\n";
	foreach my $ncog(sort { $cogaccum{$b}<=>$cogaccum{$a}; } keys %cogaccum) {
		print outcog "$ncog\t$cogaccum{$ncog}";
		foreach my $isam(sort keys %accum) {
			my $dato=$accum{$isam}{cog}{$ncog} || "0";
			print outcog "\t$dato";
			}
		print outcog "\t$cog{$ncog}{fun}\t$cog{$ncog}{path}\n";
		}

	close outcog;
	}	 

if(!$nokegg) {
	print "   KEGG table: $resultsdir/$output_all.funkegg\n";		
	print outsyslog "   KEGG table: $resultsdir/$output_all.funkegg\n";		
	open(outkegg,">$resultsdir/$output_all.funkegg");
	print outkegg "# Created by $0 from data in $equivfile", scalar localtime,"\n";
	print outkegg "KEGG\tTotal";
	foreach my $sprint(sort keys %accum) { print outkegg "\t$sprint"; }
	print outkegg "\tFunction\tClass\n";
	foreach my $nkegg(sort { $keggaccum{$b}<=>$keggaccum{$a}; } keys %keggaccum) {
		print outkegg "$nkegg\t$keggaccum{$nkegg}";
		foreach my $isam(sort keys %accum) {
			my $dato=$accum{$isam}{kegg}{$nkegg} || "0";
			print outkegg "\t$dato";
			}
		print outkegg "\t$kegg{$nkegg}{fun}\t$kegg{$nkegg}{path}\n";
		}

	close outkegg;
	}	 

if($opt_db) {
	foreach my $extdbname(keys %allext) {
		print "   $extdbname table: $resultsdir/$output_all.fun$extdbname\n";		
		print outsyslog "   $extdbname table: $resultsdir/$output_all.fun$extdbname\n";		
		open(outopt,">$resultsdir/$output_all.fun$extdbname");
		print outopt "# Created by $0 from data in $opt_db", scalar localtime,"\n";
		print outopt "$extdbname\tTotal";
		foreach my $sprint(sort keys %accum) { print outopt "\t$sprint"; }
		print outopt "\tFunction\n";
		foreach my $nopt(sort { $optaccum{$extdbname}{$b}<=>$optaccum{$extdbname}{$a}; } keys %{ $optaccum{$extdbname} }) {
			print outopt "$nopt\t$optaccum{$extdbname}{$nopt}";
			foreach my $isam(sort keys %accum) {
				my $dato=$accum{$isam}{$extdbname}{$nopt} || "0";
				print outopt "\t$dato";
				}
			 print outopt "\t$opt{$extdbname}{$nopt}{fun}\n";
			}
		close outopt;
		}
	}	 

print "   Mapping statistics: $resultsdir/$output_counts\n";
print "   Condensed annotations for mapped reads: $resultsdir/$output_all\n";
print outsyslog "   Mapping statistics: $resultsdir/$output_counts\n";
print outsyslog "   Condensed annotations for mapped reads: $resultsdir/$output_all\n";

$currtime=timediff();
print CYAN "\n[",$currtime->pretty,"]: DONE! Have fun!\n"; print RESET;
print outsyslog "\n[",$currtime->pretty,"]: DONE! Have fun!\n";
close outmet;
print "For citation purposes, you can find a summary of methods in the file $methodsfile\n";
print RESET;
print outsyslog "\nNormal termination, all good apparently\n";
close outsyslog;


sub run_blastx {

	#-- Run Diamond search
	my $queryfile=shift;
	my $blastxout=shift;
	print "  Running Diamond BlastX (Buchfink et al 2015, Nat Methods 12, 59-60)\n";
	my $blastx_command="$diamond_soft blastx -q $queryfile -p $numthreads -d $nr_db -f tab -F 15 --quiet --range-culling -b $blocksize -e $evalue --id $miniden --top 10 -o $blastxout -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen";	
	print outsyslog "Running Diamond BlastX: $blastx_command\n";
	print outmet "Additional ORFs were obtained by Diamond BlastX (Buchfink et al 2015, Nat Methods 12, 59-60)\n";
	print "Running Diamond Blastx: $blastx_command**\n" if $verbose;
	system $blastx_command;
	}

sub collapse {

	#-- Collapse hits using blastxcollapse.pl
	my($blastxout,$collapsed)=@_;
	print "  Collapsing hits with blastxcollapse.pl\n";
	my $collapse_command="$installpath/lib/SqueezeMeta/blastxcollapse.pl -n -s -f -m 50 -l 70 -p $numthreads $blastxout > $collapsed";
	print outsyslog "Collapsing hits with blastxcollapse.pl: $collapse_command\n";
	print "Collapsing hits with blastxcollapse.pl: $collapse_command\n" if $verbose;
	close outsyslog;
	system $collapse_command;
	open(outsyslog,">>$resultsdir/syslog");
	}
	
sub merge {

	#-- Merge frameshifts
	my($collapsed,$collapsedmerged)=@_;
	my $merge_command="$installpath/lib/SqueezeMeta/mergehits.pl $collapsed > $collapsedmerged";
	print "  Merging splitted hits with mergehits.pl\n";
	print outsyslog "Merging splitted hits with mergehits.pl: $merge_command\n";
	print "Merging splitted hits with mergehits.pl: $merge_command\n" if $verbose;
	close outsyslog;
	system $merge_command;
	open(outsyslog,">>$resultsdir/syslog");
	}

sub getseqs {

	#-- Get new nt sequences
	my($collapsedmerged,$fastafile,$ntseqs)=@_;
	print "  Getting nt sequences from $fastafile\n";
	my %orfstoget;
	open(infile4,$collapsedmerged) || die "Can't open $collapsedmerged\n";
	while(<infile4>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @t=split(/\t/,$_);
		my @w=split(/\_/,$t[0]);
		my $posn=pop @w;
		my $contname=join("_",@w);
		$orfstoget{$contname}{$posn}=1;
		# print "$contname*$posn*\n";
		}
	close infile4;
	my($currcontig,$newcontig,$contigseq);
	open(outfile2,">$ntseqs") || die "Can't open $ntseqs for writing\n";
	open(infile5,$fastafile) || die "Can't open $fastafile\n";
	while(<infile5>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if($_=~/^\>([^ ]+)/) {
			$newcontig=$1; 
			if($currcontig) { 
				foreach my $gorf(keys %{ $orfstoget{$currcontig} }) { 
					my($pinit,$pend)=split(/\-/,$gorf);
					my $tlen=$pend-$pinit+1;
					my $mseq=substr($contigseq,$pinit,$tlen);
					print outfile2 ">$currcontig\_$gorf\n$mseq\n";
					}
				}
			$currcontig=$newcontig;
			$contigseq="";
			}
		else { $contigseq.=$_; }
	}	
	close infile5;
	foreach my $gorf(keys %{ $orfstoget{$currcontig} }) {
		my($pinit,$pend)=split(/\-/,$gorf);
		my $tlen=length($pend-$pinit+1);
		my $mseq=substr($contigseq,$pinit,$tlen);
		print outfile2 ">$currcontig\_$gorf\n$mseq\n";
		}
	close outfile2;	
	print "  Sequences stored in $ntseqs\n";				
	}


sub lca {
	my($collapsedmerged,$thissampledir,$scriptdir,$thissample,$numthreads)=@_;
	my $lca_command="perl $auxdir/lca_collapse.pl $collapsedmerged $thissampledir $scriptdir $thissample $numthreads";
	$currtime=timediff();
	print CYAN "[",$currtime->pretty,"]: Running LCA\n"; print RESET;
	close(outsyslog);
	system($lca_command);
	open(outsyslog,">>$resultsdir/syslog");
	}

#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}





