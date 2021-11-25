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
#-- Configuration variables from conf file
our($databasepath,$installpath);

my($numthreads,$project,$equivfile,$rawseqs,$evalue,$dietext,$blocksize,$currtime,$nocog,$nokegg,$opt_db,$hel,$nodiamond,$euknofilter,$methodsfile,$printversion);

my $helpshort="Usage: SQM_reads.pl -p <project name> -s <samples file> -f <raw fastq dir> [options]\n";

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
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nodiamond" => \$nodiamond,   
		     "extdb=s" => \$opt_db, 
		     "euk" => \$euknofilter,
                     "b|block_size=i" => \$blocksize,
		     "v|version" => \$printversion,
		     "h" => \$hel
		    );

if(!$numthreads) { $numthreads=12; }
if(!$evalue) { $evalue=0.001; }
my $miniden=30;         #-- Minimum identity for the hit
my $querycover=0;	#-- Minimum coverage of hit in query
	
print BOLD "\nSqueezeMeta on Reads v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;
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

my $nr_db="$databasepath/nr.dmnd";
my $cog_db="$databasepath/eggnog";
my $kegg_db="$databasepath/keggdb";
my $diamond_soft="$installpath/bin/diamond";
my $coglist="$installpath/data/coglist.txt";    #-- COG equivalence file (COGid -> Function -> Functional class)
my $kegglist="$installpath/data/keggfun2.txt";  #-- KEGG equivalence file (KEGGid -> Function -> Functional class)
my %ranks=('k',1,'p',1,'c',1,'o',1,'f',1,'g',1,'s',1);    #-- Only these taxa will be considered for output

my $resultsdir=$project;
my @fields=split(/\//, $resultsdir);
my $project=$fields[-1];
if (-d $resultsdir) { print RED "WARNING: Project name $resultsdir already exists\n"; print RESET; } else { system("mkdir $resultsdir"); }
$methodsfile="$resultsdir/methods.txt";
open(outmet,">$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
print outmet "Analysis done with SqueezeMeta on Reads v$version (Tamames & Puente-Sanchez 2019, Frontiers in Microbiology 9, 3349)\n";
if(!$nodiamond) { print outmet "Similarity searches for"; }

my $output_all="$project.out.allreads";
open(outall,">$resultsdir/$output_all") || die;

my $output_counts="$project.out.mappingstat";
open(outcount,">$resultsdir/$output_counts") || die;

#-- Reading the sample file 

print "Now reading samples from $equivfile\n";
open(infile1,$equivfile) or do { print RED "Cannot open samples file $equivfile\n"; print RESET; die; };
while(<infile1>) {
	chomp;
	$_=~s/\r//g; # Remove windows line terminators
	next if(!$_ || ($_=~/^\#/));
	my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
	if($_=~/ /) { print RED "Please do not use blank spaces in the samples file\n"; print RESET; die; }
	if(($iden ne "pair1") && ($iden ne "pair2")) { print RED "Samples file, line $_: file label must be \"pair1\" or \"pair2\". For single reads, use \"pair1\"\n";  print RESET; die; }
	if((!$sample) || (!$file) || (!$iden)) { print RED "Bad format in samples file $equivfile. Missing fields\n"; print RESET; die; }
	if(-e "$rawseqs/$file") {} else { print RED "Cannot find sample file $rawseqs/$file for sample $sample in the samples file. Please check\n"; print RESET; die; }
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
print outall "# Created by $0 from data in $equivfile, ", scalar localtime,"\n";
print outall "# Sample\tFile\tRead\tTax";
if(!$nocog) { print outall "\tCOG"; }
if(!$nokegg) { print outall "\tKEGG"; }
if($opt_db) {  foreach my $extdb(sort keys %allext) { print outall "\t$extdb"; } }
print outall "\n";
print outcount "# Created by $0 from data in $equivfile, ", scalar localtime,"\n";
print outcount "# Sample\tFile\tTotal Reads\tReads with hits to nr\n";

	
my(%cogaccum,%keggaccum);
foreach my $thissample(keys %allsamples) {
	my %store;
	$sampnum++;
	print BOLD "\nSAMPLE $sampnum/$numsamples: $thissample\n\n"; print RESET;
	my $thissampledir="$resultsdir/$thissample";
	if(-d $thissampledir) {} else { system("mkdir $thissampledir"); }
	foreach my $thisfile(sort keys %{ $allsamples{$thissample} }) {                
		print "   File: $thisfile\n";
		my $idenf=$allsamples{$thissample}{$thisfile};
		my($seqformat,$seqcompress);
		if($thisfile=~/\.fq\.|\.fq$|\.fastq\.|\.fastq$/) { $seqformat="fastq"; }
		elsif($thisfile=~/\.fa\.|\.fa$|\.fasta\.|\.fasta$/) { $seqformat="fasta"; }
		else { die "Unrecognized file format (must be fq, fastq, fa or fasta)\n"; }
		if($thisfile=~/\.gz$/) { $seqcompress="gz"; }
		my $command;
		if($seqformat eq "fastq") {
			if($seqcompress) { $command="zcat $rawseqs/$thisfile | wc > rc.txt"; }
			else { $command="wc $rawseqs/$thisfile > rc.txt"; }
			}
		else {
			if($seqcompress) { $command="zcat $rawseqs/$thisfile | grep -c \"^>\" > rc.txt"; }
			else { $command="grep -c \"^>\" $rawseqs/$thisfile > rc.txt"; }
			}
		system($command);
		# print "$seqformat; Counting reads: $command\n";
		open(inw,"rc.txt");
		my $line=<inw>;
		$line=~s/^\s+//g;
		my @l=split(/\s+/,$line);
		my $numseqs=$l[0];
		close inw;
		chomp $numseqs;
		if($seqformat eq "fastq") { $numseqs/=4; }
		system("rm rc.txt");
		$totalseqs{$thisfile}++;
		$currtime=timediff();
		my $outfile="$thissampledir/$thisfile.tax.m8";
		my $outfile_tax="$thissampledir/$thisfile.tax.wranks";
		my $outfile_tax_nofilter="$thissampledir/$thisfile.tax_nofilter.wranks";
		my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $nr_db -e $evalue --quiet -f tab -b $blocksize -o $outfile";
		# print "Running BlastX: $blastx_command\n";
		my %iblast;
		if($nodiamond) { print "   (Skipping Diamond run because of --nodiamond flag)\n"; } 
		else { 
			if(-e $outfile) { print "Diamond result found in $outfile, skipping run\n"; }
			else {
				print CYAN "[",$currtime->pretty,"]: Running Diamond (Buchfink et al 2015, Nat Methods 12, 59-60) for taxa (GenBank nr, Clark et al 2016, Nucleic Acids Res 44, D67-D72)\n"; print RESET;
				system($blastx_command);
				print outmet "GenBank (Clark et al 2016, Nucleic Acids Res 44, D67-D72), ";
				}
			}
		open(inf,$outfile);
		while(<inf>) {
			chomp;
			next if !$_;
			my @h=split(/\t/,$_);
			$iblast{$h[0]}=1;
			}
		close inf;
		my @y=keys %iblast;
		my $numhits=($#y)+1;
		print outcount "$thissample\t$thisfile\t$numseqs\t$numhits\n";
			
		my $lca_command="perl $auxdir/lca_reads.pl $outfile $numthreads";
		if(-e $outfile_tax) { print "LCA result found in $outfile, skipping run\n"; }
		else {
			$currtime=timediff();
			print CYAN "[",$currtime->pretty,"]: Running LCA\n"; print RESET;
			system($lca_command);
			}
		open(infiletax,$outfile_tax) || die;
		while(<infiletax>) {
			chomp;
			next if(!$_ || ($_=~/^\#/));
			my @f=split(/\t/,$_);
			my $orfid="$f[0]\_$idenf";
			$store{$orfid}{tax}=$f[1];
			}
		close infiletax;
		if($euknofilter) {     #-- Drops the filters for eukaryotes
			open(infiletax,$outfile_tax_nofilter) || die;
			while(<infiletax>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				my $orfid="$f[0]\_$idenf";
				if($f[1]=~/Eukaryota/) { $store{$orfid}{tax}=$f[1]; }
				}
			close infiletax;
			}
		
		$currtime=timediff();
		if(!$nocog) {
			my $outfile="$thissampledir/$thisfile.cogs.m8";
			if(-e $outfile) { print "Diamond result for COGs found in $outfile, skipping run\n"; }
			else {
				my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $cog_db -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
				#print "Running BlastX: $blastx_command\n";
				if($nodiamond) { print "   (Skipping Diamond run because of --nodiamond flag)\n"; } 
				else { 
					print CYAN "[",$currtime->pretty,"]: Running Diamond for COGs\n"; print RESET;
					system($blastx_command); 
					print outmet "eggNOG (Huerta-Cepas et al 2016, Nucleic Acids Res 44, D286-93), ";
					}
				}
			my $outfile_cog="$thissampledir/$thisfile.cogs";
			my $func_command="perl $auxdir/func.pl $outfile $outfile_cog";
			if(-e $outfile_cog) { print "COGs results found in $outfile_cog, skipping run\n"; }
			else {
				$currtime=timediff();
				print CYAN "[",$currtime->pretty,"]: Running fun3\n"; print RESET;
				system($func_command);
				}
			open(infilecog,$outfile_cog) || die;
			while(<infilecog>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				my $orfid="$f[0]\_$idenf";
				$store{$orfid}{cog}=$f[1];
				if($f[1] eq $f[2]) { $store{$orfid}{cog}.="*"; }
				}
			close infilecog;
			}
			
		if(!$nokegg) {
			$currtime=timediff();
			my $outfile="$thissampledir/$thisfile.kegg.m8";
			if(-e $outfile) { print "Diamond result for KEGG found in $outfile, skipping run\n"; }
			else {
				my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $kegg_db -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
				#print "Running BlastX: $blastx_command\n";
				if($nodiamond) { print "   (Skipping Diamond run because of --nodiamond flag)\n"; }
				else { 
				print CYAN "[",$currtime->pretty,"]: Running Diamond for KEGG\n"; print RESET;
					system($blastx_command); 
					print outmet "KEGG (Kanehisa and Goto 2000, Nucleic Acids Res 28, 27-30), ";
					}
				}
			my $outfile_kegg="$thissampledir/$thisfile.kegg";
			if(-e $outfile_kegg) { print "KEGG results found in $outfile_kegg, skipping run\n"; }
			else {
				my $func_command="perl $auxdir/func.pl $outfile $outfile_kegg";
				$currtime=timediff();
				print CYAN "[",$currtime->pretty,"]: Running fun3\n"; print RESET;
				system($func_command);
				}
			open(infilekegg,$outfile_kegg) || die;
			while(<infilekegg>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				my $orfid="$f[0]\_$idenf";
				$store{$orfid}{kegg}=$f[1];
				if($f[1] eq $f[2]) { $store{$orfid}{kegg}.="*"; }
				}
			close infilekegg;
			}
		if($opt_db) {	
			open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
			while(<infile0>) {
				chomp;
				next if(!$_ || ($_=~/\#/));
				my($extdbname,$extdb,$dblist)=split(/\t/,$_);
				$currtime=timediff();
				my $outfile="$thissampledir/$thisfile.$extdbname.m8";
				if(-e $outfile) { print "Diamond result for $extdbname found in $outfile, skipping run\n"; }
				else {
					my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $extdb -e $evalue --query-cover $querycover --id $miniden --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
					#print "Running BlastX: $blastx_command\n";
					if($nodiamond) { print "   (Skipping Diamond run because of --nodiamond flag)\n"; }
					else { 
						print CYAN "[",$currtime->pretty,"]: Running Diamond for $extdbname\n"; print RESET;
						system($blastx_command); 
						print outmet "$extdbname, ";
						}
					}
				my $outfile_opt="$thissampledir/$thisfile.$extdbname";
				if(-e $outfile_opt) { print "$extdbname results found in $outfile_opt, skipping run\n"; }
				else {
					my $func_command="perl $auxdir/func.pl $outfile $outfile_opt";
					$currtime=timediff();
					print CYAN "[",$currtime->pretty,"]: Running fun3\n"; print RESET;
					system($func_command);
					}
				open(infileopt,$outfile_opt) || die;
				while(<infileopt>) {
					chomp;
					next if(!$_ || ($_=~/^\#/));
					my @f=split(/\t/,$_);
					my $orfid="$f[0]\_$idenf";
					$store{$orfid}{$extdbname}=$f[1];
					if($f[1] eq $f[2]) { $store{$orfid}{$extdbname}.="*"; }
					}
				close infileopt;
				}
			}
		}
		
	if(!$nodiamond) { print outmet " were done using Diamond (Buchfink et al 2015, Nat Methods 12, 59-60)\n"; }		
	foreach my $k(sort keys %store) {
		my @tfields=split(/\;/,$store{$k}{tax});	#-- As this will be a huge file, we do not report the full taxonomy, just the deepest taxon
		my $lasttax=$tfields[$#tfields];
		my @j=split(/\_/,$k);
		my $rpair=pop @j;
		my $rread=join("_",@j);
		print outall "$thissample\t$rread\t$rpair\t$lasttax\t";
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

open(infile2,$kegglist) or do { print RED "WARNING: Missing KEGG equivalence file\n"; print RESET; };
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
	open(infile2,$allext{$idb}) or do { print RED "WARNING: Missing $idb equivalence file\n"; print RESET; };
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
print CYAN "\n[",$currtime->pretty,"]: Creating global tables\n";
print "   Tax table: $resultsdir/$output_all.mcount\n";		
open(outtax,">$resultsdir/$output_all.mcount");
print outtax "# Created by $0 from data in $equivfile", scalar localtime,"\n";
print outtax "Rank\tTax\tTotal";
foreach my $sprint(sort keys %accum) { print outtax "\t$sprint"; }
print outtax "\n";
my %taxaccum;
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
	}
foreach my $ntax(sort { $taxaccum{total}{$b}<=>$taxaccum{total}{$a}; } keys %{ $taxaccum{total} }) {
	my @stx=split(/\;/,$ntax);
	my($lastrank,$lasttax)=split(/\_/,$stx[$#stx]);
	next if(!$ranks{$lastrank});
	print outtax "$lastrank\t$ntax\t$taxaccum{total}{$ntax}";
	foreach my $isam(sort keys %accum) {
		my $dato=$taxaccum{$isam}{$ntax} || "0";
		print outtax "\t$dato";
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

$currtime=timediff();
print CYAN "\n[",$currtime->pretty,"]: DONE! Have fun!\n";
close outmet;
print "For citation purposes, you can find a summary of methods in the file $methodsfile\n";
print RESET;

#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}


