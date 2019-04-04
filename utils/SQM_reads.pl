#!/usr/bin/perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use lib ".";
use strict;
use File::Basename;
our $utilsdir = dirname(__FILE__);
our $installpath = "$utilsdir/..";
our $scriptdir = "$installpath/scripts";
our $auxdir = "$installpath/lib/SQM_reads";

my $version="0.1.0, Feb 2019";
my $start_run = time();

do "$scriptdir/SqueezeMeta_conf.pl";
#-- Configuration variables from conf file
our($databasepath);


my($numthreads,$project,$equivfile,$rawseqs,$evalue,$dietext,$blocksize,$currtime,$nocog,$nokegg);

my $helptext = <<END_MESSAGE;
Usage: SQM_reads.pl -p <project name> -s <samples file> -f <raw fastq dir> <options>

Arguments:

 Mandatory parameters:
   -p: Project name (REQUIRED)
   -s|-samples: Samples file (REQUIRED)
   -f|-seq: Fastq read files' directory (REQUIRED)
   
 Options:
   --nocog: Skip COG assignment (Default: no)
   --nokegg: Skip KEGG assignment (Default: no)
   -e|-evalue: max evalue for discarding hits diamond run  (Default: 1e-03)
   -t: Number of threads (Default: 12)
   -b|-block-size: block size for diamond against the nr database (Default: 8)


END_MESSAGE

my $result = GetOptions ("t=i" => \$numthreads,
                     "p=s" => \$project,
                     "s|samples=s" => \$equivfile,
                     "f|seq=s" => \$rawseqs, 
		     "e|evalue=f" => \$evalue,   
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
                     "b|block_size=i" => \$blocksize,
		    );

if(!$numthreads) { $numthreads=12; }
if(!$evalue) { $evalue=0.001; }

print "\nSqueezeMeta on Reads v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n";

if(!$project) { $dietext.="MISSING ARGUMENT: -p: Project name\n"; }
if(!$rawseqs) { $dietext.="MISSING ARGUMENT: -f|-seq:Read files' directory\n"; }
if(!$equivfile) { $dietext.="MISSING ARGUMENT: -s|-samples: Samples file\n"; }
if($dietext) { die "$dietext\n$helptext\n"; }

my(%allsamples,%ident,%noassembly,%accum,%totalseqs);
my($sample,$file,$iden,$mapreq);
tie %allsamples,"Tie::IxHash";

my $nr_db="$databasepath/nr.dmnd";
my $cog_db="$databasepath/eggnog";
my $kegg_db="$databasepath/keggdb";
my $diamond_soft="$installpath/bin/diamond";
my $coglist="$installpath/data/coglist.txt";    #-- COG equivalence file (COGid -> Function -> Functional class)
my $kegglist="$installpath/data/keggfun2.txt";  #-- KEGG equivalence file (KEGGid -> Function -> Functional class)
my %ranks=('k',1,'phylum',1,'c',1,'o',1,'f',1,'g',1,'s',1);    #-- Only these taxa will be considered for output

my $resultsdir=$project;
if (-d $resultsdir) { die "Project name $project already exists\n"; } else { system("mkdir $resultsdir"); }

my $output_all="$project.out.allreads";
open(outall,">$resultsdir/$output_all") || die;

my $output_counts="$project.out.counts";
open(outcount,">$resultsdir/$output_counts") || die;

#-- Reading the sample file 

print "Now reading samples from $equivfile\n";
open(infile1,$equivfile) || die "Cannot open samples file $equivfile\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
	if($_=~/ /) { die "Please do not use blank spaces in the samples file\n"; }
	if(($iden ne "pair1") && ($iden ne "pair2")) { die "Samples file, line $_: file label must be \"pair1\" or \"pair2\". For single reads, use \"pair1\"\n"; }
	if((!$sample) || (!$file) || (!$iden)) { die "Bad format in samples file $equivfile. Missing fields\n"; }
	if(-e "$rawseqs/$file") {} else { die "Cannot find sample file $rawseqs/$file for sample $sample in the samples file. Please check\n"; }
	$allsamples{$sample}{$file}=1;
	$ident{$sample}{$file}=$iden;
}
close infile1;

my @nmg=keys %allsamples;
my $numsamples=$#nmg+1;
my $sampnum;
print "$numsamples metagenomes found";
print "\n";
print outall "# Created by $0 from data in $equivfile", scalar localtime,"\n";
print outall "# Sample\tRead\tTax\tCOG\tKEGG\n";
print outcount "# Created by $0 from data in $equivfile", scalar localtime,"\n";
print outcount "# Sample\tFile\tReads\tReads with hits\n";

my(%cogaccum,%keggaccum);
foreach my $thissample(keys %allsamples) {
	my %store;
	$sampnum++;
	print "\nSAMPLE $sampnum/$numsamples: $thissample\n\n"; 
	my $thissampledir="$resultsdir/$thissample";
	system("mkdir $thissampledir");
	foreach my $thisfile(sort keys %{ $allsamples{$thissample} }) {                
		print "   File: $thisfile\n";
		if($thisfile=~/fastq.gz/) { system("zcat $thisfile | wc > rc.txt"); }
		elsif($thisfile=~/fastq/) { system("wc $thisfile > rc.txt"); }
		elsif($thisfile=~/fasta.gz/) { system("zcat $thisfile | grep -c \"^>\" > rc.txt"); }
		elsif($thisfile=~/fasta/) { system("grep -c \"^>\" $thisfile > rc.txt"); }
		open(inw,"rc.txt");
		my $numseqs=<inw>;
		close inw;
		chomp $numseqs;
		if($thisfile=~/fastq/) { $numseqs/=4; }
		system("rm rc.txt");
		$totalseqs{$thisfile}++;
		$currtime=timediff();
		print "   [",$currtime->pretty,"]: Running Diamond for taxa\n";
		my $outfile="$thissampledir/$thisfile.tax.m8";
		my $outfile_tax="$thissampledir/$thisfile.tax.wranks";
		my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $nr_db -e $evalue --quiet -f tab -b 8 -o $outfile";
		# print "Running BlastX: $blastx_command\n";
		my %iblast;
		system($blastx_command);
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
			
		my $lca_command="perl $auxdir/lca_reads.pl $outfile";
		$currtime=timediff();
		print "   [",$currtime->pretty,"]: Running LCA\n";
		system($lca_command);
		open(infiletax,$outfile_tax) || die;
		while(<infiletax>) {
			chomp;
			next if(!$_ || ($_=~/^\#/));
			my @f=split(/\t/,$_);
			$store{$f[0]}{tax}=$f[1];
			}
		close infiletax;
		
		
		$currtime=timediff();
		if(!$nocog) {
			print "   [",$currtime->pretty,"]: Running Diamond for COGs\n";
			my $outfile="$thissampledir/$thisfile.cogs.m8";
			my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $cog_db -e $evalue --id 30 --quiet -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
			#print "Running BlastX: $blastx_command\n";
			system($blastx_command);
			my $outfile_cog="$thissampledir/$thisfile.cogs";
			my $func_command="perl $auxdir/func.pl $outfile $outfile_cog";
			$currtime=timediff();
			print "   [",$currtime->pretty,"]: Running fun3\n";
			system($func_command);
			open(infilecog,$outfile_cog) || die;
			while(<infilecog>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				$store{$f[0]}{cog}=$f[1];
				if($f[1] eq $f[2]) { $store{$f[0]}{cog}.="*"; }
				}
			close infilecog;
			}
			
		if(!$nokegg) {
			$currtime=timediff();
			print "   [",$currtime->pretty,"]: Running Diamond for KEGG\n";
			my $outfile="$thissampledir/$thisfile.kegg.m8";
			my $blastx_command="$diamond_soft blastx -q $rawseqs/$thisfile -p $numthreads -d $kegg_db -e $evalue --id 30 --quiet -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o $outfile";
			#print "Running BlastX: $blastx_command\n";
			system($blastx_command);
			my $outfile_kegg="$thissampledir/$thisfile.kegg";
			my $func_command="perl $auxdir/func.pl $outfile $outfile_kegg";
			$currtime=timediff();
			print "   [",$currtime->pretty,"]: Running fun3\n";
			system($func_command);
			open(infilekegg,$outfile_kegg) || die;
			while(<infilekegg>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @f=split(/\t/,$_);
				$store{$f[0]}{kegg}=$f[1];
				if($f[1] eq $f[2]) { $store{$f[0]}{kegg}.="*"; }
				}
			close infilekegg;
			}
		}
		
		
	foreach my $k(sort keys %store) {
		my @tfields=split(/\;/,$store{$k}{tax});	#-- As this is a huge file, we do not report the full taxonomy, just the deepest taxon
		my $lasttax=$tfields[$#tfields];
		print outall "$thissample\t$k\t$lasttax\t$store{$k}{cog}\t$store{$k}{kegg}\n";
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
		}
	}
		
print "Output in $output_all\n";
close outall;	
close outcount;


#------------ Global tables --------------#

my(%cog,%kegg);

	#-- Reading data for KEGGs (names, pathways)

open(infile2,$kegglist) || warn "Missing KEGG equivalence file\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$kegg{$t[0]}{name}=$t[1];
	$kegg{$t[0]}{fun}=$t[2];
	$kegg{$t[0]}{path}=$t[3];
	}
close infile2;


$currtime=timediff();
print "\n[",$currtime->pretty,"]: Creating global tables\n";
print "Tax table: $resultsdir/$output_all.mcount\n";		
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

	print "COG table: $resultsdir/$output_all.funcog\n";		
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
	print "KEGG table: $resultsdir/$output_all.funkegg\n";		
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

$currtime=timediff();
print "\n[",$currtime->pretty,"]: DONE! Have fun!\n";

#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}

