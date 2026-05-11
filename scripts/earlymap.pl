#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/01/2019 for version 0.4.3, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage/RPKM for genes/contigs by mapping back reads to the contigs and count how many fall in each gene/contig
#-- Uses bowtie2 for mapping, and sqmapper for counting. 

$|=1;

use File::Basename;
use Cwd 'abs_path';
our $scriptdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $scriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $scriptdir = abs_path(dirname(__FILE__));
        }
our $installpath = abs_path("$scriptdir/..");

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";
use threads;
use Scalar::Util qw(looks_like_number);

my $pwd=cwd();

my $projectdir=$ARGV[0];
my $force_overwrite=$ARGV[1];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$userdir,$bowtieref,$bowtie2_build_soft,$project,$samtools_soft,$contigsfna,$mappingfile,$mapcountfile,$mode,$resultpath,$contigcov,$bowtie2_x_soft, $mappingstat,$nobins,
    $mapper, $mapping_options, $cleaning, $bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$mincontiglen,$doublepass,$contigslen,$gff_file_blastx,$methodsfile,$syslogfile,$keepsam10);

my $verbose=0;

my $fastqdir="$datapath/raw_fastq";
my $bamdir="$datapath/bam";

my $outfile=$mapcountfile;

my $warnmes;

if(-d $bamdir) {} else { system("mkdir $bamdir"); }


	#-- Read the sample's file names

my %allsamples;
tie %allsamples,"Tie::IxHash";
open(infile1,$mappingfile) || die "Can't open mappingfile $mappingfile\n";
print "  Reading samples from $mappingfile\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	next if(($mode eq "sequential") && ($t[0] ne $projectname));
	if($cleaning) { $userdir=$fastqdir; }
	if($t[2] eq "pair1") { $allsamples{$t[0]}{"$userdir/$t[1]"}=1; } 
	elsif ($t[2] eq "pair2") { $allsamples{$t[0]}{"$userdir/$t[1]"}=2; }
	}
close infile1;

my @f=keys %allsamples;
my $numsamples=$#f+1;
my $nums;
print "  Metagenomes found: $numsamples\n";


        #-- Count contig length from 01.lon file (using SAM headers gives trouble when using minimap2
my %lencontig;
print "  Reading contig length from $contigslen\n";
open(infilelen,$contigslen) || die "Cannot read contig length from $contigslen\n";
while(<infilelen>) {
        chomp;
        next if !$_;
        my($contigid,$tlen)=split(/\t/,$_);
        $lencontig{$contigid}=$tlen;
        }
close infilelen;


        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

my $allsamplesmapped = 1;
foreach my $thissample(keys %allsamples) {
	my $bamfile="$bamdir/$projectname.$thissample.bam";
	my $baifile="$bamdir/$projectname.$thissample.bam.bai";
	if(! -e $bamfile or ! -e $baifile) { $allsamplesmapped = 0; }
	}

if($mapper eq "bowtie") {
	print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
	print outmet "Read mapping against contigs was performed using Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n"; 
	if( (-e "$bowtieref.1.bt2" or $allsamplesmapped) and !$force_overwrite) {
		print "  Mapping reference is present or all samples are already mapped, skipping indexing\n";
		}
        else {
        	print("  Creating reference from contigs\n");
                my $bowtie_command="$bowtie2_build_soft --quiet $contigsfna $bowtieref";
		print outsyslog "Creating Bowtie reference: $bowtie_command\n";
		system($bowtie_command);
                }
        }
elsif($mapper eq "bwa") {
	print "  Mapping with BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
	print outmet "Read mapping against contigs was performed using BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
        if( (-e "$bowtieref.bwt" or $allsamplesmapped ) and !$force_overwrite) {
		print "  Mapping reference is present or all samples are already mapped, skipping indexing\n";
		}
	else {
        	print("Creating reference.\n");
                my $bwa_command="$bwa_soft index -p $bowtieref $contigsfna";
		print outsyslog "Creating BWA reference: $bwa_command\n";
                system($bwa_command);
                }
        }
elsif($mapper=~/minimap/i) { 
	print "  Mapping with Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	print outmet "Read mapping against contigs was performed using Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	}


	#-- Prepare output files

if(-e $contigcov) { system("rm $contigcov"); }
open(outfile1,">$mappingstat") || die "Can't open mappingstat file $mappingstat for writing\n";	#-- File containing mapping statistics
print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";
open(outfile3,">$mapcountfile") || die "Can't open mapcount file $mapcountfile for writing\n";
print outfile3 "# Created by $0 from $gff_file, ",scalar localtime,". SORTED TABLE\n";
print outfile3 "Gen\tLength\tReads\tBases\tRPKM\tCoverage\tTPM\tSample\n";

	#-- Now we start mapping the reads of each sample against the reference

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$samfile,$bamfile,$baifile,$formatoption);
	$nums++;
	my (@pair1,@pair2)=();
	print "  Working with sample $nums: $thissample\n";
	foreach my $ifile(sort keys %{ $allsamples{$thissample} }) {
		if(!$formatseq) {
			if($ifile=~/fasta|fa$/) { $formatseq="fasta"; }
			else { $formatseq="fastq"; }
			}
		
	        #-- Get reads from samples
		if($allsamples{$thissample}{$ifile}==1) { push(@pair1,$ifile); } else { push(@pair2,$ifile); }
		}
	my($par1name,$par2name);
	if($pair1[0]=~/gz/) { $par1name="$projectname.$thissample.current_1.gz"; } 
	else { $par1name="$projectname.$thissample.current_1"; }
	if($pair2[0]=~/gz/) { $par2name="$projectname.$thissample.current_2.gz"; }
	else { $par2name="$projectname.$thissample.current_2";}
	my $a1=join(" ",@pair1);	
	if($#pair1>0) {	$command="cat $a1 > $tempdir/$par1name; "; } else { $command="cp $a1 $tempdir/$par1name; "; }	
	if($#pair2>=0) { 
		my $a2=join(" ",@pair2);	
		if($#pair2>0) {	$command.="cat $a2 > $tempdir/$par2name; "; } else { $command.="cp $a2 $tempdir/$par2name; "; }	
		}
	$samfile="$bamdir/$projectname.$thissample.sam";
	$bamfile="$bamdir/$projectname.$thissample.bam";
	$baifile="$bamdir/$projectname.$thissample.bam.bai";
	if(-e $bamfile and -e $baifile and !$force_overwrite) { print "  BAM file already found in $bamfile, skipping\n"; }
	else {
		print "  Getting raw reads\n";
		#print "$command\n";
		print outsyslog "Getting raw reads for $thissample: $command\n";
		system $command; 
	
		#-- Now we start mapping reads against contigs
	
		print "  Aligning to reference with $mapper\n";
		#-- Support for single reads
       		if(!$mapper || ($mapper eq "bowtie")) {
           		if($formatseq eq "fasta") { $formatoption="-f"; }
			if($mapping_options eq "") { $mapping_options = "--very-sensitive-local"; } # very-sensitive-local would interfere with custom mapping options so add it here 
    	    		if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --quiet -p $numthreads -S $samfile $mapping_options"; }
	    		else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name  --quiet -p $numthreads -S $samfile $mapping_options"; } }
        	elsif($mapper eq "bwa") {
            		#Apparently bwa works seamlesly with fasta files as input.
            		if(-e "$tempdir/$par2name") { $command="$bwa_soft mem $bowtieref $tempdir/$par1name $tempdir/$par2name -v 1 -t $numthreads $mapping_options > $samfile"; }
            		else { $command="$bwa_soft mem $bowtieref $tempdir/$par1name -v 1 -t $numthreads $mapping_options > $samfile"; } }
        	elsif($mapper eq "minimap2-ont") {
            		#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            		if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $samfile"; }
            		else { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $samfile"; } }
        	elsif($mapper eq "minimap2-pb") {
            		#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            		if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $samfile"; }
            		else { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $samfile"; } }
        	elsif($mapper eq "minimap2-sr") {
            		#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            		if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $samfile"; }
            		else { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $samfile"; } 
			}
		else { die "Mapper $mapper not implemented!" }
                                  
		#print "$command\n";
		print outsyslog "Aligning with $mapper: $command\n";
        	my $ecode = 0;
		$ecode = system $command;
        	if($ecode!=0)     { die "An error occurred during mapping with command       $command!"; }
		
		#-- Transform to sorted bam
                $command = "$samtools_soft sort $samfile -o $bamfile.temp -@ $numthreads > /dev/null 2>&1";
        	$ecode = system("$command");
		if($ecode!=0) { die "Error running samtools: $command"; }
		system("rm $samfile");
		#-- Add read group tags identifying the sample from which the reads come from
		$command = "$samtools_soft addreplacerg -r ID:$thissample -r SM:$thissample $bamfile.temp -o $bamfile -@ $numthreads > /dev/null 2>&1; rm $bamfile.temp";
		$ecode = system("$command");
		if($ecode!=0) { die "Error running samtools: $command"; }
		#-- Index the bam file
		$command = "$samtools_soft index $bamfile -@ $numthreads > /dev/null 2>&1";
		$ecode = system("$command");
        	if($ecode!=0) { die "Error running samtools: $command"; }
		}

	#-- Calculating contig coverage/RPKM

	my $totalreads=contigcov(\%lencontig,$thissample,$bamfile);
	
	system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files

	}
if($warnmes) { 
	print outfile1 "#\n# Notice that mapping percentage is low (<50%) for some samples. This is a potential problem,  meaning that most reads are not represented in the assembly\n";
	if($mincontiglen>200) { 
		print outfile1 "# Notice also that you set the minimum contig length to $mincontiglen. In this way you are removing the contigs shorter than that size. This can be, at least partially, the cause of this low mapping percentage\n";
		print outfile1 "# It is likely that redoing the analysis with the default minimum contig length (200) will improve the results\n";
		print outfile1 "# If not, you could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n";
		}
	else { print outfile1 "# You could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n"; }
	}

close outfile1;



#----------------- Contig coverage


sub contigcov {
	print "  Calculating contig coverage\n";
	my %lencontig = %{$_[0]};
	my $thissample = $_[1];
	my $bamfile = $_[2];
	my %readcount;
	my($mappedreads,$totalreadcount,$totalreadlength)=0;
	open(outfile4,">>$contigcov") || die "Can't open contigcov file $contigcov for writing\n";

	#-- Count bases mapped from the sam file
	
	my($readid,%seenreads);
	open infile4,"samtools view $bamfile |" || die "Can't open bam file $bamfile\n";
	while(<infile4>) {
		chomp;
		my @t=split(/\t/,$_);
		next if($_=~/^\@/); # not really needed anymore since we don't read the headers with samtools view
	
		#-- Use the mapped reads to sum base coverage

		if($t[5]!~/\*/) { 			#-- If the read mapped, accum reads and bases
			$readid=$t[0];
	                if($mapper=~/minimap2/){        #-- Minimap2 can output more than one alignment per read
        	                if($seenreads{$readid}) { next; }
                	        else { $seenreads{$readid} = 1 }
                	}
			$readcount{$t[2]}{reads}++;
			$readcount{$t[2]}{lon}+=length $t[9];
			$mappedreads++;
		}       
		$totalreadcount++;
		$totalreadlength+=length $t[9];
	}
	close infile4;

	my $mapperc = 0; # avoid divisions by zero if the BAM has no reads (can happen with custom flags)
	if($totalreadcount) { $mapperc=($mappedreads/$totalreadcount)*100; }
	if($mapperc<50) { $warnmes=1; }
	printf outfile1 "$thissample\t$totalreadcount\t$mappedreads\t%.2f\t$totalreadlength\n",$mapperc;		#-- Mapping statistics

	#-- Output RPKM/coverage values

	my $accumrpk;
	my %rp;
	foreach my $rc(sort keys %readcount) { 
		my $longt=$lencontig{$rc};
		$rp{$rc}=$readcount{$rc}{reads}/$longt;
		$accumrpk+=$rp{$rc};
		}
	$accumrpk/=1000000;

	print outfile4 "#-- Created by $0, ",scalar localtime,"\n";
	print outfile4 "# Contig ID\tAv Coverage\tRPKM\tTPM\tContig length\tRaw reads\tRaw bases\tSample\n";
	foreach my $rc(sort keys %readcount) { 
		my $longt=$lencontig{$rc};
		next if(!$longt);
		my $coverage=$readcount{$rc}{lon}/$longt;
		my ($rpkm, $tpm) = (0,0);
		if($totalreadcount) { # avoid divisions by zero if the BAM has no reads (can happen with custom flags)
			$rpkm=($readcount{$rc}{reads}*1000000)/(($longt/1000)*$totalreadcount); #-- Length of contig in Kbs
			$tpm=$rp{$rc}/$accumrpk;
			}
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.2f\t%.1f\t%.1f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm,$tpm; }
		}
	close outfile4;	
	return $totalreadcount;
}



