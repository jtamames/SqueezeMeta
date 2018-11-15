#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage/RPKM for genes/contigs by mapping back reads to the contigs and count how many fall in each gene/contig
#-- Uses bowtie2 for mapping, and bedtools for counting. 
#-- WARNING! Bedtools version must be <0.24!

$|=1;

use strict;
use Cwd;

my $pwd=cwd();
my $project=$ARGV[0];

do "$project/SqueezeMeta_conf.pl";

	#-- Configuration variables from conf file

our($datapath,$bowtieref,$bowtie2_build_soft,$contigsfna,$mappingfile,$mode,$resultpath,$rpkmfile,$contigcov,$coveragefile,$bowtie2_x_soft,
    $mapper, $bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$bedtools_soft);

my $keepsam=1;  #-- Set to one, it keeps SAM files. Set to zero, it deletes them when no longer needed

my $fastqdir="$datapath/raw_fastq";
my $samdir="$datapath/sam";

if(-d $samdir) {} else { system("mkdir $samdir"); }

 
	#-- Read the sample's file names

my %allsamples;
open(infile1,$mappingfile) || die "Cannot find mappingfile $mappingfile\n";
print "Reading mapping file from $mappingfile\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	next if(($mode eq "sequential") && ($t[0] ne $project));
	if($t[2] eq "pair1") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=1; } 
	elsif ($t[2] eq "pair2") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=2; }
	}
close infile1;

my @f=keys %allsamples;
my $numsamples=$#f+1;
my $nums;
print "Metagenomes found: $numsamples\n";


        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

if($mapper eq "bowtie") {
        if(-e "$bowtieref.1.bt2") {}
        else {
        	print("Creating reference.\n");
                my $bowtie_command="$bowtie2_build_soft --quiet $contigsfna $bowtieref";
                system($bowtie_command);
                }
        }
elsif($mapper eq "bwa") {
        if(-e "$bowtieref.bwt") {}
        else {
        	print("Creating reference.\n");
                my $bwa_command="$bwa_soft index -p $bowtieref $contigsfna";
                system($bwa_command);
                }
        }


	#-- Prepare output files

if(-e "$resultpath/09.$project.rpkm") { system("rm $resultpath/09.$project.rpkm"); }
if(-e $rpkmfile) { system("rm $rpkmfile"); }
if(-e $contigcov) { system("rm $contigcov"); }
open(outfile1,">$resultpath/09.$project.mappingstat") || die;	#-- File containing mapping statistics
print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";

	#-- Now we start mapping the reads of each sample against the reference

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$outsam,$formatoption);
	$nums++;
	my (@pair1,@pair2)=();
	print "Working with $nums: $thissample\n";
	foreach my $ifile(sort keys %{ $allsamples{$thissample} }) {
		if(!$formatseq) {
			if($ifile=~/fasta/) { $formatseq="fasta"; }
			else { $formatseq="fastq"; }
		}
		
	#-- Get reads from samples
		
	if($allsamples{$thissample}{$ifile}==1) { push(@pair1,$ifile); } else { push(@pair2,$ifile); }
	}
	my($par1name,$par2name);
	if($pair1[0]=~/gz/) { $par1name="$project.$thissample.current_1.gz"; } 
	else { $par1name="$project.$thissample.current_1"; }
	if($pair2[0]=~/gz/) { $par2name="$project.$thissample.current_2.gz"; }
	else { $par2name="$project.$thissample.current_2";}
	if($#pair1==0) { $command="ln -s $pair1[0] $tempdir/$par1name;"; } 
	elsif($#pair1>0) { 
		my $a1=join(" ",@pair1);					
		$command="cat $a1 > $tempdir/$par1name; ";	
		}
	if($#pair2==0) { $command.="ln -s $pair2[0] $tempdir/$par2name;"; }
	elsif($#pair2>0) { 
		my $a2=join(" ",@pair2);	
		$command.="cat $a2 > $tempdir/$par2name; cat $a2 > $tempdir/$par2name;";	
		}
	print "  Getting raw reads\n";
	print "$command\n";
	system $command;
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference...\n";
	if($keepsam) { $outsam="$samdir/$project.$thissample.sam"; } else { $outsam="$samdir/$project.$thissample.current.sam"; }
	
	#-- Support for single reads
        if($mapper eq "bowtie") {
            if($formatseq eq "fasta") { $formatoption="-f"; }
    	    if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --quiet -p $numthreads -S $outsam"; }
	    else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name --quiet -p $numthreads -S $outsam"; } }
        elsif($mapper eq "bwa") {
            #Apparently bwa works seemlesly with fasta files as imput.
            if(-e "$tempdir/$par2name") { $command="$bwa_soft mem $bowtieref $tempdir/$par1name $tempdir/$par2name -v 1 -t $numthreads > $outsam"; }
            else { $command="$bwa_soft mem $bowtieref $tempdir/$par1name -v 1 -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-ont") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-pb") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }
        elsif($mapper eq "minimap2-sr") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
            else { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name -t $numthreads > $outsam"; } }

                                  
	print "$command\n";
	system $command;
	
	#-- And then we call bedtools for counting
	
	# htseq();
	system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files
	bedtools($thissample,$outsam);
	contigcov($thissample,$outsam);
}
close outfile1;
system("rm $samdir/current.sam");   


#----------------- htseq counting (deprecated)

#sub htseq {
#	print "  Counting with HTSeq\n";
#	my $command="htseq-count -m intersection-nonempty -s no -t CDS -i \"ID\" $outsam $gff_file > $project.$thissample.current.htseq";
#	system $command;
#	print "  Calculating RPKM from HTseq\n";
#	$command="perl $scriptdir/09.rpkm.pl $project.$thissample.current.htseq $gff_file $thissample >> $resultpath/06.$project.rpkm";
#	system $command;
#}

#----------------- bedtools counting 

sub bedtools {
	print "  Counting with Bedtools\n";
	my($thissample,$outsam)=@_;

	#-- Creating reference for bedtools from the gff file
	#-- Reference has the format: <contig_id> <gen init pos> <gen end pos> <gen ID>

	open(infile2,$gff_file) || die "Cannot find gff file $gff_file\n";  
	# print "Reading gff file from $gff_file\n";
	my $bedreference=$gff_file;
	$bedreference=~s/gff/refbed/;
	print "    Generating reference: $bedreference\n";
	
	open(outfile2,">$bedreference") || die;
	while(<infile2>) {
		my $gid;
		chomp;
		next if(!$_ || ($_=~/^\#/));
		if($_=~/ID\=([^;]+)\;/) { $gid=$1; }	#-- Orf's ID
		my @k=split(/\t/,$_);
		print outfile2 "$k[0]\t$k[3]\t$k[4]\t$gid\n";		# <contig_id> <gen init pos> <gen end pos> <gen ID>
		}
	close infile2;
	close outfile2;

	#-- Creating bedfile from the sam file. It has the format <read id> <init pos match> <end pos match>

	my $bedfile="$tempdir/$project.$thissample.current.bed";
	print "    Generating Bed file: $bedfile\n";
	open(outfile3,">$bedfile") || die;
	open(infile3,$outsam) || die;

	#-- Reading sam file

	while(<infile3>) {
		next if($_=~/^\@/);
		my @k=split(/\t/,$_);
		next if($k[2]=~/\*/);
		my $cigar=$k[5];                       
		my $end=$k[3];

		#-- Calculation of the length match end using CIGAR string

		while($cigar=~/^(\d+)([IDM])/) {
			my $mod=$1;
			my $type=$2;
			if($type=~/M|D/) { $end+=$mod; }	#-- Update end position according to the match found
			elsif($type=~/I/) { $end-=$mod; }
			$cigar=~s/^(\d+)([IDM])//g;
			}
		print outfile3 "$k[2]\t$k[3]\t$end\n";		#-- <read id> <init pos match> <end pos match>
		}
	close infile3;
	close outfile3;

	#-- Call bedtools for counting reads
	
	my $command="$bedtools_soft coverage -a $bedfile -b $bedreference > $tempdir/$project.$thissample.current.bedcount";
	print "    Counting reads: $command\n";
	system $command;	

	#-- Call bedtools for counting bases

	$command="$bedtools_soft coverage -a $bedfile -b $bedreference -d > $tempdir/$project.$thissample.currentperbase.bedcount";
	print "    Counting bases: $command\n";
	system $command;
	
	#-- Run RPKM calculation (rpkm.pl)
	
	print "  Calculating RPKM from Bedtools\n";
	$command="perl $scriptdir/09.rpkm.pl $tempdir/$project.$thissample.current.bedcount $gff_file $thissample >> $rpkmfile";
	system $command;

	#-- Run coverage calculation (coverage.pl)	

	print "  Calculating Coverage from Bedtools\n";
	$command="perl $scriptdir/09.coverage.pl $tempdir/$project.$thissample.currentperbase.bedcount $gff_file $thissample >> $coveragefile";
	system $command;
	
	#-- Remove files
	
	print "  Removing files\n";
	#system("rm $tempdir/$project.$thissample.current_1.fastq.gz");
	#system("rm $tempdir/$project.$thissample.current_2.fastq.gz");
	system("rm $tempdir/$project.$thissample.current.bedcount");
	system("rm $tempdir/$project.$thissample.currentperbase.bedcount"); 
	system("rm $tempdir/$project.$thissample.current.bed"); 
}


#----------------- Contig coverage


sub contigcov {
	print "  Calculating contig coverage\n";
	my($thissample,$outsam)=@_;
	my(%lencontig,%readcount)=();
	my($mappedreads,$totalreadcount,$totalreadlength)=0;
	open(outfile4,">>$contigcov") || die;

	#-- Count length of contigs and bases mapped from the sam file

	open(infile4,$outsam);
	while(<infile4>) {
		chomp;
		my @t=split(/\t/,$_);

		#-- Use the headers to extract contig length

		if($_=~/^\@/) {
		$t[1]=~s/SN\://;
		$t[2]=~s/LN\://;
		$lencontig{$t[1]}=$t[2];
		}
	
		#-- And the mapped reads to sum base coverage

		else {
			if($t[2]!~/\*/) { 			#-- If the read mapped, accum reads and bases
				$readcount{$t[2]}{reads}++;
				$readcount{$t[2]}{lon}+=length $t[9];
				$mappedreads++;
			}       
			$totalreadcount++;
			$totalreadlength+=length $t[9];
		} 
	}
	close infile4;
	
	my $mapperc=($mappedreads/$totalreadcount)*100;
	printf outfile1 "$thissample\t$totalreadcount\t$mappedreads\t%.2f\t$totalreadlength\n",$mapperc;		#-- Mapping statistics

	#-- Output RPKM/coverage values

	print outfile4 "#-- Created by $0, ",scalar localtime,"\n";
	print outfile4 "# Contig ID\tAv Coverage\tRPKM\tContig length\tRaw reads\tRaw bases\tSample\n";
	foreach my $rc(sort keys %readcount) { 
		my $longt=$lencontig{$rc};
		next if(!$longt);
		my $coverage=$readcount{$rc}{lon}/$longt;
		my $rpkm=($readcount{$rc}{reads}*1000000000)/($longt*$totalreadcount);
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.3f\t%.3f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm; }
		}
	close outfile4;	
}

