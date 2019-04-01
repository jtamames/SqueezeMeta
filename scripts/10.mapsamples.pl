#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 28/01/2019 for version 0.4.3, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage/RPKM for genes/contigs by mapping back reads to the contigs and count how many fall in each gene/contig
#-- Uses bowtie2 for mapping, and sqmapper for counting. 

$|=1;

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 

do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

	#-- Configuration variables from conf file

our($datapath,$bowtieref,$bowtie2_build_soft,$contigsfna,$mappingfile,$mapcountfile,$mode,$resultpath,$contigcov,$bowtie2_x_soft,
    $mapper, $bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$doublepass,$gff_file_blastx,$keepsam10);

my $verbose=0;

my $fastqdir="$datapath/raw_fastq";
my $samdir="$datapath/sam";

my $outfile="$resultpath/10.$project.mapcount";

if(-d $samdir) {} else { system("mkdir $samdir"); }

if($doublepass) { $gff_file=$gff_file_blastx; }

	#-- Read the sample's file names

my(%allsamples,%rpk);
tie %allsamples,"Tie::IxHash";
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

#if(-e "$resultpath/09.$project.rpkm") { system("rm $resultpath/09.$project.rpkm"); }
#if(-e $rpkmfile) { system("rm $rpkmfile"); }
if(-e $contigcov) { system("rm $contigcov"); }
open(outfile1,">$resultpath/10.$project.mappingstat") || die;	#-- File containing mapping statistics
print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";
open(outfile3,">$mapcountfile") || die;
print outfile3 "# Created by $0 from $gff_file, ",scalar localtime,"\n";
print outfile3 "Gen\tLength\tReads\tBases\tRPKM\tCoverage\tTPM\tSample\n";

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
	my $a1=join(" ",@pair1);					
	$command="cat $a1 > $tempdir/$par1name; ";	
	if($#pair2>=0) { 
		my $a2=join(" ",@pair2);	
		$command.="cat $a2 > $tempdir/$par2name;";	
		}
	print "  Getting raw reads\n";
	# print "$command\n";
	system $command; 
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference with $mapper\n";
	if($keepsam10) { $outsam="$samdir/$project.$thissample.sam"; } else { $outsam="$samdir/$project.$thissample.current.sam"; }
	
	#-- Support for single reads
        if(!$mapper || ($mapper eq "bowtie")) {
            if($formatseq eq "fasta") { $formatoption="-f"; }
    	    if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --quiet -p $numthreads -S $outsam"; }
	    else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name --quiet -p $numthreads -S $outsam"; } }
        elsif($mapper eq "bwa") {
            #Apparently bwa works seamlesly with fasta files as input.
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

                                  
	# print "$command\n";
	system($command);
        my $ecode = 0;
	#if(-e $outsam) {} else { $ecode = system $command; }
        #if($ecode!=0)     { die "An error occurred during mapping!"; }

	#-- Calculating contig coverage/RPKM

	 my $totalreads=contigcov($thissample,$outsam);
	
	#-- And then we call the counting
	
	 system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files
	 sqm_counter($thissample,$outsam,$totalreads,$gff_file); 
}
close outfile1;

print "Output in $mapcountfile\n";
close outfile3;


#----------------- sqm_counter counting 

sub sqm_counter {
	print "  Counting with sqm_counter\n";
	my($thissample,$samfile,$totalreadcount,$gff_file)=@_;
	my(%genesincontigs,%accum,%long_gen);
	my $countreads;
	open(infile2,$gff_file) || die;
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		my $posinit=$k[3];
		my $posend=$k[4];
		my $genid;
  		my @e=split(/\;/,$k[8]);
 		my @n=split(/\_/,$e[0]);
  		my $ord=$n[$#n];
 		$genid="$k[0]\_$ord";
		$genesincontigs{$k[0]}{$posinit}="$posend:$genid";
		# print "$k[0]*$posinit*$posend*$genid\n";
		$long_gen{$genid}=$posend-$posinit+1;
		}
	close infile2;

	open(infile3,$samfile) || die "Cannot open sam file $samfile\n"; ;
	while(<infile3>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@SQ/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		next if($k[2]=~/\*/);
		my $cigar=$k[5];
		next if($cigar eq "*");
		# print "\n$_\n";
		my $initread=$k[3];                     
		my $incontig=$k[2];
		my $endread=$initread;
		$countreads++;
		if($countreads%1000000==0) { print "    $countreads reads counted\r"; }
		#-- Calculation of the length match end using CIGAR string

		while($cigar=~/^(\d+)([IDM])/) {
			my $mod=$1;
			my $type=$2;
			if($type=~/M|D/) { $endread+=$mod; }	#-- Update end position according to the match found
			elsif($type=~/I/) { $endread-=$mod; }
			$cigar=~s/^(\d+)([IDM])//g;
			}
		# print "*$incontig*$init*$end\n";
		foreach my $initgen(sort { $a<=>$b; } keys %{ $genesincontigs{$incontig} }) {
			my $basesingen;
			last if($endread<$initgen);
			my($endgen,$genname)=split(/\:/,$genesincontigs{$incontig}{$initgen});		
			# print "  $incontig*$initread-$endread*$initgen-$endgen*\n"; 
			if((($initread>=$initgen) && ($initread<=$endgen)) && (($endread>=$initgen) && ($endread<=$endgen))) {   #-- El read esta contenido en el gen
				$basesingen=$endread-$initread;
				if($verbose) { print "Read contenido: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			elsif(($initread>=$initgen) && ($initread<=$endgen)) {   #-- El read empieza dentro de este gen
				$basesingen=$endgen-$initread;
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				if($verbose) {  print "Inicio read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
 			elsif(($endread>=$initgen) && ($endread<=$endgen)) {   #-- El read termina dentro de este gen
				$basesingen=$endread-$initgen;
				if($verbose) {  print "Final read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{bases}+=$basesingen;
				$accum{$genname}{reads}++;
				}
			elsif(($initread<=$initgen) && ($endread>=$endgen)) {  #-- El gen esta contenido en el read
				if($verbose) {  print "Gen contenido: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$basesingen=$endgen-$initgen;
				# print outfile2 "$readid\t$genname\t$basesingen\n";
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			}

		
		}
	close infile3;
	print "\n";

	my $accumrpk;
	foreach my $print(sort keys %accum) { 
		my $longt=$long_gen{$print};
		$rpk{$print}=$accum{$print}{reads}/$longt;
		$accumrpk+=$rpk{$print};
		}
	$accumrpk/=1000000;

	foreach my $print(sort keys %accum) { 
		my $longt=$long_gen{$print};
		my $coverage=$accum{$print}{bases}/$longt;
		my $rpkm=($accum{$print}{reads}*1000000)/($longt*$totalreadcount);
		my $tpm=$rpk{$print}/$accumrpk;
		printf outfile3 "$print\t$longt\t$accum{$print}{reads}\t$accum{$print}{bases}\t%.3f\t%.3f\t%.3f\t$thissample\n",$rpkm,$coverage,$tpm; 
		}
}


#----------------- Contig coverage


sub contigcov {
	print "  Calculating contig coverage\n";
	my($thissample,$outsam)=@_;
	my(%lencontig,%readcount)=();
	my($mappedreads,$totalreadcount,$totalreadlength)=0;
	open(outfile4,">>$contigcov") || die;

	#-- Count length of contigs and bases mapped from the sam file

	open(infile4,$outsam) || die "Cannot open $outsam\n"; ;
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
		my $rpkm=($readcount{$rc}{reads}*1000000000)/($longt*$totalreadcount);
		my $tpm=$rp{$rc}/$accumrpk;
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.3f\t%.3f\t%.3f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm,$tpm; }
		}
	close outfile4;	
	return $totalreadcount;
}

