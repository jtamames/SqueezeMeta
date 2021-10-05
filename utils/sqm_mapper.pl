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
use threads;

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my($symlinkpath,$symlinkdest,$utilsdir);
if(-l __FILE__)
        {
        $symlinkpath = dirname(__FILE__);
        $symlinkdest = readlink(__FILE__);
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
print BOLD "\nSQM_mapper v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;


	#-- Configuration variables from conf file

my($mappingfile,$fastqdir,$gff_file,$mapper,$outdir,$reference,$hel,$project,$numthreads,$funfile);
my $verbose=0;
my $fullmap=1;		#-- Creates a file with the mapping of all reads
my $bowtie2_build_soft = "$installpath/bin/bowtie2/bowtie2-build";
my $bowtie2_x_soft     = "$installpath/bin/bowtie2/bowtie2";
my $bwa_soft           = "$installpath/bin/bwa";
my $minimap2_soft      = "$installpath/bin/minimap2";

my $helpshort="Usage: SQM_mapper.pl -r <reference file> -s <samples file> -f <raw fastq dir> -g <GFF file> [options]\n";

my $helptext = <<END_MESSAGE;
Usage: SQM_mapper.pl -r <reference file> -g <GFF file> -s <samples file> -f <raw fastq dir> -o <output directory> [options]

Mandatory parameters:
   -r|-reference: Reference (meta)genome (REQUIRED)
   -g: GFF file with ORF positions for the reference (REQUIRED)
   -s: Samples file, same format than standard SqueezeMeta (REQUIRED)
   -f|-seq: Read files (fastq/fasta) directory (REQUIRED)
   -o|-outdir: Output directory (REQUIRED)
   
 Options:
   -n|-name: Prefix name for the results (Default: sqm)
   -t: Number of threads (Default: 12)
   -m: Mapper to be used (Bowtie, BWA, minimap2-ont, minimap2-pb, minimap2-sr) (Default: Bowtie)
   -fun: File containing functional annotations for the genes in the reference
   -h: this help

END_MESSAGE

my $result = GetOptions ("t=i" => \$numthreads,
                     "s|samples=s" => \$mappingfile,
                     "f|seq=s" => \$fastqdir,
		     "g=s" => \$gff_file, 
		     "m=s" => \$mapper,
		     "o|outdir=s" => \$outdir,
		     "r|reference=s"  => \$reference,
		     "n|name=s"  => \$project,
		     "fun=s" => \$funfile,
		     "h" => \$hel
		    );

if($hel) { print "$helptext\n"; die; }
if(!$numthreads) { $numthreads=12; }
if(!$mapper) { $mapper="bowtie"; }
if(!$project) { $project="sqm"; }
if(!$reference) { die "Missing argument -r\n $helpshort\n"; }
if(!$gff_file) { die "Missing argument -g\n $helpshort\n"; } 
if(!$fastqdir) { die "Missing argument -f\n $helpshort\n"; } 
if(!$outdir) { die "Missing argument -o\n $helpshort\n"; }
$mapper=~tr/A-Z/a-z/;
if(($mapper ne "bowtie") && ($mapper ne "bwa") && ($mapper ne "minimap2-ont") && ($mapper ne "minimap2-pb") && ($mapper ne "minimap2-sr")) { die "Unrecognized mapper $mapper. Valid options are Bowtie, BWA, minimap2-ont, minimap2-pb, minimap2-sr\n"; } 

if(-e $outdir) { print "WARNING: Output directory $outdir already exists\n"; }
else { system("mkdir $outdir"); }
my $samdir="$outdir/sam";
if(-e $samdir) {} else { system("mkdir $samdir"); }
my $tempdir="$outdir/temp";
if(-e $tempdir) {} else { system("mkdir $tempdir"); }


my $syslogfile="$outdir/syslog";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


	#-- Read the sample's file names

my %allsamples;
tie %allsamples,"Tie::IxHash";
open(infile1,$mappingfile) || die "Can't open mappingfile $mappingfile\n";
print "  Reading mapping file from $mappingfile\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	if($t[2] eq "pair1") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=1; } 
	elsif ($t[2] eq "pair2") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=2; }
	}
close infile1;

my @f=keys %allsamples;
my $numsamples=$#f+1;
my $nums;
print "  Samples found: $numsamples\n";


        #-- Creates Bowtie2, BWA reference for mapping (index the contigs)

my $bowtieref="$tempdir/$project";
if($mapper eq "bowtie") {
	print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
        if(-e "$bowtieref.1.bt2") { print "   Found reference in $bowtieref.1.bt2, Skipping\n"; }
        else {
                my $bowtie_command="$bowtie2_build_soft --quiet $reference $bowtieref";
    	   	print("  Creating reference in $bowtieref\n");
 		print outsyslog "Creating Bowtie reference: $bowtie_command\n";
                system($bowtie_command);
                }
        }
elsif($mapper eq "bwa") {
	print "  Mapping with BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
        if(-e "$bowtieref.bwt") { print "Found reference in $bowtieref.bwt, Skipping\n";}
        else {
        	print("Creating reference in $bowtieref\n");
                my $bwa_command="$bwa_soft index -p $bowtieref $reference";
		print outsyslog "Creating BWA reference: $bwa_command\n";
                system($bwa_command);
                }
        }
elsif($mapper=~/minimap/i) { 
	print "  Mapping with Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	print outsyslog "Read mapping against contigs was performed using Minimap2\n"; 
	}

	#-- Prepare output files

my $mappingstat="$outdir/$project.mappingstat";
my $mapcountfile="$outdir/$project.mapcount";
my $contigcov="$outdir/$project.contigcov";
my $mapfunfile="$outdir/$project.mapcount.fun";
my $fullmapfile="$outdir/$project.fullmap";

open(outfile1,">$mappingstat") || die "Can't open mappingstat file $mappingstat for writing\n";	#-- File containing mapping statistics
print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";
if($fullmap) {
	open(outfile2,">$fullmapfile") || die "Can't open fullmap file $fullmapfile for writing\n";	
	print outfile2 "#-- Created by $0, ",scalar localtime,"\n";
	print outfile2 "# Read\tGene\tBases\n";
	}
	
open(outfile3,">$mapcountfile") || die "Can't open mapcount file $mapcountfile for writing\n";
print outfile3 "# Created by $0 from $gff_file, ",scalar localtime,". SORTED TABLE\n";
print outfile3 "Gen\tTotal length\tReads\tBases\tRPKM\tCoverage\tTPM\tSample\n";
if($funfile) {
	open(outfilefun,">$mapfunfile") || die "Can't open mapcount file $mapcountfile for writing\n";
	print outfilefun "# Created by $0 from $gff_file and $funfile, ",scalar localtime,"\n";
	print outfilefun "Function\tSample\tCopy number\tLength\tReads\tBases\tRPKM\tCoverage\tTPM\n";
	}

	#-- Now we start mapping the reads of each sample against the reference

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$outsam,$formatoption);
	$nums++;
	my (@pair1,@pair2)=();
	print "  Working with sample $nums: $thissample\n";
	foreach my $ifile(sort keys %{ $allsamples{$thissample} }) {
		if(!$formatseq) {
			if($ifile=~/fasta|fa$|fa.gz$/) { $formatseq="fasta"; }
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
	if($#pair1>0) {	$command="cat $a1 > $tempdir/$par1name; "; } else { $command="cp $a1 $tempdir/$par1name; "; }	
	if($#pair2>=0) { 
		my $a2=join(" ",@pair2);	
		if($#pair2>0) {	$command.="cat $a2 > $tempdir/$par2name; "; } else { $command.="cp $a2 $tempdir/$par2name; "; }	
		}
	print "  Getting raw reads\n";
	# print "$command\n";
	print outsyslog "Getting raw reads for $thissample: $command\n";
	system $command; 
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference with $mapper\n";
	$outsam="$samdir/$project.$thissample.sam";
	if(-e $outsam) { print "  SAM file already found in $outsam, skipping\n"; }
	else {
	
		#-- Support for single reads
     	   if(!$mapper || ($mapper eq "bowtie")) {
      	      if($formatseq eq "fasta") { $formatoption="-f"; }
    		    if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --very-sensitive-local --quiet -p $numthreads -S $outsam"; }
	  	  else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name --very-sensitive-local --quiet -p $numthreads -S $outsam"; } }
       	 elsif($mapper eq "bwa") {
       	     #Apparently bwa works seamlesly with fasta files as input.
        	    if(-e "$tempdir/$par2name") { $command="$bwa_soft mem $bowtieref $tempdir/$par1name $tempdir/$par2name -v 1 -t $numthreads > $outsam"; }
        	    else { $command="$bwa_soft mem $bowtieref $tempdir/$par1name -v 1 -t $numthreads > $outsam"; } }
       	 elsif($mapper eq "minimap2-ont") {
        	    #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
         	   if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $reference $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
          	  else { $command="$minimap2_soft -ax map-ont $reference $tempdir/$par1name -t $numthreads > $outsam"; } }
       	 elsif($mapper eq "minimap2-pb") {
       	     #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
        	    if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $reference $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
        	    else { $command="$minimap2_soft -ax map-pb $reference $tempdir/$par1name -t $numthreads > $outsam"; } }
       	 elsif($mapper eq "minimap2-sr") {
         	   #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
         	   if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $reference $tempdir/$par1name $tempdir/$par2name -t $numthreads > $outsam"; }
         	   else { $command="$minimap2_soft -ax sr $reference $tempdir/$par1name -t $numthreads > $outsam"; } }

                                  
		# print "$command\n";
		print outsyslog "Aligning with $mapper: $command\n";
		system($command);
       		my $ecode = 0;
		if(-e $outsam) {} else { $ecode = system $command; }
       		if($ecode!=0)     { die "An error occurred during mapping!"; }
	}

	#-- Calculating contig coverage/RPKM

	  my $totalreads=contigcov($thissample,$outsam);
	
	#-- And then we call the counting
	
	 # system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files
	 print outsyslog "Calling sqm_counter: Sample $thissample, SAM $outsam, Number of reads $totalreads, GFF $gff_file\n";
	 sqm_counter($thissample,$outsam,$totalreads,$gff_file,$funfile); 
	}

close outfile1;

print "  Output in $mapcountfile\n";
close outfile3;
close outfilefun;
# system("rm $bowtieref.*");	#-- Deleting bowtie references
system("rm $tempdir/count.*");

	#-- Sorting the mapcount table is needed for reading it with low memory consumption in step 13
	
my $command="sort -t _ -k 2 -k 3 -n $mapcountfile > $tempdir/mapcount.temp; mv $tempdir/mapcount.temp $mapcountfile";
print outsyslog "Sorting mapcount table: $command\n";
system($command);	


#----------------- sqm_counter counting 

sub sqm_counter {
	print "  Counting with sqm_counter: Opening $numthreads threads\n";
	my($thissample,$samfile,$totalreadcount,$gff_file,$funfile)=@_;
	my(%genesincontigs,%accum,%accumfun,%long_gen,%fun,%funlong,%copynumber);
	my($countreads,$lastread);
	open(infile2,$gff_file) || die "Can't open gff file $gff_file for writing\n";
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		my $posinit=$k[3];
		my $posend=$k[4];
		if($posend<$posinit) { my $tpos=$posinit; $posinit=$posend; $posend=$posinit; }  
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

	#-- Read the function file (optional)
	
	if($funfile) {
		if(-e $funfile) { 
			print "  Reading functions from $funfile\n";
			open(infile2,$funfile) || warn "WARNING: Cannot open function file $funfile!\n";
 			while(<infile2>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @l=split(/\t/,$_);
				$fun{$l[0]}=$l[1];
				$funlong{$l[1]}+=$long_gen{$l[0]};
				$copynumber{$l[1]}++;
				}
			close infile2;
			}
		else { print RED "WARNING: Cannot open function file $funfile!\n"; }
		}

	
	my($tolines,$thread);
	$tolines=int($totalreadcount/$numthreads);

	for(my $thread=1; $thread<=$numthreads; $thread++) {
		my $thr=threads->create(\&current_thread,\%genesincontigs,\%long_gen,$thread,$samfile,$tolines,$totalreadcount,$gff_file);
		}
	$_->join() for threads->list();
	
	for(my $thread=1; $thread<=$numthreads; $thread++) {
		my $provfile="$tempdir/count.$thread";
		open(intemp,$provfile) || warn "Cannot open $provfile\n";
		while(<intemp>) {
			chomp;
			my @li=split(/\t/,$_);
			my $tfun=$fun{$li[0]};
			$accum{$li[0]}{reads}+=$li[1];
			$accum{$li[0]}{bases}+=$li[2];
			$accumfun{$tfun}{reads}+=$li[1];
			$accumfun{$tfun}{bases}+=$li[2];
			}
		close intemp;
		}
	
	my $accumrpk;
	my %rpk;
	foreach my $print(sort keys %accum) { 
		my $longt=$long_gen{$print};
		next if(!$longt);
		$rpk{$print}=$accum{$print}{reads}/$longt;
		$accumrpk+=$rpk{$print};
		}
	$accumrpk/=1000000;

	#-- Reading genes from gff for: 1) include all ORFs, even these with no counts, and 2) in the fixed order

	my $currentgene;
	open(infilegff,$gff_file) || die "Can't open gff file $gff_file\n";
	while(<infilegff>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		if($_=~/ID\=([^;]+)/) { $currentgene=$1; }	
		my $longt=$long_gen{$currentgene};
		next if(!$longt);
		my $coverage=$accum{$currentgene}{bases}/$longt;
		my $rpkm=($accum{$currentgene}{reads}*1000000)/(($longt/1000)*$totalreadcount);  #-- Length of gene in Kbs
		my $tpm=$rpk{$currentgene}/$accumrpk;
		printf outfile3 "$currentgene\t$longt\t$accum{$currentgene}{reads}\t$accum{$currentgene}{bases}\t%.3f\t%.3f\t%.3f\t$thissample\n",$rpkm,$coverage,$tpm;
		}
	close infilegff;
	
	#-- Now calculating for functions
	
	if(-e $funfile) {
		my($accumrpkfun,$mapped);
		my %rpkfun;
		foreach my $print(sort keys %accumfun) { 
			my $longt=$funlong{$print};
			next if(!$longt);
			my $averlong=$longt/$copynumber{$print};
			$rpkfun{$print}=$accumfun{$print}{reads}/$averlong;
			$accumrpkfun+=$rpkfun{$print};
			$mapped+=$accumfun{$print}{reads};
			}
		$accumrpkfun/=1000000;
		foreach my $currentgene(sort keys %accumfun) { 
			my $longt=$funlong{$currentgene}; 
			next if(!$longt); 
			my $coverage=$accumfun{$currentgene}{bases}/$longt;
			my $averlong=$longt/$copynumber{$currentgene};
			my $rpkm=($accumfun{$currentgene}{reads}*1000000)/(($longt/1000)*$totalreadcount);  #-- Length of gene in Kbs
			my $tpm=$rpkfun{$currentgene}/$accumrpkfun;
			printf outfilefun "$currentgene\t$thissample\t$copynumber{$currentgene}\t$longt\t$accumfun{$currentgene}{reads}\t$accumfun{$currentgene}{bases}\t%.3f\t%.3f\t%.3f\n",$rpkm,$coverage,$tpm;
			# printf "$currentgene\t$longt\t$accumfun{$currentgene}{reads}\t$accumfun{$currentgene}{bases}\t%.3f\t%.3f\t%.3f\t$thissample\n",$rpkm,$coverage,$tpm;
			}
		}
		
	}	
	
sub current_thread {	
	my %genesincontigs=%{$_[0]}; shift;
	my %long_gen=%{$_[0]}; shift;
	my $thread=shift;
	my $samfile=shift;
	my $tolines=shift;
	my $totalreadcount=shift;
	my $gff_file=shift;
	my($countreads,$lastread);
	my %accum;
	my $initline=$tolines*($thread-1);
	my $endline=$tolines*($thread);
	# print "   Thread $thread opening $samfile\n";
	open(infile3,$samfile) || die "Can't open sam file $samfile\n"; ;
	while(<infile3>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		next if(($k[0] eq $lastread) && ($mapper=~/minimap2/));       #-- Minimap2 can output more than one alignment per read
		$countreads++;
		last if($countreads>$endline);
		next if($countreads<$initline);
		$lastread=$readid;
		my $cigar=$k[5];
		next if($k[2]=~/\*/);
		next if($cigar eq "*");
		# print "\n$_\n";
		my $initread=$k[3];                     
		my $incontig=$k[2];
		my $endread=$initread;
		#-- Calculation of the length match end using CIGAR string

		$cigar=~s/^\d+S//;
		$cigar=~s/\d+S$//;
		while($cigar=~/^(\d+)([IDMNSHPX])/) {
			my $mod=$1;
			my $type=$2;
			if($type=~/M|D|N/) { $endread+=$mod; }	#-- Update end position according to the match found
			elsif($type=~/I|S/) { $endread-=$mod; }
			$cigar=~s/^(\d+)([IDMNSHPX])//g;
			}
		# print "*$incontig*$init*$end\n";
		
		foreach my $initgen(sort { $a<=>$b; } keys %{ $genesincontigs{$incontig} }) {
			my $basesingen;
			last if($endread<$initgen);
			my($endgen,$genname)=split(/\:/,$genesincontigs{$incontig}{$initgen});		
			# print "  $incontig*$initread-$endread*$initgen-$endgen*\n"; 
			if((($initread>=$initgen) && ($initread<=$endgen)) && (($endread>=$initgen) && ($endread<=$endgen))) {   #-- Read is fully contained in the gene
				$basesingen=$endread-$initread;
				if($verbose) { print "Read within gene: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				if($fullmap) { print outfile2 "$readid\t$genname\t$basesingen\n"; }
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			elsif(($initread>=$initgen) && ($initread<=$endgen)) {   #-- Read starts within this gene
				$basesingen=$endgen-$initread;
				if($fullmap) { print outfile2 "$readid\t$genname\t$basesingen\n"; }
				if($verbose) {  print "Start of read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
 			elsif(($endread>=$initgen) && ($endread<=$endgen)) {   #-- Read ends within this gene
				$basesingen=$endread-$initgen;
				if($verbose) {  print "End of read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				if($fullmap) { print outfile2 "$readid\t$genname\t$basesingen\n"; }
				$accum{$genname}{bases}+=$basesingen;
				$accum{$genname}{reads}++;
				}
			elsif(($initread<=$initgen) && ($endread>=$endgen)) {  #-- Gen is fully contained in the read
				if($verbose) {  print "Gene within read: $readid $initread-$endread $incontig $initgen-$endgen $basesingen\n"; }
				$basesingen=$endgen-$initgen;
				if($fullmap) { print outfile2 "$readid\t$genname\t$basesingen\n"; }
				$accum{$genname}{reads}++;
				$accum{$genname}{bases}+=$basesingen;
				}
			}

		
		}
	close infile3;
	print "  $countreads reads counted\n";
	open(outtemp,">$tempdir/count.$thread") || die;
	foreach my $h(sort keys %accum) { print outtemp "$h\t$accum{$h}{reads}\t$accum{$h}{bases}\n"; }
	close outtemp;
}


#----------------- Contig coverage


sub contigcov {
	print "  Calculating contig coverage\n";
	my($thissample,$outsam)=@_;
	my(%lencontig,%readcount)=();
	my($mappedreads,$totalreadcount,$totalreadlength)=0;
	open(outfile4,">>$contigcov") || die "Can't open contigcov file $contigcov for writing\n";

	#-- Count length of contigs and bases mapped from the sam file

	my($thisr,$lastr);
	open(infile4,$outsam) || die "Can't open $outsam\n"; ;
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
				$thisr=$t[0];
				next if(($thisr eq $lastr) && ($mapper=~/minimap2/));
				$lastr=$thisr;
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
		my $rpkm=($readcount{$rc}{reads}*1000000000)/(($longt/1000)*$totalreadcount); #-- Length of contig in Kbs
		my $tpm=$rp{$rc}/$accumrpk;
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.2f\t%.1f\t%.1f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm,$tpm; }
		}
	return $totalreadcount;
}

