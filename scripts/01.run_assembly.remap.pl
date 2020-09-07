#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/08/2018 for version 0.3.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs (currently megahit or spades). Uses prinseq to filter out contigs by length (excluding small ones).

use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 
use Tie::IxHash;

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$megahit_soft,$mappingfile,$assembler_options,$extassembly,$numthreads,$spades_soft,$prinseq_soft,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$methodsfile,$syslogfile);
our($datapath,$bowtieref,$bowtie2_build_soft,$project,$contigsfna,$mappingfile,$mapcountfile,$mode,$resultpath,$contigcov,$bowtie2_x_soft,
    $mapper, $bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$mincontiglen,$doublepass,$gff_file_blastx,$methodsfile,$syslogfile,$keepsam10);

my($seqformat,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name);
my $fastqdir="$datapath/raw_fastq";
my $samdir="$datapath/sam";

if(-e "$datapath/raw_fastq/par1.fastq.gz") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fasta.gz") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta.gz"; $par2name="$datapath/raw_fastq/par2.fasta.gz"; }
elsif(-e "$datapath/raw_fastq/par1.fastq") { $seqformat="fastq"; $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
elsif(-e "$datapath/raw_fastq/par1.fasta") { $seqformat="fasta"; $par1name="$datapath/raw_fastq/par1.fasta"; $par2name="$datapath/raw_fastq/par2.fasta"; }
else { die "Can't find read files in $datapath/raw_fastq\n"; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- trimmomatic commands

if($cleaning) {
	my $orig1=$par1name;
	my $orig2=$par2name;
	$orig1=~s/\.fastq/\.original.fastq/;
	$orig1=~s/\.fasta/\.original.fasta/;
	$orig2=~s/\.fastq/\.original.fastq/;
	$orig2=~s/\.fasta/\.original.fasta/;
	my $tcommand="mv $par1name $orig1; mv $par2name $orig2";
	system $tcommand; 
	if(-e $orig2) { $trimmomatic_command="$trimmomatic_soft PE -threads $numthreads -phred33 $orig1 $orig2 $par1name $par1name.removed $par2name $par2name.removed $cleaningoptions > /dev/null 2>&1"; }
	else { $trimmomatic_command="$trimmomatic_soft SE -threads $numthreads -phred33 $orig1 $par1name $cleaningoptions > /dev/null 2>&1"; }

	if($cleaning) {
		print "  Running trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20) for quality filtering\n";
		print outsyslog "Running trimmomatic: $trimmomatic_command";
		my $ecode = system $trimmomatic_command;
		if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
		print outmet "Quality filtering was done using Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20)\n";
		}
	}

if($extassembly) {
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
	$outassembly=$extassembly; 
	}
else {

	#-- Checks the assembler

	if($assembler=~/megahit/i) { 
		system("rm -r $datapath/megahit > /dev/null 2>&1"); 
		$outassembly="$datapath/megahit/final.contigs.fa";
		if(-e $par2name) { $command="$megahit_soft $assembler_options -1 $par1name -2 $par2name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }
		else {  $command="$megahit_soft $assembler_options -r $par1name -t $numthreads -o $datapath/megahit >> $syslogfile 2>&1"; }  #-- Support for single reads
		print outmet "Assembly was done using Megahit (Li et al 2015, Bioinformatics 31(10):1674-6)\n";
		}
	elsif($assembler=~/spades/i) { 
		system("rm -r $datapath/spades > /dev/null 2>&1"); 
		$outassembly="$datapath/spades/contigs.fasta";
		if(-e $par2name) { $command="$spades_soft --meta --pe1-1 $par1name --pe1-2 $par2name -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile 2>&1"; }
		else { $command="$spades_soft --meta --s1 $par1name  -m 400 -k 21,33,55,77,99,127 $assembler_options -t $numthreads -o $datapath/spades >> $syslogfile"; } #-- Support for single reads
		print outmet "Assembly was done using SPAdes (Bankevich et al 2012, J Comp Biol 19(5):455-77)\n";
		}
        elsif($assembler=~/canu/i) {
                system("rm -r $datapath/canu/* > /dev/null 2>&1");
                $outassembly="$datapath/canu/contigs.fasta";
                if($canumem eq "NF") {
                        print "  Setting available memory for Canu\n";
                        my %mem=get_mem_info;
                        my $ram=$mem{"MemAvailable"}/(1024*1024);
                        my $ramstr=sprintf('%.2f',$ram);
                        $canumem=sprintf('%.0f',int($ram));
			print "  AVAILABLE (free) RAM memory: $ramstr Gb. We will set canu to use $canumem Gb.\n  You can override this setting using the -canumem option when calling SqueezeMeta.pl\n";
			print outsyslog "canumem set to $canumem (Free Mem $ramstr Gb)\n";
			}
     	   	$command="$canu_soft -p $projectname -d $datapath/canu genomeSize=5m corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=$canumem oeaMemory=$canumem batMemory=$canumem mhapThreads=$numthreads mmapThreads=$numthreads ovlThreads=$numthreads ovbThreads=$numthreads ovsThreads=$numthreads corThreads=$numthreads oeaThreads=$numthreads redThreads=$numthreads batThreads=$numthreads gfaThreads=$numthreads merylThreads=$numthreads $assembler_options -nanopore-raw $par1name > $syslogfile 2>&1; "; 
	   	$command.="mv $datapath/canu/$project.contigs.fasta $outassembly"; 
 		print outmet "Assembly was done using Canu (Koren et al 2017, Genome Res 27(5):722-36)\n";
     	  }


	        elsif($assembler=~/flye/i) {
                system("rm -r $datapath/flye > /dev/null 2>&1");
                $outassembly="$datapath/flye/assembly.fasta";
                print outmet "Assembly was done using Flye (Kolmogorov et al 2019, Nat Biotech 37:540–546)\n";
                $command="/home/tamames/Flye/bin/flye -o $datapath/flye --plasmids --meta --genome-size 5m --threads $numthreads --nano-raw $par1name >> $syslogfile";
                }
        elsif($assembler=~/raven/i) {
                system("rm -r $datapath/raven > /dev/null 2>&1");
                system("mkdir $datapath/raven > /dev/null 2>&1");
                $outassembly="$datapath/raven/contigs.fasta";
                print outmet "Assembly was done using Raven (Vaser and Sikic 2019, BioRxiv 10.1101/656306\n";
                $command="/home/tamames/raven/build/bin/raven --graphical-fragment-assembly $datapath/raven/graph.gfa --threads $numthreads $par1name > $outassembly";
                }
        elsif($assembler) { die "Unrecognized assembler\n"; }
	
	#-- Run assembly

	print "  Running assembly with $assembler\n";
	print outsyslog "Running assembly with $assembler: $command\n";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	}

#-- Run prinseq_lite for removing short contigs

if($extassembly) { system("cp $outassembly $contigsfna"); }
else {
	$command="$prinseq_soft -fasta $outassembly -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna > /dev/null 2>&1";
	print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen \n";
	print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n  ";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	}

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$projectname.stats";
print outsyslog "Running prinseq for contig statistics: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";

#-- Standardization of contig names

print "  Renaming contigs\n";
open(infile1,$contigsfna) || die "Can't open $contigsfna\n";
my $provcontigs="$tempdir/contigs.prov";
open(outfile1,">$provcontigs") || die "Can't open $provcontigs for writing\n";
my $cocount;
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^\>/) {
		$cocount++;
		my $newcontigname="$assembler\_$cocount";
		print outfile1 ">$newcontigname\n";
		}
	else { print outfile1 "$_\n"; }
	}
close infile1;

#---------- Mapping reads

        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

if($mapper eq "bowtie") {
	print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
        if(-e "$bowtieref.1.bt2") {}
        else {
        	print("  Creating reference from contigs\n");
                my $bowtie_command="$bowtie2_build_soft --quiet $contigsfna $bowtieref";
                system($bowtie_command);
                }
        }
elsif($mapper eq "bwa") {
	print "  Mapping with BWA (Li and Durbin 2009, Bioinformatics 25(14), 1754-60)\n"; 
       if(-e "$bowtieref.bwt") {}
        else {
        	print("Creating reference.\n");
                my $bwa_command="$bwa_soft index -p $bowtieref $contigsfna";
                system($bwa_command);
                }
        }
elsif($mapper=~/minimap/i) { 
	print "  Mapping with Minimap2 (Li 2018, Bioinformatics 34(18), 3094-3100)\n"; 
	}

my(%allsamples);
tie %allsamples,"Tie::IxHash";
open(infile1,$mappingfile) || die "Can't open mappingfile $mappingfile\n";
print "  Reading mapping file from $mappingfile\n";
while(<infile1>) {
	chomp;
	next if !$_;
	my @t=split(/\t/,$_);
	next if(($mode eq "sequential") && ($t[0] ne $projectname));
	if($t[2] eq "pair1") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=1; } 
	elsif ($t[2] eq "pair2") { $allsamples{$t[0]}{"$fastqdir/$t[1]"}=2; }
	}
close infile1;

my @f=keys %allsamples;
my $numsamples=$#f+1;
my $nums;

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$outsam,$formatoption);
	$nums++;
	my (@pair1,@pair2)=();
	print "  Working with sample $nums: $thissample\n";
	foreach my $ifile(sort keys %{ $allsamples{$thissample} }) {
		if(!$formatseq) {
			if($ifile=~/fasta/) { $formatseq="fasta"; }
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
	$command="cat $a1 > $tempdir/$par1name; ";	
	if($#pair2>=0) { 
		my $a2=join(" ",@pair2);	
		$command.="cat $a2 > $tempdir/$par2name;";	
		}
	print "  Getting raw reads\n";
	# print "$command\n";
	print outsyslog "Getting raw reads for $thissample: $command\n";
	system $command; 
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference with $mapper\n";
	system("mkdir $samdir");
	if($keepsam10) { $outsam="$samdir/$projectname.$thissample.sam"; } else { $outsam="$samdir/$projectname.$thissample.current.sam"; }
	
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

        print "****$command***\n";                          
	system($command);
system("cp $provcontigs $provcontigs.old");
	open(outsingletons,">>$tempdir/singletons.fasta");
	open(infilesam,$outsam) || die;
		while(<infilesam>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@SQ/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		if($k[2]=~/\*/) { 
			print outsingletons ">$k[0]\n$k[9]\n";
			$cocount++;
                	my $newcontigname="$assembler\_$cocount $readid";
                	print outfile1 ">$newcontigname\n$k[9]\n";
			}
		}
	close infilesam;
	close outsingletons;

 }
close outfile1;
system("mv $provcontigs $contigsfna");

#-- Counts length of the contigs (we will need it later)

my $numc;
print "  Counting length of contigs\n";
open(outfile2,">$contigslen") || die "Can't open $contigslen for writing\n";
open(infile2,$contigsfna) || die "Can't open $contigsfna\n";
while(<infile2>) {
        chomp;
        next if !$_;
        if($_=~/^\>([^ ]+)/) {
                $numc++;
                $thisname=$1;
                if($contigname) {
                        $len=length $seq;
                        print outfile2 "$contigname\t$len\n";
                        }
                $seq="";
                $contigname=$thisname;
                }
        else { $seq.=$_;}
        }
close infile2;
if($contigname) { $len=length $seq; print outfile2 "$contigname\t$len\n"; }
close outfile2;
close outmet;
close outsyslog;

print "  Contigs stored in $contigsfna\n  Number of contigs: $numc\n";
#system("rm $datapath/raw_fastq/par1.$format.gz; rm $datapath/raw_fastq/par2.$format.gz");


