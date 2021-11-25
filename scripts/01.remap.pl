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

our($datapath,$assembler,$outassembly,$megahit_soft,$mappingfile,$assembler_options,$extassembly,$numthreads,$spades_soft,$prinseq_soft,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$methodsfile,$syslogfile,$mincotiglen);
our($bowtieref,$bowtie2_build_soft,$project,$contigsfna,$mappingfile,$mapcountfile,$mode,$resultpath,$contigid,$contigcov,$bowtie2_x_soft,$mapper,$bwa_soft, $minimap2_soft, $gff_file,$tempdir,$numthreads,$scriptdir,$mincontiglen,$doublepass,$gff_file_blastx,$methodsfile,$syslogfile,$keepsam10);

my($seqformat,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name);
my $fastqdir="$datapath/raw_fastq";
my $samdir="$datapath/sam";
my $outsam;

my $samcommand="mkdir $samdir > /dev/null 2>&1";
system $samcommand;
system("rm $interdir/01.$projectname.singletons  > /dev/null 2>&1");
open(outsyslog,">>$syslogfile");

my $assemblybck="$interdir/01.$projectname.assembly_orig.fasta";    #-- Make a copy of the original assembly
print "  Copying contigs from $contigsfna to $assemblybck\n";
system("cp $contigsfna $assemblybck");

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

my(%allsamples,%singletonlist,%mapped_original,%readcount);
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

my $cocount;
open(infile0,$contigsfna) || die;       #-- Count contigs, to start numbering after the last one
while(<infile0>) {
	chomp;
	if($_=~/\>([^ ]+)/) {
	my @l=split(/\_/,$1);
	$cocount=pop @l;
	# print "  Last contig number: $cocount\n";
	}
}
close infile0;

foreach my $thissample(keys %allsamples) {
	my($formatseq,$command,$formatoption);
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
	if($pair1[0]=~/gz/) { $par1name="$projectname.$thissample\_1.gz"; } 
	else { $par1name="$projectname.$thissample\_1"; }
	if($pair2[0]=~/gz/) { $par2name="$projectname.$thissample\_2.gz"; }
	else { $par2name="$projectname.$thissample\_2";}
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
	
	# print "  Aligning to reference with $mapper\n";
	$outsam="$samdir/$projectname.$thissample.sam";
	
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
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads --secondary=no > $outsam"; }
            else { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name -t $numthreads --secondary=no > $outsam"; } }
        elsif($mapper eq "minimap2-pb") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads --secondary=no > $outsam"; }
            else { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name -t $numthreads --secondary=no > $outsam"; } }
        elsif($mapper eq "minimap2-sr") {
            #Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads --secondary=no > $outsam"; }
            else { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name -t $numthreads --secondary=no > $outsam"; } }

	# print "  Mapping: $command\n";
        print outsyslog "  Mapping: $command\n";                          
	system($command);
	
	my $provcontigs="$tempdir/contigs.$nums.fasta";
	my $singletonfile="$interdir/01.$projectname.singletons";
	my $lastid;
	if(!$contigid) { $contigid="$assembler"; }
	# print "  Reading SAM in $outsam to create $provcontigs and $singletonfile\n";
	open(outsingletons,">$provcontigs");
	open(infilesam,$outsam) || die "Cannot open SAM file $outsam\n";
		while(<infilesam>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		next if($readid eq $lastid);	#-- Minimap2 can output more than one alignment per read
		$readcount{$thissample}++;
		if($k[2]=~/\*/) { 
			$cocount++;
			# print "   Contignumber $cocount\n";
                	my $newcontigname="$contigid\_$cocount";
                	print outsingletons ">$newcontigname singleton $readid\n$k[9]\n";
			$singletonlist{$newcontigname}="$thissample\t$readid";
			}
		else { $mapped_original{$thissample}++; }
		$lastid=$readid;
		}
	close infilesam;
	close outsingletons;
 }
my $command="cat $contigsfna $tempdir/contigs.*.fasta > $tempdir/prov.fasta; mv $tempdir/prov.fasta $contigsfna";
system ($command); 

#-- Make original mapping statistics

open(outmap,">$interdir/01.$projectname.mappingstat_orig") || die;
print outmap "#-- Created by $0, ",scalar localtime,"\n";
print outmap "# Sample\tTotal reads\tMapped reads\tMapping perc\n";
foreach my $tsamp(sort keys %readcount) { 
	my $perc=($mapped_original{$tsamp}/$readcount{$tsamp})*100;
	printf outmap "$tsamp\t$readcount{$tsamp}\t$mapped_original{$tsamp}\t%.2f\n",$perc;
	}
print outmap "\n# This result corresponds to the mapping of reads to the original assembly, before the introduction of singletons\n";
close outmap;

#-- Run prinseq_lite for removing short contigs

$command="$prinseq_soft -fasta $contigsfna -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna > /dev/null 2>&1";
print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen \n";
print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$projectname.stats";
print outsyslog "Running prinseq for contig statistics: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }

#-- Counts length of the contigs (we will need it later)

my $numc;
my $singletonfile="$interdir/01.$projectname.singletons";
print "  Counting length of contigs\n";
open(outfile2,">$contigslen") || die "Can't open $contigslen for writing\n";
open(outfile3,">$singletonfile") || die;
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
        		if($singletonlist{$contigname}) { print outfile3 "$contigname\t$singletonlist{$contigname}\n"; }     
		       }
                $seq="";
                $contigname=$thisname;
                }
        else { $seq.=$_;}
        }
close infile2;
if($contigname) { $len=length $seq; print outfile2 "$contigname\t$len\n"; }
close outfile2;
close outfile3;
close outmet;
close outsyslog;

print "  Contigs stored in $contigsfna\n  Number of contigs: $numc\n";
#system("rm $datapath/raw_fastq/par1.$format.gz; rm $datapath/raw_fastq/par2.$format.gz");
system("rm $bowtieref.1.bt2; rm $bowtieref.bwt; rm $outsam");

