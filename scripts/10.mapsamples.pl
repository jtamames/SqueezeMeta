#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 28/01/2019 for version 0.4.3, (c) Javier Tamames, CNB-CSIC
#-- Calculates coverage/RPKM for genes/contigs by mapping back reads to the contigs and count how many fall in each gene/contig
#-- Uses bowtie2 for mapping, and sqmapper for counting. 

$|=1;

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";
use threads;

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

if($doublepass) { $gff_file=$gff_file_blastx; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";


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


        #-- Get ORF info from gff file
my(%genesincontigs,%genespos,%long_gen,@genes_ordered);
print "  Reading orf info from $gff_file\n";
open(infile2,$gff_file) || die "Can't open gff file $gff_file for reading\n";
while(<infile2>) {
        chomp;
        next if(!$_ || ($_=~/^\#/));
        my @k=split(/\t/,$_);
	my $contigid=$k[0];
        my $initgen=$k[3];
        my $endgen=$k[4];
        if($endgen<$initgen) { my $tpos=$initgen; $initgen=$endgen; $endgen=$initgen; }
        my $genid;
        my @e=split(/\;/,$k[8]);
        my @n=split(/\_/,$e[0]);
        my $ord=$n[$#n];
        $genid="$contigid\_$ord";
	if(defined $genespos{$genid}) { next; } # avoid counting duplicated entries in the gff file twice
        push @genes_ordered, $genid;
        if(not defined $genesincontigs{$contigid}) {
                my @arr;
                $genesincontigs{$contigid} = [ @arr ];
                }
        push @{$genesincontigs{$contigid}}, $genid;
        $genespos{$genid} = [$initgen, $endgen];
        $long_gen{$genid}=$endgen-$initgen+1;
        }

	#-- Split the contigs into $numthreads chunks and create bed files for them
my @contigchunks; # create 2D array
for(my $thread=1; $thread<=$numthreads; $thread++) {
        my @chunk;
        push @contigchunks, [@chunk];
}

my $i = 0; # populate each subarray with contigs
foreach my $contigid (keys %lencontig) {
	push @{ @contigchunks[$i % $numthreads] }, $contigid;
	$i++;
	}

my @bed_chunk_files_contigs;
my @bed_chunk_files_orfs;
my @count_chunk_files;
for(my $thread=1; $thread<=$numthreads; $thread++) { # write each subarray to a different bed file
	my $bedfilectgs = "$tempdir/count.$thread.contigs.bed";
	my $bedfileorfs = "$tempdir/count.$thread.orfs.bed";
	my $countfile = "$tempdir/count.$thread";
	push @bed_chunk_files_contigs, $bedfilectgs;
	push @bed_chunk_files_orfs, $bedfileorfs;
	push @count_chunk_files, $countfile;
	open(outbedctgs, ">$bedfilectgs") || die "Cannot open $bedfilectgs for writing\n";
	open(outbedorfs, ">$bedfileorfs") || die "Cannot open $bedfileorfs for writing\n";
	foreach my $contigid (@{$contigchunks[$thread-1]}) { # our array is zero-indexed
		print outbedctgs "$contigid\t1\t$lencontig{$contigid}\n";
		foreach my $genid (@{$genesincontigs{$contigid}}) {
			my($initgen, $endgen) = @{$genespos{$genid}};
			print outbedorfs "$contigid\t$initgen\t$endgen\t$genid\n";	
			}
		}
	close outbedctgs;
	close outbedorfs;
	}
	undef @contigchunks;
	undef %genesincontigs;
	undef %genespos;


        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

if($mapper eq "bowtie") {
	print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
	print outmet "Read mapping against contigs was performed using Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n"; 
	#if(-e "$bowtieref.1.bt2") { print "  Found reference in $bowtieref.1.bt2, skipping\n"; } #This will only trigger if index building was incomplete, which beats the purpose
	if(0){}
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
	#if(-e "$bowtieref.bwt") { print "Found reference in $bowtieref.bwt, Skipping\n"; }
        if(0){}
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

#if(-e "$resultpath/09.$project.rpkm") { system("rm $resultpath/09.$project.rpkm"); }
#if(-e $rpkmfile) { system("rm $rpkmfile"); }
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
	print "  Getting raw reads\n";
	#print "$command\n";
	print outsyslog "Getting raw reads for $thissample: $command\n";
	system $command; 
	
	#-- Now we start mapping reads against contigs
	
	print "  Aligning to reference with $mapper\n";
	$samfile="$bamdir/$projectname.$thissample.sam";
	$bamfile="$bamdir/$projectname.$thissample.bam";
	$baifile="$bamdir/$projectname.$thissample.bam.bai";
	if(-e $bamfile and -e $baifile and !$force_overwrite) { print "  BAM file already found in $bamfile, skipping\n"; }
	else {

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

        	my $ecode = system("$samtools_soft sort $samfile -o $bamfile -@ $numthreads; $samtools_soft index $bamfile -@ $numthreads > /dev/null 2>&1");
        	if($ecode!=0) { die "Error running samtools"; }
	 	system("rm $samfile");
	}

	#-- Calculating contig coverage/RPKM

	my $totalreads=contigcov(\%lencontig,$thissample,$bamfile);
	
	#-- And then we call the counting
	
	system("rm $tempdir/$par1name $tempdir/$par2name");   #-- Delete unnecessary files
	print outsyslog "Calling sqm_counter: Sample $thissample, BAM $bamfile, Number of reads $totalreads\n";
	sqm_counter(\@contigchunks,\@genes_ordered,\%long_gen,$thissample,$bamfile,\@bed_chunk_files_contigs,\@count_chunk_files,$totalreads,$mapper,$verbose);


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

print "  Output in $mapcountfile\n";
close outfile3;
if($mapper eq "bowtie" or $mapper eq "bwa") {
	system("rm $bowtieref.*");	#-- Deleting bowtie references
}
system("rm $tempdir/count.*");

	#-- Sorting the mapcount table is needed for reading it with low memory consumption in step 13
	
my $command="sort -T $tempdir -t _ -k 2 -k 3 -n $mapcountfile > $tempdir/mapcount.temp; mv $tempdir/mapcount.temp $mapcountfile";
print outsyslog "Sorting mapcount table: $command\n";
system($command);	



#----------------- sqm_counter counting 

sub sqm_counter {
	print "  Counting with sqm_counter: Opening $numthreads threads\n";
	my($contigchunks,$genes_ordered,$long_gen,$thissample,$bamfile,$bed_chunk_files_contigs,$count_chunk_files,$totalreadcount,$mapper)=@_;
	@contigchunks = @{$contigchunks};
	@genes_ordered = @{$genes_ordered};
	%long_gen = %{$long_gen};
	@bed_chunk_files_contigs = @{$bed_chunk_files_contigs};
	@count_chunk_files = @{$count_chunk_files};
	my %accum;

	my ($bedfilectgs, $bedfileorfs, $countfile);

	my $use_fork=0;

	for(my $thread=1; $thread<=$numthreads; $thread++) {
		$bedfilectgs = @bed_chunk_files_contigs[$thread-1]; # zero indexing
		$bedfileorfs = @bed_chunk_files_orfs[$thread-1];
		$countfile = @count_chunk_files[$thread-1];
		if($use_fork) {
			my $pid = fork;
			die if not defined $pid;
			if (not $pid) { # $pid will be 0 if this is a child process, the parent will don't run this
				#system("perl gencount.pl $thread $bamfile $bedfilectgs $bedfileorfs $countfile $mapper '$samtools_soft' $verbose");
                        	current_thread($thread,$bamfile,$bedfilectgs,$bedfileorfs,$countfile,$mapper,$samtools_soft,$verbose);
				exit;
                        	}
			}
		else {
			my $thr=threads->create(\&current_thread,$thread,$bamfile,$bedfilectgs,$bedfileorfs,$countfile,$mapper,$samtools_soft,$verbose);
			}	
		}

	# Thread/children cleanup.
	if($use_fork) {
		for(my $thread=1; $thread<=$numthreads; $thread++) {
			my $finished = wait(); # wait till all the children have stopped
			}
		}
	else { $_->join() for threads->list(); }


	# Process results	
	for(my $thread=1; $thread<=$numthreads; $thread++) {
		my $countfile = @count_chunk_files[$thread-1];
		open(intemp,$countfile) || warn "Cannot open $countfile\n";
		while(<intemp>) {
			chomp;
			my @li=split(/\t/,$_);
			$accum{$li[0]}{reads}+=$li[1];
			$accum{$li[0]}{bases}+=$li[2];
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

	#-- Go through genes as they appeared in the gff for: 1) include all ORFs, even these with no counts, and 2) in the fixed order

	foreach my $currentgene (@genes_ordered) {
		my $longt=$long_gen{$currentgene};
		next if(!$longt);
		my $coverage=$accum{$currentgene}{bases}/$longt;
		my $rpkm=($accum{$currentgene}{reads}*1000000)/(($longt/1000)*$totalreadcount);  #-- Length of gene in Kbs
		my $tpm=$rpk{$currentgene}/$accumrpk;
		printf outfile3 "$currentgene\t$longt\t$accum{$currentgene}{reads}\t$accum{$currentgene}{bases}\t%.3f\t%.3f\t%.3f\t$thissample\n",$rpkm,$coverage,$tpm;
		}
	close infilegff;

}

sub current_thread {	
	my $thread=shift;
	my $bamfile=shift;
	my $bedfilectgs=shift;
	my $bedfileorfs=shift;
	my $countfile=shift;
	my $mapper=shift;
	my $samtools_soft=shift;
	my $verbose=shift;

	my(%genesincontigs,%genespos);
	open infile,"$bedfileorfs" || die "Can't open file $bedfileorfs\n";
	my($contigid,$initgen,$endgen,$genid);
	while(<infile>) {
		chomp;
		($contigid,$initgen,$endgen,$genid) = split("\t");
        	if(not defined $genesincontigs{$contigid}) {
                	my @arr;
			$genesincontigs{$contigid} = [ @arr ];
                	}
        	push @{$genesincontigs{$contigid}}, $genid;
        	$genespos{$genid} = [$initgen, $endgen];
		}
	close infile;

	my($countreads,%seenreads,$readid);
	my %accum;
	open(my $infile3,"-|", "$samtools_soft view $bamfile -M -L $bedfilectgs") || die "Can't open bam file $bamfile\n";
	while(<$infile3>) { 
		chomp;
		next if(!$_ || ($_=~/^\#/)|| ($_=~/^\@/));
		my @k=split(/\t/,$_);
		my $readid=$k[0];
		if($mapper=~/minimap2/){	#-- Minimap2 can output more than one alignment per read
			if($seenreads{$readid}) { next; }
			else { $seenreads{$readid} = 1 }
		}
		$countreads++;
		my $cigar=$k[5];
		next if($k[2]=~/\*/);
		next if($cigar eq "*");
		# print "\n$_\n";
		my $initread=$k[3];                     
		my $contigid=$k[2];
		my $endread=$initread;
		#-- Calculation of the length match end using CIGAR string

		$cigar=~s/^\d+S//;  #-- Elimination of soft clipping in minimap2 mapping
		$cigar=~s/\d+S$//;  

		while($cigar=~/^(\d+)([IDMNSHPX])/) {
			my $mod=$1;
			my $type=$2;
			if($type=~/M|D|N/) { $endread+=$mod; }	#-- Update end position according to the match found
			elsif($type=~/I|S/) { $endread-=$mod; }
			$cigar=~s/^(\d+)([IDMNSHPX])//g;
			}
		if(not defined $genesincontigs{$contigid}) { next; } # no genes for this contig
		foreach $genid (@{$genesincontigs{$contigid}}) {
			my($initgen, $endgen) = @{$genespos{$genid}};
			my $basesingen;
			next if(($endread<$initgen) || ($initread>$endgen)); # the read doesn't map over this orf
			if((($initread>=$initgen) && ($initread<=$endgen)) && (($endread>=$initgen) && ($endread<=$endgen))) {   #-- Read is fully contained in the gene
				$basesingen=$endread-$initread;
				if($verbose) { print "Read within gene: $readid $initread-$endread $contigid $initgen-$endgen $basesingen\n"; }
				}
			elsif(($initread>=$initgen) && ($initread<=$endgen)) {   #-- Read starts within this gene
				$basesingen=$endgen-$initread;
				if($verbose) {  print "Start of read: $readid $initread-$endread $contigid $initgen-$endgen $basesingen\n"; }
				}
 			elsif(($endread>=$initgen) && ($endread<=$endgen)) {   #-- Read ends within this gene
				$basesingen=$endread-$initgen;
				if($verbose) {  print "End of read: $readid $initread-$endread $contigid $initgen-$endgen $basesingen\n"; }
				}
			elsif(($initread<=$initgen) && ($endread>=$endgen)) {  #-- Gen is fully contained in the read
				$basesingen=$endgen-$initgen;
				if($verbose) {  print "Gene within read: $readid $initread-$endread $contigid $initgen-$endgen $basesingen\n"; }
				}
			else { die "Failed sanity check"; }
			$accum{$genid}{reads}++;
                        $accum{$genid}{bases}+=$basesingen;
			}
		}
	close $infile3;
	print "  $countreads reads counted\n";
	open(outtemp,">$countfile") || die;
	foreach my $h(sort keys %accum) { print outtemp "$h\t$accum{$h}{reads}\t$accum{$h}{bases}\n"; }
	close outtemp;
}


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
	
	my $mapperc=($mappedreads/$totalreadcount)*100;
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
		my $rpkm=($readcount{$rc}{reads}*1000000)/(($longt/1000)*$totalreadcount); #-- Length of contig in Kbs
		my $tpm=$rp{$rc}/$accumrpk;
		if(!$rpkm) { print outfile4 "$rc\t0\t0\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n"; } 
		else { printf outfile4 "$rc\t%.2f\t%.1f\t%.1f\t$longt\t$readcount{$rc}{reads}\t$readcount{$rc}{lon}\t$thissample\n",$coverage,$rpkm,$tpm; }
		}
	close outfile4;	
	return $totalreadcount;
}

