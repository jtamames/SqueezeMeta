#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 27/02/2022 for version 1.6.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs. Uses prinseq to filter out contigs by length (excluding small ones).

#-- This piece of software was developed during JT confinment on board of ship Hesperides, on transit to Antarctica. Plenty of free time.


use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 
use Term::ANSIColor qw(:constants);

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
our(%assemblers);
my $project=$projectname;
do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$installpath,$userdir,$assembler,$outassembly,$mode,$megahit_soft,$mapper,$bowtie2_x_soft,$bowtieref,$bowtie2_build_soft,$mapping_options,$mappingstat, $bwa_soft, $minimap2_soft,$assembler_options,$extassembly,$contigid,$numthreads,$keepsam10,$spades_soft,$flye_soft,$prinseq_soft,$mappingfile,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$contigsfna,$contigslen,$cleaning,$cleaningoptions,$scriptdir,$singletons,$methodsfile,$syslogfile,$norename,$force_overwrite);

my($seqformat,$gzipped,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name,$warnmes,%extassemblies,%datasamples);

my $assemblerdir="$installpath/lib/SqueezeMeta";

if((-e $contigsfna) && (-e $contigslen) && (!$force_overwrite)) { print "  Assembly results already present in file $contigsfna, skipping\n"; exit; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Read samples from samples file

	my %allsamples;
	open(infile1,$mappingfile);  #-- To check for extassemblies in sequential mode
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		$_=~s/\r//g;
		my($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		if($mapreq!~/noassembly/) {
			if($mapreq=~/extassembly\=(.*)/) { $extassemblies{$sample}=$1; }  #-- Store external assemblies if specified in the samples file
			elsif(($mode eq "sequential") && ($sample eq $projectname)) {   #-- If in sequential mode, only assemble the current sample
				$datasamples{$sample}{$iden}{$file}=1; 
				if($iden eq "pair1") { $allsamples{$sample}{"$userdir/$file"}=1; } 
				elsif ($iden eq "pair2") { $allsamples{$sample}{"$userdir/file"}=2; }
				}  
			else { 
				$datasamples{$sample}{$iden}{$file}=1; 
				if($iden eq "pair1") { $allsamples{$sample}{"$userdir/$file"}=1; } 
				elsif ($iden eq "pair2") { $allsamples{$sample}{"$userdir/$file"}=2; }				
				}
			}
		}
	close infile1;
	

#-- If extassembly option was specified, skip all assembly steps

if($extassembly) {
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	print outsyslog "  External assembly provided: $extassembly. Overriding assembly\n";
	if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
	$outassembly=$extassembly; 
	}

#-- Otherwise, assembly all samples

elsif($mode ne "coassembly") {		#-- Sequential, merged or seqmerge: Assembly all samples individually
	my ($seqformat,$p2,$filen1,$filen2)="";
	for my $asamples(sort keys %datasamples) {
		next if(($mode eq "sequential") && ($projectname ne $asamples));
		if($extassemblies{$asamples}) { 
			$outassembly=$extassembly; 
			print "  External assembly provided: $extassembly. Overriding assembly\n";
			print outsyslog "  External assembly provided: $extassembly. Overriding assembly\n";
			if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
			}
		else {	
		
			#-- Concatenate all pair1 and pair2 files
		
			print "  Working for sample $asamples\n";
			my $cat1="cat ";
			my $cat2="cat ";
			foreach my $tp("pair1","pair2") {
			foreach my $par(sort keys %{ $datasamples{$asamples}{$tp} }) { 
				 if($tp eq "pair1") { $cat1.="$userdir/$par "; }
				 if($tp eq "pair2") { $cat2.="$userdir/$par "; $p2=1; }
					if($par=~/\.fastq|\.fq/) { 
						if($seqformat eq "fasta") { die "Mixing formats not allowed in sequence files: fastq/fasta\n"; } 
						$seqformat="fastq"; 
						if($par=~/\.gzip|\.gz/) { $gzipped=1; }
						}
					elsif($par=~/\.fasta|\.fa/) { 
						if($seqformat eq "fastq") { die "Mixing formats not allowed in sequence files: fastq/fasta\n"; } 
						$seqformat="fasta"; 
						if($par=~/\.gzip|\.gz/) { $gzipped=1; }
						}
					else { die "Unrecognized format in file $par: Not fastq or fasta\n"; }			
					}
				}
				
			#--Put the adequate extension (fasta/fastq)	
				
			if($seqformat eq "fastq") { 
				$filen1="par1.fastq";
				$filen2="par2.fastq";
				}
			else {	
				$filen1="par1.fasta";
				$filen2="par2.fasta";
				}
				
			#-- Consider if the files were zipped	
				
			if($gzipped) { $filen1.=".gz"; $filen2.=".gz"; }	
			$cat1.="> $datapath/raw_fastq/$filen1";
			$cat2.="> $datapath/raw_fastq/$filen2";
		
			print "  Preparing files for pair1: $cat1\n";
			print outsyslog "  Preparing files for pair1: $cat1\n";
			system($cat1); 
			if($p2) {			
				print "  Preparing files for pair2: $cat2\n";
				print outsyslog "  Preparing files for pair2: $cat2\n";
				system($cat2);
				}
				
			#-- Call the assemblers	
			
			my $provname="$interdir/01.$asamples.fasta";
			if((-e $provname) && (!$force_overwrite)) { print "  Assembly results already present in file $provname, skipping\n"; }
			else {
				if($mode eq "sequential") { $projectname=$asamples; }
				if($p2) { assembly($projectname,$asamples,$filen1,$filen2); } else { assembly($projectname,$asamples,$filen1); }
				}
					
					#-- Now we need to rename the contigs for minimus2, otherwise there will be contigs with same names in different assemblies

				if($mode ne "sequential") {
					print "  Renaming contigs\n"; 
					system("cp $provname $provname.prov");
					open(outfile1,">$provname") || die "Can't open $contigsfna for writing\n";
					open(infile2,"$provname.prov") || die "Can't open $contigsfna.prov\n";
					while(<infile2>) {
						chomp;
						if($_=~/^\>([^ ]+)/) { 
						my $tname=$1; 
						$_=~s/$tname/$tname\_$asamples/; 
						}
					print outfile1 "$_\n";
					}
				close infile2;
				close outfile1;
				# system("rm $provname.prov");
				}

			}
		last if	($mode eq "sequential"); #-- This avoids iterating for samples that are not the current one
		}
	}
else {      #-- For coassembly: join all samples to assembly them as an unique set

	my($seqformat,$p2,$filen1,$filen2)="";
	my $cat1="cat ";
	my $cat2="cat ";
		
			#-- Concatenate all pair1 and pair2 files
		
	print "  Concatenating all samples: ";
	for my $asamples(sort keys %datasamples) {
		foreach my $tp("pair1","pair2") {
		foreach my $par(sort keys %{ $datasamples{$asamples}{$tp} }) { print " $par ";
			if($tp eq "pair1") { $cat1.="$userdir/$par "; }
			if($tp eq "pair2") { $cat2.="$userdir/$par "; $p2=1; }
			if($par=~/\.fastq|\.fq/) { 
				if($seqformat eq "fasta") { die "Mixing formats not allowed in sequence files: fastq/fasta\n"; } 
				$seqformat="fastq"; 
				if($par=~/\.gzip|\.gz/) { $gzipped=1; }
				}
			elsif($par=~/\.fasta|\.fa/) { 
				if($seqformat eq "fastq") { die "Mixing formats not allowed in sequence files: fastq/fasta\n"; } 
				$seqformat="fasta"; 
				if($par=~/\.gzip|\.gz/) { $gzipped=1; }
				}
			else { die "Unrecognized format in file $par: Not fastq or fasta\n"; }			
			}
		}
	}	
	print "\n";			
	
	#--Put the adequate extension (fasta/fastq)	
				
	if($seqformat eq "fastq") { 
		$filen1="par1.fastq";
		$filen2="par2.fastq";
	}
	else {	
		$filen1="par1.fasta";
		$filen2="par2.fasta";
	}
				
	#-- Call the assemblers	
				
	if($gzipped) { $filen1.=".gz"; $filen2.=".gz"; }	
	$cat1.="> $datapath/raw_fastq/$filen1";
	$cat2.="> $datapath/raw_fastq/$filen2";
		
	print "    pair1: $cat1\n";
	print outsyslog "  Preparing files for pair1: $cat1\n";
	system($cat1); 
	if($p2) {			
		print "    pair2: $cat2\n";
		print outsyslog "  Preparing files for pair2: $cat2\n";
		system($cat2);
		}
				
			#-- Call the assemblers	
				
	if($p2) { assembly($projectname,$projectname,$filen1,$filen2); } else { assembly($projectname,$projectname,$filen1); }
	}
	
#-- If we are in merged or seqmerge modes, do the merging now

if(($mode eq "merged") || ($mode eq "seqmerge")) {
	my $scriptname;
	if($mode eq "merged") { $scriptname="01.merge_assemblies.pl"; }
	elsif($mode eq "seqmerge") { $scriptname="01.merge_sequential.pl"; }
	print outsyslog "\n STEP1 -> $scriptname\n";
	print BLUE " STEP1 -> MERGING ASSEMBLIES: $scriptname\n"; print RESET;
	#if($verbose) { print " (This will take all the individual assemblies and merge them in a single, combined assembly)\n"; }
	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
	if($ecode!=0)   { error_out(1,$scriptname); }
	}

#-- Run prinseq_lite for removing short contigs

if($extassembly) { system("cp $extassembly $contigsfna"); }
else {
	$command="$prinseq_soft -fasta $contigsfna -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna > /dev/null 2>&1";
	print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen \n";
	print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n  ";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	}
	
	
#-- Standardization of contig names

if(!$norename) {
	print "  Renaming contigs\n";
	open(infile1,$contigsfna) || die "Can't open $contigsfna\n";
	my $provcontigs="$tempdir/contigs.prov";
	open(outfile1,">$provcontigs") || die "Can't open $provcontigs for writing\n";
	my $cocount;
	if(!$contigid) { $contigid="$assembler"; }
	while(<infile1>) {
		chomp;
		next if !$_;
		if($_=~/^\>/) {
			$cocount++;
			my $newcontigname="$contigid\_$cocount";
			print outfile1 ">$newcontigname\n";
			}
		else { print outfile1 "$_\n"; }
		}
	close infile1;
	close outfile1;
	system("mv $provcontigs $contigsfna");
	}

#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
print outsyslog "Running prinseq for contig statistics: $command\n  ";
print "Running prinseq for contig statistics: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";

mapping();
	
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

print "  Contigs stored in $contigsfna\n  Number of contigs: $numc\n";
#system("rm $datapath/raw_fastq/par1.$format.gz; rm $datapath/raw_fastq/par2.$format.gz");



my $scriptname;
my $wc=qx(wc -l $contigsfna);
my($wsize,$rest)=split(/\s+/,$wc);
if($singletons) {
	$scriptname="01.remap.pl";
	print outfile3 "1\t$scriptname\n";
	print outsyslog "  STEP1 -> $scriptname\n";
	print BLUE "  STEP1 ->  ADDING SINGLETONS: $scriptname\n"; print RESET;
	# if($verbose) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
	if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname\n"; print RESET; die; }
	my $wc=qx(wc -l $contigsfna);
	my($wsize,$rest)=split(/\s+/,$wc);
	if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname. File $contigsfna is empty!\n"; print RESET; die; }
}
if($wsize<2)         { error_out(1,$scriptname,$contigsfna); }

close outsyslog;
close outmet;


sub mapping {
	print "NOW MAPPING DATA\n";
	my $fastqdir="$datapath/raw_fastq";
	my $samdir="$datapath/sam";
	if(-d $samdir) {} else { system("mkdir $samdir"); }


        #-- Creates Bowtie2 or BWA reference for mapping (index the contigs)

	if($mapper eq "bowtie") {
		print "  Mapping with Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n";
		print outmet "Read mapping against contigs was performed using Bowtie2 (Langmead and Salzberg 2012, Nat Methods 9(4), 357-9)\n"; 
        	if(-e "$bowtieref.1.bt2") { print "  Found reference in $bowtieref.1.bt2, skipping\n"; }
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
        	if(-e "$bowtieref.bwt") { print "Found reference in $bowtieref.bwt, Skipping\n"; }
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

	open(outfile1,">$mappingstat") || die "Can't open mappingstat file $mappingstat for writing\n";	#-- File containing mapping statistics
	print outfile1 "#-- Created by $0, ",scalar localtime,"\n";
	print outfile1 "# Sample\tTotal reads\tMapped reads\tMapping perc\tTotal bases\n";

	my $nums;
	foreach my $thissample(keys %allsamples) {
		next if(($mode eq "sequential") && ($thissample ne $projectname));   #-- If in sequential mode, only assemble the current sample
		my($formatseq,$command,$outsam,$formatoption);
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
		if($keepsam10) { $outsam="$samdir/$projectname.$thissample.sam"; } else { $outsam="$samdir/$projectname.$thissample.current.sam"; }
		if(-e $outsam) { print "  SAM file already found in $outsam, skipping\n"; }
		else {

			#-- Support for single reads
       			if(!$mapper || ($mapper eq "bowtie")) {
           			if($formatseq eq "fasta") { $formatoption="-f"; }
				if($mapping_options eq "") { $mapping_options = "--very-sensitive-local"; } # very-sensitive-local would interfere with custom mapping options so add it here 
    	    			if(-e "$tempdir/$par2name") { $command="$bowtie2_x_soft -x $bowtieref $formatoption -1 $tempdir/$par1name -2 $tempdir/$par2name --quiet -p $numthreads -S $outsam $mapping_options"; }
	    			else { $command="$bowtie2_x_soft -x $bowtieref $formatoption -U $tempdir/$par1name  --quiet -p $numthreads -S $outsam $mapping_options"; } }
        		elsif($mapper eq "bwa") {
            			#Apparently bwa works seamlesly with fasta files as input.
            			if(-e "$tempdir/$par2name") { $command="$bwa_soft mem $bowtieref $tempdir/$par1name $tempdir/$par2name -v 1 -t $numthreads $mapping_options > $outsam"; }
            			else { $command="$bwa_soft mem $bowtieref $tempdir/$par1name -v 1 -t $numthreads $mapping_options > $outsam"; } }
        		elsif($mapper eq "minimap2-ont") {
            			#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            			if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $outsam"; }
            			else { $command="$minimap2_soft -ax map-ont $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $outsam"; } }
        		elsif($mapper eq "minimap2-pb") {
            			#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            			if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $outsam"; }
            			else { $command="$minimap2_soft -ax map-pb $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $outsam"; } }
        		elsif($mapper eq "minimap2-sr") {
            			#Minimap2 does not require to create a reference beforehand, and work seamlessly with fasta as an input.
            			if(-e "$tempdir/$par2name") { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name $tempdir/$par2name -t $numthreads $mapping_options > $outsam"; }
            			else { $command="$minimap2_soft -ax sr $contigsfna $tempdir/$par1name -t $numthreads $mapping_options > $outsam"; } 
				}
			}
                                  
		# print "$command\n";
		print outsyslog "Aligning with $mapper: $command\n";
		system($command);
        	my $ecode = 0;
		if(-e $outsam) {} else { $ecode = system $command; }
        	if($ecode!=0)     { die "An error occurred during mapping!"; }
		
			#-- Counting mapping percentage
			
		my($thisr,$lastr,$mappedreads,$totalreadcount,$totalreadlength);
		my %readcount;
		open(infile4,$outsam) || die "Can't open $outsam\n"; 
		while(<infile4>) {
			chomp;
			my @t=split(/\t/,$_);
			next if($_=~/^\@/);
	
			#-- Use the mapped reads to sum base coverage

			if($t[5]!~/\*/) { 			#-- If the read mapped, accum reads and bases
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
		close infile4;
	
		my $mapperc=($mappedreads/$totalreadcount)*100;
		printf outfile1 "$thissample\t$totalreadcount\t$mappedreads\t%.2f\t$totalreadlength\n",$mapperc;		#-- Mapping statistics
		}
		 
		if($warnmes) { 
			print outfile1 "\n# Notice that mapping percentage is low (<50%) for some samples. This is a potential problem,  meaning that most reads are not represented in the assembly\n";
			print RED "\n# Notice that mapping percentage is low (<50%) for some samples. This is a potential problem,  meaning that most reads are not represented in the assembly\n";
			if($mincontiglen>200) { 
			print outfile1 "# Notice also that you set the minimum contig length to $mincontiglen. In this way you are removing the contigs shorter than that size. This can be, at least partially, the cause of this low mapping percentage\n";
			print outfile1 "# It is likely that redoing the analysis with the default minimum contig length (200) will improve the results\n";
			print outfile1 "# If not, you could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n";
			print RED "# Notice also that you set the minimum contig length to $mincontiglen. In this way you are removing the contigs shorter than that size. This can be, at least partially, the cause of this low mapping percentage\n";
			print RED "# It is likely that redoing the analysis with the default minimum contig length (200) will improve the results\n";
			print RED "# If not, you could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n";
			}
		else { 
			print outfile1 "# You could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n"; 
			print RED "# You could redo your analysis using assignment of the raw reads instead of relying on the assembly. Use sqm_reads.pl or sqm_longreads.pl for this purpose. That strategy will not provide bins\n"; 
			}
		print RESET;	
		}

	
	}

	

sub assembly {	

	my($prname,$samplename,$p1name,$p2name)=@_;
	my $par1name="$datapath/raw_fastq/$p1name";
	my $par2name="$datapath/raw_fastq/$p2name";

	#-- Checks the assembler and call to the appropriate script (as specified in SqueezeMeta_conf.pl)

	my $command="perl $assemblerdir/$assemblers{$assembler} $prname $samplename $par1name $par2name";
	print "  Running assembly with $assembler: $command\n";
	print outsyslog "Running assembly with $assembler: $command\n";
	system($command);
	
	
	}


