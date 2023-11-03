#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 27/02/2022 for version 1.6.0, (c) Javier Tamames, CNB-CSIC
#-- Runs assembly programs. Uses prinseq to filter out contigs by length (excluding small ones).

#-- This piece of software was developed during JT confinment on board of ship Hesperides, on transit to Antarctica. Plenty of free time.


use strict;
use Cwd;
use Linux::MemInfo;
use lib "."; 
use Term::ANSIColor qw(:constants);
use File::Basename;

$|=1;

my $pwd=cwd();
my $projectdir=$ARGV[0];
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
our(%assemblers);
my $project=$projectname;

#-- Configuration variables from conf file

our($datapath,$installpath,$userdir,$assembler,$outassembly,$mode,$megahit_soft,$assembler_options,$extassembly,$extbins,$contigid,$numthreads,$spades_soft,$flye_soft,$prinseq_soft,$mappingfile,$trimmomatic_soft,$canu_soft,$canumem,$mincontiglen,$resultpath,$interdir,$tempdir,$binresultsdir, $contigsfna,$contigslen,$cleaning,$cleaningoptions,$scriptdir,$singletons,$methodsfile,$syslogfile,$norename,$force_overwrite);

my($seqformat,$gzipped,$outassemby,$trimmomatic_command,$command,$thisname,$contigname,$seq,$len,$par1name,$par2name,%extassemblies,%datasamples);

my $assemblerdir="$installpath/lib/SqueezeMeta";

if((-e $contigsfna) && (-e $contigslen) && (!$force_overwrite)) { print "  Assembly results already present in file $contigsfna, skipping\n"; exit; }

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

#-- Read samples from samples file

	open(infile1,$mappingfile);  #-- To check for extassemblies in sequential mode
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		$_=~s/\r//g;
		my($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		if($mapreq!~/noassembly/) {
			if($mapreq=~/extassembly\=(.*)/) { $extassemblies{$sample}=$1; }  #-- Store external assemblies if specified in the samples file
			elsif(($mode eq "sequential") && ($sample eq $projectname)) { $datasamples{$sample}{$iden}{$file}=1; }  #-- If in sequential mode, only assemble the current sample
			else { $datasamples{$sample}{$iden}{$file}=1; }
			}
		}
	close infile1;
	

#-- If extassembly option was specified, skip all assembly steps

if($extassembly) {
	print "  External assembly provided: $extassembly. Overriding assembly\n";
	print outsyslog "  External assembly provided: $extassembly. Overriding assembly\n";
	if(-e $extassembly) {} else { die "Can't find assembly file $extassembly\n"; }
	}

#-- If extbins option was specified, skip all assembly steps
elsif($extbins) {
        print "  Directory with external bins provided: $extbins. Overriding assembly\n";
        print outsyslog "  Directory with external bins provided: $extbins. Overriding assembly\n";
        if(-d $extbins) {} else { die "Can't find bin directory $extbins\n"; }
        }


#-- Otherwise, assembly all samples

elsif($mode ne "coassembly") {		#-- Sequential, merged, seqmerge or clustered: Assembly all samples individually
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
		
			$p2=0;
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
				if($p2) { assembly($projectdir,$asamples,$filen1,$filen2); } else { assembly($projectdir,$asamples,$filen1); }
				}
					
					#-- Now we need to rename the contigs for minimus2, otherwise there will be contigs with same names in different assemblies

				if($mode ne "sequential") {
					print "  Renaming contigs in $provname\n"; 
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
				
	if($p2) { assembly($projectdir,$projectname,$filen1,$filen2); } else { assembly($projectdir,$projectname,$filen1); }
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
	if($ecode!=0)   { 
		 print RED; print "Stopping in STEP1 -> $scriptname. Program finished abnormally\n"; print RESET; print outsyslog "Stopping in STEP1 -> $scriptname. Program finished abnormally\n"; die; 
		 }
	}

#-- Run prinseq_lite for removing short contigs

if($mode eq "clustered") { system("cat $interdir/01*fasta > $contigsfna"); }

if($extassembly) { system("cp $extassembly $contigsfna"); }
elsif($extbins) {}
elsif(-e $contigsfna) {
	$command="$prinseq_soft -fasta $contigsfna -min_len $mincontiglen -out_good $resultpath/prinseq; mv $resultpath/prinseq.fasta $contigsfna > /dev/null 2>&1";
	print "  Running prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4) for selecting contigs longer than $mincontiglen \n";
	print outsyslog "Running prinseq for selecting contigs longer than $mincontiglen: $command\n  ";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }
	print outmet "Short contigs (<$mincontiglen bps) were removed using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	}
else { die "Assembly not present in $contigsfna, exiting\n"; }	
	
#-- Standardization of contig names (and creating assembly file and bin files if using $extbins)

if(!$norename) { print "  Renaming contigs in $contigsfna\n"; }
my $provcontigs="$tempdir/contigs.prov";
open(outfile1,">$provcontigs") || die "Can't open $provcontigs for writing\n";

# What are our input files? Only $contigsfna, unless we are running in -extbins mode

my @fastafiles;
if(!$extbins) { @fastafiles = ($contigsfna); }
else {
	@fastafiles = glob("$extbins/*");
	if (scalar @fastafiles < 1) { die "No files were found in $extbins\n"; }
}

# Now process them, changing the contig names if needed
my $cocount;
my ($binfile, $path, $suffix);
my $assembler2 = $assembler;
if($assembler2 eq "spades-base") { $assembler2 = "spades"; }
for(@fastafiles) {
	open(infile1,$_) || die "Can't open $_\n";
	if($extbins) {
		($binfile,$path,$suffix) = fileparse($_, qr/\.[^.]*/);
		$binfile = "$binfile.fa"; # Ensure fa extension for compatibility with step 17
		$binfile = "$binresultsdir/$binfile";
		open(outfileB,">$binfile") || die "Can't open $binfile for writing\n";
		}
	while(<infile1>) {
		chomp;
		next if !$_;
		if($_=~/^\>([^ ]+)/) {
			my @fd=split(/\_/,$1); 
			$cocount++;
			my $newcontigname;
			if($norename)                 { $newcontigname = $_;                      }
			elsif($contigid)              { $newcontigname=">$contigid\_$cocount";    }
			#elsif($extassembly||$extbins) { $newcontigname=">contig\_$cocount";       }
			#elsif($mode eq "clustered")   { $newcontigname=">$fd[$#fd]\_$cocount";    }
			#elsif($mode eq "merged")      { $newcontigname=">merged\_$cocount";       }
			#elsif($mode eq "seqmerge")    { $newcontigname=">seqmerge\_$cocount";     }
			#else                          { $newcontigname=">$assembler2\_$cocount";  }
			else                          { $newcontigname=">$projectname\_$cocount"; }
			print outfile1 "$newcontigname\n";
			if($extbins) { print outfileB "$newcontigname\n"; }
			}
		else {
			print outfile1 "$_\n"; 
			if($extbins) { print outfileB "$_\n"; }
			}
		}
	close infile1;
	if($extbins) { close outfileB; }
	}
close outfile1;
system("mv $provcontigs $contigsfna");


#-- Run prinseq_lite for statistics

$command="$prinseq_soft -fasta $contigsfna -stats_len -stats_info -stats_assembly > $interdir/01.$project.stats";
print outsyslog "Running prinseq for contig statistics: $command\n  ";
print "Running prinseq for contig statistics: $command\n  ";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "Contig statistics were done using prinseq (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";

	
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

#if($singletons) {			#-- REMOVING THIS SINCE SINGLETONS ARE CALLED IN MAIN SQUEEZEMETA.PL SCRIPT!
#	$scriptname="01.remap.pl";
#	print outfile3 "1\t$scriptname\n";
#	print outsyslog "  STEP1 -> $scriptname\n";
#	print BLUE "  STEP1 ->  ADDING SINGLETONS: $scriptname\n"; print RESET;
#	# if($verbose) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
#	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
#	if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname\n"; print RESET; die; }
#	my $wc=qx(wc -l $contigsfna);
#	my($wsize,$rest)=split(/\s+/,$wc);
#	if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname. File $contigsfna is empty!\n"; print RESET; die; }
#}

if($wsize<2)      { 
		 print RED; print "Stopping in STEP1 -> $scriptname. Program finished abnormally\n"; print RESET; print outsyslog "Stopping in STEP1 -> $scriptname. Program finished abnormally\n"; die; 
		 }

close outsyslog;
close outmet;


sub assembly {	

	my($prname,$samplename,$p1name,$p2name)=@_;
	my $par1name="$datapath/raw_fastq/$p1name";
	my $par2name="";
	if($p2name) { $par2name="$datapath/raw_fastq/$p2name"; } 

	#-- Checks the assembler and call to the appropriate script (as specified in SqueezeMeta_conf.pl)

	my $command="perl $assemblerdir/$assemblers{$assembler} $prname $samplename $par1name $par2name";
	print "  Running assembly with $assembler: $command\n";
	print outsyslog "Running assembly with $assembler: $command\n";
	system($command);
	
	
	}


