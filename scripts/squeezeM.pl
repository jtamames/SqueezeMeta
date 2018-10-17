#!/usr/bin/perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use strict;

my $version="0.3.0, Ago 2018";
my $start_run = time();

###scriptdir patch, Fernando Puente-Sánchez, 29-V-2018
use File::Basename;
our $scriptdir = dirname(__FILE__);
our $installpath = "$scriptdir/..";
###

our $pwd=cwd();
our($nocog,$nokegg,$nopfam,$nobins,$nomaxbin,$nometabat)="0";
our($numsamples,$numthreads,$mode,$mincontiglen,$assembler,$mapper,$counter,$project,$equivfile,$rawfastq,$evalue,$miniden,$spadesoptions,$megahitoptions,$assembler_options,$ver,$hel);
our($databasepath,$extdatapath,$softdir,$basedir,$datapath,$resultpath,$tempdir,$mappingfile,$contigsfna,$contigslen,$mcountfile,$rnafile,$gff_file,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$rpkmfile,$coveragefile,$contigcov,$contigtable,$mergedfile,$bintax,$checkmfile,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bwa_soft,$minimap2_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft);
our %bindirs;  

#-- Handle variables from command line

my $result = GetOptions ("t=i" => \$numthreads,
                     "m|mode=s" => \$mode,
                     "c|contiglen=i" => \$mincontiglen,
                     "a=s" => \$assembler,
                     "map=s" => \$mapper,
#                     "count=s" => \$counter,
                     "p=s" => \$project,
                     "s|samples=s" => \$equivfile,
                     "f|seq=s" => \$rawfastq, 
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nopfam" => \$nopfam,   
		     "nobins" => \$nobins,   
		     "nomaxbin" => \$nomaxbin,   
		     "nometabat" => \$nometabat,   
		     "e|evalue=f" => \$evalue,   
		     "minidentity=f" => \$miniden,   
		     "spades_options=s" => \$spadesoptions,
		     "megahit_options=s" => \$megahitoptions,
		     "v" => \$ver,
		     "h" => \$hel
		    );

#-- Set some default values

if(!$numthreads) { $numthreads=12; }
if(!$mincontiglen) { $mincontiglen=1200; }
if(!$assembler) { $assembler="megahit"; }
if(!$mapper) { $mapper="bowtie"; }
if(!$counter) { $counter="bedtools"; }
if(!$evalue) { $evalue=1e-03; }
if(!$miniden) { $miniden=50; }
if(!$nocog) { $nocog=0; }
if(!$nokegg) { $nokegg=0; }
if(!$nopfam) { $nopfam=0; }
if(!$nobins) { $nobins=0; }
if(!$nomaxbin) { $nomaxbin=0; }
if(!$nometabat) { $nometabat=0; }

#-- Check if we have all the needed options

my $helptext="Usage: squeezeM.pl -m <mode> -p <projectname> -s <equivfile> -f <raw fastq dir> <options>\n\nArguments:\n\n Mandatory parameters:\n  -m: Mode (sequential, coassembly, merged) (REQUIRED)\n  -s|-samples: Samples file (REQUIRED)\n  -f|-seq: Fastq read files' directory (REQUIRED)\n  -p: Project name (REQUIRED in coassembly and merged modes)\n\n Assembly:\n  -a: assembler [megahit,spades] (Default: $assembler)\n  --megahit_options: Options for megahit assembler\n  --spades_options: Options for spades assembler\n  -map: mapping software [bowtie,bwa,minimap2-ont,minimap2-pb,minimap2-sr]\n  -c|-contiglen: Minimum length of contigs (Default:$mincontiglen)\n\n Functional & taxonomic assignments:\n  --nocog: Skip COG assignment (Default: no)\n  --nokegg: Skip KEGG assignment (Default: no)\n  --nopfam: Skip Pfam assignment  (Default: no)\n  -e|-evalue: max evalue for discarding hits diamond run  (Default: 1e-03)\n  -miniden: identity perc for discarding hits in diamond run  (Default: 50)\n\n Binning:\n  --nobins: Skip all binning  (Default: no)\n  --nomaxbin: Skip MaxBin binning  (Default: no)\n  --nometabat: Skip MetaBat2 binning  (Default: no)\n\n Performance:\n  -t: Number of threads (Default:$numthreads)\n\n Information\n  -v: Version number\n  -h: help\n";

print "\nSqueezeM v$version - (c) J. Tamames, CNB-CSIC\n\n";

if($ver) { exit; }
if($hel) { die "$helptext\n"; } 
if((!$rawfastq) || (!$equivfile) || (!$mode)) { die "$helptext\n"; }
if(($mode!~/sequential/i) && (!$project)) { die "$helptext\n"; }
if(($mode=~/sequential/i) && ($project)) { die "$helptext\nPlease DO NOT specify project name in sequential mode. The name will be read from the samples in $equivfile\n"; }
if($mode!~/sequential|coassembly|merged/i) { die "$helptext\n"; }
if($mapper!~/bowtie|bwa|minimap2-ont|minimap2-pb|minimap2-sr/i) { die "$helptext\n"; }
if($rawfastq=~/^\//) {} else { $rawfastq="$pwd/$rawfastq"; }

my $currtime=timediff();
my $commandline = join " ", $0, @ARGV;
print "Run started ",scalar localtime," in $mode mode\n";


#--------------------------------------------------------------------------------------------------
#----------------------------------- SEQUENTIAL MODE ----------------------------------------------

if($mode=~/sequential/i) { 
        
	my(%allsamples,%ident,%noassembly);
	my($sample,$file,$iden,$mapreq);
	tie %allsamples,"Tie::IxHash";

	#-- Reading the sample file given by the -s option, to locate the sample files

	print "Now reading samples\n";
	open(infile1,$equivfile) || die;
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		if((!$sample) || (!$file) || (!$iden)) { die "Bad format in samples file $equivfile\n"; }
		$allsamples{$sample}{$file}=1;
		$ident{$sample}{$file}=$iden;
		if($mapreq && ($mapreq=~/noassembly/i)) { $noassembly{$file}=1; }   #-- Files marked for no assembly (but they will be mapped)
	}
	close infile1;

	my @nmg=keys %allsamples;
	$numsamples=$#nmg+1;
	print "$numsamples metagenomes found";
	print "\n";

	open(outfile1,">$pwd/global_progress") || die;  	#-- An index indicating where are we and which parts of the method finished already. For the global process
	open(outfile2,">$pwd/global_syslog") || die; 		 #-- A log file for the global proccess
	print outfile2 "Run started ",scalar localtime," in SEQUENTIAL mode (it will proccess all metagenomes sequentially)\n";
	$commandline = join " ", $0, @ARGV;
	print outfile2 "Command: $commandline\n"; 
	print outfile2 "Options: threads=$numthreads; contiglen=$mincontiglen; assembler=$assembler; sample file=$equivfile; raw fastq=$rawfastq\n";

	#-- Now we start processing each sample individually

	foreach my $thissample(keys %allsamples) {

		#-- We start creating directories, progress and log files
         
		$project=$thissample;
		my $projectdir="$pwd/$thissample";
		if (-d $projectdir) { die "Project name $projectdir already exists\n"; } else { system("mkdir $projectdir"); }
		print "Working with $thissample\n";
		print outfile1 ">$thissample\n";
	
		open(outfile3,">$projectdir/progress") || die;  	#-- An index indicating where are we and which parts of the method finished already. For the global process
		open(outfile4,">$projectdir/syslog") || die;  	#-- A log file for the global proccess
		$currtime=timediff();
		print outfile4 "Run started ",scalar localtime," in SEQUENTIAL mode (it will proccess all metagenomes sequentially)\n";
		print "Run started ",scalar localtime," in SEQUENTIAL mode\n";
                my $params = join(" ", @ARGV);
                print outfile2 "$0 $params\n";
		print outfile2 "Run started for $thissample, ",scalar localtime,"\n";
		print outfile4 "Project: $project\n";
		print outfile4 "Map file: $equivfile\n";
		print outfile4 "Fastq directory: $rawfastq\n";
                print outfile4 "[",$currtime->pretty,"]: STEP0 -> squeezeM.pl\n";
                print outfile2 "[",$currtime->pretty,"]: STEP0 -> squeezeM.pl\n";
		print "Now creating directories\n";
		open(infile2,"$scriptdir/squeezeM_conf.pl") || die;
	
		#-- Creation of the new configuration file for this sample
	
		open(outfile5,">$projectdir/squeezeM_conf.pl") || die;

		print outfile5 "\$mode=\"$mode\";\n\n";
                print outfile5 "\$installpath=\"$installpath\";\n";

		while(<infile2>) {
                        chomp;
                        next if !$_;
			if($_=~/^\$basedir/) { print outfile5 "\$basedir=\"$pwd\";\n"; }
			elsif($_=~/^\$projectname/) { print outfile5 "\$projectname=\"$project\";\n"; }
			elsif($_=~/^\$evalue/) { print outfile5 "\$evalue=$evalue;\n"; }
			elsif($_=~/^\$miniden/) { print outfile5 "\$miniden=$miniden;\n"; }
			elsif($_=~/^\$nocog/) { print outfile5 "\$nocog=$nocog;\n"; }
			elsif($_=~/^\$nokegg/) { print outfile5 "\$nokegg=$nokegg;\n"; }
			elsif($_=~/^\$nopfam/) { print outfile5 "\$nopfam=$nopfam;\n"; }
			elsif($_=~/^\$nobins/) { print outfile5 "\$nobins=$nobins;\n"; }
			elsif($_=~/^\$nomaxbin/) { print outfile5 "\$nomaxbin=$nomaxbin;\n"; }
			elsif($_=~/^\$nometabat/) { print outfile5 "\$nometabat=$nometabat;\n"; }
                        elsif($_=~/^\$mapper/) { print outfile5 "\$mapper=\"$mapper\";\n"; }
			else { print outfile5 "$_\n"; }
        	}
	 	close infile2; 

                if($counter=~/featurecounts/i) { $counter="featureCounts"; }
                if($assembler eq "megahit") { $assembler_options=$megahitoptions; } else { $assembler_options=$spadesoptions; }
                print outfile5 "\n#-- Options\n\n\$numthreads=$numthreads;\n\$mincontiglen=$mincontiglen;\n\$assembler=\"$assembler\";\n";
                if($assembler_options) { print outfile5 "\$assembler_options=$assembler_options"; }
                close outfile5;
        
		#-- Creation of directories
	    
		print "Reading configuration from $projectdir/squeezeM_conf.pl\n";
		do "$projectdir/squeezeM_conf.pl";
		system ("mkdir $datapath");
 		system ("mkdir $resultpath");
 		system ("mkdir $tempdir");
 		system ("mkdir $datapath/raw_fastq"); 
	
		#-- Linkage of files to put them into our data directories
	
		print "Now linking read files\n";
 		foreach my $file(sort keys %{ $allsamples{$thissample} }) {
 			 if(-e "$rawfastq/$file") { 
  				my $tufile="$datapath/raw_fastq/$file";
 				# system("cp $rawfastq/$file $tufile");
 				 system("ln -s $rawfastq/$file $tufile");
 				# if($tufile!~/\.gz$/) { system("gzip $datapath/raw_fastq/$file"); }
			}
 			 else { die "Cannot find read file $file (Sample $sample)\n"; }

	system("cp $equivfile $mappingfile");
			#open(out2,">$mappingfile") || die;     #-- If we gzipped the files, we have to change the mapping file. This is no loger needed with linkage
 			#open(in2,$equivfile) || die;
 			#while(<in2>) {
  				#chomp;
 				#next if(!$_ || ($_=~/^\#/));
  				#@u=split(/\t/,$_);
  				#next if($u[0] ne $thissample);   #-- The new mapping file contains only the files for the current sample
  				#if($u[1]!~/gz$/) { $u[1].=".gz"; }
  				#$line=join("\t",@u);
  				#print out2 "$line\n";
			#}
			#close in2;
			#close out2;
	
		}

		#-- Preparing the files for the assembly, merging all files corresponding to each pair

 		my($par1files,$par2files)=0;
		my($par1name,$par2name);
 		print "Now preparing files\n";
 		my($ca1,$ca2)="";
 		foreach my $afiles(sort keys %{ $ident{$thissample} }) {
 			next if($noassembly{$afiles});
  			my $gzfiles=$afiles;
  			#if($gzfiles!~/gz$/) { $gzfiles.=".gz"; }
 			if($ident{$thissample}{$afiles} eq "pair1") { $ca1.="$datapath/raw_fastq/$gzfiles "; $par1files++; } 
			else { $ca2.="$datapath/raw_fastq/$gzfiles "; $par2files++; } 
			if($ident{$thissample}{$afiles}=~/gz$/) { $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }
			else { $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
		}
		if($par1files>1) { 
			my $command="cat $ca1 > $par1name"; 
			print "$command\n"; 
			system($command); 
		} 
		else {
                        #my $command="cp $ca1 $par1name"; 
                        my $command="ln -s $ca1 $par1name";
                        print "$command\n";
                        system($command);
 
		}
 		if($par2files>1) { system("cat $ca2 > $par2name"); } 
                elsif ($par2files==1) { system("ln -s $ca2 $par2name"); }    #-- Support for single reads
                #else { system("cp $ca2 $par2name"); }
		#-- CALL TO THE STANDARD PIPELINE
		
		pipeline();
		
		
 		close outfile4;		#-- Closing log file for the sample
 		close outfile3;		#-- Closing progress file for the sample
							       
	}    #-- End of this sample, go for the next

}            #-- END


#-----------------------------------------------------------------------------------------------------
#----------------------------------- COASSEMBLY AND MERGED MODES -------------------------------------

else {      

	my $projectdir="$pwd/$project";
	if (-d $projectdir) { die "Project name $projectdir already exists\n"; } else { system("mkdir $projectdir"); }
		
	#-- We start creating directories, progress and log files
                        
	open(outfile3,">$pwd/$project/progress") || die;  #-- Un indice que indica en que punto estamos (que procedimientos han terminado)
	open(outfile4,">$pwd/$project/syslog") || die;
        my $params = join(" ", @ARGV);
        print outfile4 "$0 $params\n";
	print outfile4 "Run started ",scalar localtime," in $mode mode\n";
	print outfile4 "Command: $commandline\n"; 
	print outfile4 "Project: $project\n";
	print outfile4 "Map file: $equivfile\n";
	print outfile4 "Fastq directory: $rawfastq\n";
	print outfile4 "Options: threads=$numthreads; contiglen=$mincontiglen; assembler=$assembler;\n";
	print outfile4 "[",$currtime->pretty,"]: STEP0 -> squeezeM.pl\n";
     
	print "Now creating directories\n";
	
	#-- Creation of the new configuration file for this sample
	open(infile3,"$scriptdir/squeezeM_conf.pl") || die "Cannot open $scriptdir/squeezeM_conf.pl\n";
	open(outfile6,">$projectdir/squeezeM_conf.pl") || die;

	print outfile6 "\$mode=\"$mode\";\n\n";
        print outfile6 "\$installpath=\"$installpath\";\n";
	while(<infile3>) {
		if($_=~/^\$basedir/) { print outfile6 "\$basedir=\"$pwd\";\n"; }
		elsif($_=~/^\$projectname/) { print outfile6 "\$projectname=\"$project\";\n"; }
		elsif($_=~/^\$evalue/) { print outfile6 "\$evalue=$evalue;\n"; }
		elsif($_=~/^\$miniden/) { print outfile6 "\$miniden=$miniden;\n"; }
		elsif($_=~/^\$nocog/) { print outfile6 "\$nocog=$nocog;\n"; }
		elsif($_=~/^\$nokegg/) { print outfile6 "\$nokegg=$nokegg;\n"; }
		elsif($_=~/^\$nopfam/) { print outfile6 "\$nopfam=$nopfam;\n"; }
		elsif($_=~/^\$nobins/) { print outfile6 "\$nobins=$nobins;\n"; }
		elsif($_=~/^\$nomaxbin/) { print outfile6 "\$nomaxbin=$nomaxbin;\n"; }
		elsif($_=~/^\$nometabat/) { print outfile6 "\$nometabat=$nometabat;\n"; }
                elsif($_=~/^\$mapper/) { print outfile6 "\$mapper=\"$mapper\";\n"; }
		elsif(($_=~/^\%bindirs/) && ($nomaxbin)) { print outfile6 "\%bindirs=(\"metabat2\",\"\$resultpath/metabat2\");\n"; }
		elsif(($_=~/^\%bindirs/) && ($nometabat)) { print outfile6 "\%bindirs=(\"maxbin\",\"\$resultpath/maxbin\");\n"; }
		else { print outfile6 $_; }
	 }
	close infile3;

	if($assembler eq "megahit") { $assembler_options=$megahitoptions; } else { $assembler_options=$spadesoptions; }
	print outfile6 "\n#-- Options\n\n\$numthreads=$numthreads;\n\$mincontiglen=$mincontiglen;\n\$assembler=$assembler;\n";
	if($assembler_options) { print outfile6 "\$assembler_options=$assembler_options"; }
	close outfile6;

	print "Reading configuration from $projectdir/squeezeM_conf.pl\n";
	do "$projectdir/squeezeM_conf.pl" || die;

	
	#-- Creation of directories
	
	system ("mkdir $datapath");
	system ("mkdir $resultpath");
	system ("mkdir $tempdir");
	system ("mkdir $datapath/raw_fastq"); 
 
	#-- Preparing the files for the assembly
	   
	moving();
	
	#-- CALL TO THE STANDARD PIPELINE
	
	pipeline();
	
	close outfile4;  #-- Closing log file for the sample
	close outfile3;	  #-- Closing progress file for the sample

}                        #------ END


#----------------------------PREPARING FILES FOR MERGING AND COASSEMBLY MODES-------------------------------------------

sub moving {

	print "Now moving read files\n";
	
	#-- Reading samples from the file specified with -s option
	
	my(%allsamples,%ident,%noassembly);
	open(infile4,$equivfile) || die;
	while(<infile4>) {
 		chomp;
 		next if(!$_ || ($_=~/^\#/));
		my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		$allsamples{$sample}=1;
		$ident{$file}=$iden;
		if(($mapreq) && ($mapreq=~/noassembly/i)) { $noassembly{$file}=1; }    #-- Files flagged for no assembly (but they will be mapped)
		if(-e "$rawfastq/$file") { 
			my $tufile="$datapath/raw_fastq/$file";
			# system("cp $rawfastq/$file $tufile");
			system("ln -s $rawfastq/$file $tufile");
			# if($tufile!~/\.gz$/) { system("gzip $datapath/raw_fastq/$file"); }
		}
		else { die "Cannot find read file $file (Sample $sample)\n"; }
	}
	close infile4;

	#-- Setting number of samples for running binning or not

	my @nmg=keys %allsamples;
	$numsamples=$#nmg+1;
	print outfile3 "Samples:$numsamples\nMode:$mode\n0\n";
	if($numsamples==1) { print "$numsamples sample found: Skipping all binning methods\n"; }
	else { print "$numsamples samples found\n"; }

	system("cp $equivfile $mappingfile");
	#open(out2,">$mappingfile") || die;     #-- If we gzipped the files, we have to change the mapping file. No longer needed
	#open(in2,$equivfile) || die;
	#while(<in2>) {
	# chomp;
	# next if(!$_ || ($_=~/^\#/));
 	#@u=split(/\t/,$_);
	# if($u[1]!~/gz$/) { $u[1].=".gz"; }
	# $line=join("\t",@u);
	# print out2 "$line\n";
	#             }
	#close in2;
	#close out2;

	#-- For coassembly mode, we merge all individual files for each pair

	if($mode=~/coassembly/) {
		print "Now merging files\n";
		my($par1name,$par2name,$par1files,$par2files,$ca1,$ca2);
		foreach my $afiles(sort keys %ident) {
			next if($noassembly{$afiles});
			my $gzfiles=$afiles;
			# if($gzfiles!~/gz$/) { $gzfiles.=".gz"; }
			if($gzfiles=~/gz$/) { $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }
			else { $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
			if($ident{$afiles} eq "pair1") { $ca1.="$datapath/raw_fastq/$gzfiles "; $par1files++; } else { $ca2.="$datapath/raw_fastq/$gzfiles "; $par2files++; } 
		}
				     
		if($par1files>1) { system("cat $ca1 > $par1name"); } else { system("ln -s $ca1 $par1name"); }
		if($par2files>1) { system("cat $ca2 > $par2name"); } elsif($par2files==1) { system("ln -s $ca2 $par2name"); }  #-- Support for single reads
	}
}               #-- END


#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}

#---------------------------------------- PIPELINE --------------------------------

sub pipeline {

	if(-e "$tempdir/$project.log") { system("rm $tempdir/$project.log"); }
	my $rpoint=0;

    #-------------------------------- STEP1: Run assembly

		#-- In coassembly mode

	if($mode=~/coassembly/) {
		my $scriptname="01.run_assembly.pl";
		print outfile3 "1\t$scriptname ($assembler)\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
		print "[",$currtime->pretty,"]: STEP1 -> RUNNING CO-ASSEMBLY: $scriptname ($assembler)\n";
		system("perl $scriptdir/$scriptname $project ");
		if(-s $contigsfna<1000) { die "Stopping in STEP1 -> $scriptname ($assembler)\n"; }
	}

		#-- In merged mode. Includes merging assemblies

	elsif($mode=~/merged/) {
		my $scriptname="01.run_assembly_merged.pl";
		print outfile3 "1\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
	
			#-- Merging individual assemblies
 
		my $scriptname="01.merge_assemblies.pl";
		print outfile3 "1.5\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1.5 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP1.5 -> MERGING ASSEMBLIES: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $contigsfna<1000) { die "Stopping in STEP1.5 -> $scriptname\n"; }
	}
	
		#-- In sequential mode. 
     
	elsif($mode=~/sequential/) {
		my $scriptname="01.run_assembly.pl";
 		print outfile3 "1\t$scriptname\n";
 		$currtime=timediff();
 		print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname\n";
 		print "[",$currtime->pretty,"]: STEP1 ->  RUNNING ASSEMBLY: $scriptname\n";
 		system("perl $scriptdir/$scriptname $project");
		if(-s $contigsfna<1000) { die "Stopping in STEP1 -> $scriptname\n"; }
	}		
			
    #-------------------------------- STEP2: Run RNA prediction

	if($rpoint<=2) {
		my $scriptname="02.run_barrnap.pl";
		print outfile3 "2\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP2 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP2 -> RNA PREDICTION: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		my $masked="$resultpath/02.$project.maskedrna.fasta";
		if(-s $masked<1000) { die "Stopping in STEP2 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if($rpoint<=3) {
		my $scriptname="03.run_prodigal.pl";
 		print outfile3 "3\t$scriptname\n";
 		$currtime=timediff();
 		print outfile4 "[",$currtime->pretty,"]: STEP3 -> $scriptname\n";
 		print "[",$currtime->pretty,"]: STEP3 -> ORF PREDICTION: $scriptname\n";
 		system("perl $scriptdir/$scriptname $project");
		if(-s $aafile<1000) { die "Stopping in STEP3 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if($rpoint<=4) {
		my $scriptname="04.rundiamond.pl";
		print outfile3 "4\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP4 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP4 -> HOMOLOGY SEARCHES: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $taxdiamond<1000) { die "Stopping in STEP4 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if($rpoint<=5) {
		if(!$nopfam) {
			my $scriptname="05.run_hmmer.pl";
			print outfile3 "5\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP5 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP5 -> HMMER/PFAM: $scriptname\n";
			system("perl $scriptdir/$scriptname $project");
			if(-s $pfamhmmer<1000) { die "Stopping in STEP5 -> $scriptname\n"; }
		}
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if($rpoint<=6) {
		my $scriptname="06.lca.pl";
		print outfile3 "6\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP6 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP6 -> TAXONOMIC ASSIGNMENT: $scriptname\n";
		my $ecode = system("perl $scriptdir/$scriptname $project");
                if($ecode!=0)        { die "Stopping in STEP6 -> $scriptname\n"; }
		if(-s $fun3tax<1000) { die "Stopping in STEP6 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if($rpoint<=7) {
		my $scriptname="07.fun3assign.pl";
		if((!$nocog) || (!$nokegg) || (!$nopfam)) {
		print outfile3 "7\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP7 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP7 -> FUNCTIONAL ASSIGNMENT: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if((-s $fun3cog<1000) && (-s $fun3kegg<1000) && (-s $fun3pfam<1000)) { die "Stopping in STEP7 -> $scriptname\n"; }
		}
	}
			
    #-------------------------------- STEP8: Taxonomic annotation for the contigs (consensus of gene annotations)

	if($rpoint<=8) {
		my $scriptname="08.summarycontigs3.pl";
		print outfile3 "8\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP8 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP8 -> CONTIG TAX ASSIGNMENT: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $alllog<1000) { die "Stopping in STEP8 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP9: Mapping of reads onto contigs for abundance calculations
	
	if($rpoint<=9) {
		my $scriptname="09.mapbamsamples.pl";
		print outfile3 "9\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP9 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP9 -> MAPPING READS: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $rpkmfile<1000) { die "Stopping in STEP9 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP10: Count of taxa abundances
	
	if($rpoint<=10) {
		my $scriptname="10.mcount.pl";
		print outfile3 "10\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP10 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP10 -> COUNTING TAX ABUNDANCES: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $mcountfile<1000) { die "Stopping in STEP10 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP11: Count of function abundances
	
	if($rpoint<=11) {
		my $scriptname="11.funcover.pl";
		print outfile3 "11\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP11 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP11 -> COUNTING FUNCTION ABUNDANCES: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		my $cogfuncover="$resultpath/11.$project.cog.funcover";
		my $keggfuncover="$resultpath/11.$project.kegg.funcover";
		if((-s $cogfuncover<1000) && (-s $keggfuncover<1000)) { die "Stopping in STEP11 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP12: Generation of the gene table
		
	if($rpoint<=12) {
		my $scriptname="12.mergeannot2.pl";
		print outfile3 "12\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP12 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP12 -> CREATING GENE TABLE: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $mergedfile<1000) { die "Stopping in STEP12 -> $scriptname\n"; }
	}
			
    #-------------------------------- STEP13: Running Maxbin (only for merged or coassembly modes)		
	
	if(($mode!~/sequential/i) && ($numsamples>1) && (!$nobins)) {	       
		if(($rpoint<=13) && (!$nomaxbin)) {
			my $scriptname="13.bin_maxbin.pl";
			print outfile3 "13\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP13 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP13 -> MAXBIN BINNING: $scriptname\n";
			system("perl $scriptdir/$scriptname $project >> $tempdir/$project.log");
			my $dirbin=$bindirs{maxbin};
			open(indir1,$dirbin);
			my @binfiles=grep(/maxbin.*fasta/,readdir indir1);
			closedir indir1;
			my $firstfile="$dirbin/$binfiles[0]";
			if(-s $firstfile<1000) { die "Stopping in STEP13 -> $scriptname\n"; }
		}
			
    #-------------------------------- STEP14: Running Metabat (only for merged or coassembly modes)		
	
		if(($rpoint<=14) && (!$nometabat)) {
			my $scriptname="14.bin_metabat2.pl";
			print outfile3 "14\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP14 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP14 -> METABAT BINNING: $scriptname\n";
			system("perl $scriptdir/$scriptname $project >> $tempdir/$project.log");
			my $dirbin=$bindirs{metabat2};
			open(indir2,$dirbin);
			my @binfiles=grep(/fasta/,readdir indir2);
			closedir indir2;
			my $firstfile="$dirbin/$binfiles[0]";
			if(-s $firstfile<1000) { die "Stopping in STEP14 -> $scriptname\n"; }
		}
			
    #-------------------------------- STEP15: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if($rpoint<=15) {
			my $scriptname="15.addtax2.pl";
			print outfile3 "15\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP15 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP15 -> BIN TAX ASSIGNMENT: $scriptname\n";
			system("perl $scriptdir/$scriptname $project >> $tempdir/$project.log");
			if(-s $bintax<1000) { die "Stopping in STEP15 -> $scriptname\n"; }
		}
			
    #-------------------------------- STEP16: Checking of bins for completeness and contamination (checkM)		
	
		if($rpoint<=16) {
			my $scriptname="16.checkM_batch.pl";
			print outfile3 "16\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP16 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP16 -> CHECKING BINS: $scriptname\n";
			system("perl $scriptdir/$scriptname $project >> $tempdir/$project.log");
			foreach my $binmethod(keys %bindirs) {
				$checkmfile="$resultpath/16.$project.$binmethod.checkM";
				if(-s $checkmfile<1000) { die "Cannot find $checkmfile\nStopping in STEP16 -> $scriptname\n"; }
				}
		}
			
    #-------------------------------- STEP17: Make bin table		
	
		if($rpoint<=17) {
			my $scriptname="17.getbins.pl";
			print outfile3 "17\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP17 -> $scriptname\n";
			print "[",$currtime->pretty,"]: STEP17 -> CREATING BIN TABLE: $scriptname\n";
			system("perl $scriptdir/$scriptname $project");
			if(-s $bintable<1000) { die "Stopping in STEP17 -> $scriptname\n"; }
		}
	}

    #-------------------------------- STEP18: Make contig table		

	if($rpoint<=18) {
		my $scriptname="18.getcontigs.pl";
		print outfile3 "18\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP18 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP18 -> CREATING CONTIG TABLE: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		if(-s $contigtable<1000) { die "Stopping in STEP18 -> $scriptname\n"; }
	}

    #-------------------------------- STEP19: Make stats		

	if($rpoint<=19) {
		my $scriptname="19.stats.pl";
		print outfile3 "19\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP19 -> $scriptname\n";
		print "[",$currtime->pretty,"]: STEP19 -> MAKING FINAL STATISTICS: $scriptname\n";
		system("perl $scriptdir/$scriptname $project");
		my $statfile="$resultpath/19.$project.stats";
		if(-s $statfile<1000) { die "Stopping in STEP19 -> $scriptname\n"; }
	}

    #-------------------------------- STEP20: Pathways in bins          

        if($rpoint<=20) {
                my $scriptname="20.minpath.pl";
                print outfile3 "20\t$scriptname\n";
                $currtime=timediff();
                print outfile4 "[",$currtime->pretty,"]: STEP20 -> $scriptname\n";
                print "[",$currtime->pretty,"]: STEP20 -> CREATING TABLE OF PATHWAYS IN BINS: $scriptname\n";
                system("perl $scriptdir/$scriptname $project");
                my $statfile="$resultpath/20.$project.kegg.pathways";
                if(-s $statfile<1000) { die "Stopping in STEP20 -> $scriptname\n"; }
        }

    #-------------------------------- END OF PIPELINE		

	print outfile3 "END\n";
	$currtime=timediff();
	print outfile4 "[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
	print "[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
}



