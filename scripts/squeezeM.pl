#!/usr/bin/perl

# v0.1.0 27/04/2018 Original version, (c) Javier Tamames, CNB-CSIC
# v0.1.1 07(05/2018 Added dynamic path detection. FPS.

$|=1;

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use strict;

my $version="0.1.2";
my $start_run = time();

###scriptdir patch, Fernando Puente-Sánchez, 07-V-2018
use File::Basename;
our $scriptdir = dirname(__FILE__);
our $installpath = "$scriptdir/..";
###

our $pwd=cwd();
our($nocog,$nokegg,$nopfam,$nobins)="0";
our($numsamples,$numthreads,$mode,$mincontiglen,$assembler,$project,$equivfile,$rawfastq,$evalue,$miniden,$spadesoptions,$megahitoptions,$assembler_options);
our($databasepath,$extdatapath,$softdir,$basedir,$datapath,$resultpath,$tempdir,$mappingfile,$contigsfna,$contigslen,$mcountfile,$rnafile,$gff_file,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$rpkmfile,$coveragefile,$contigcov,$contigtable,$mergedfile,$bintax,$checkmfile,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft);
our %bindirs;  

#-- Handle variables from command line

my $result = GetOptions ("t=i" => \$numthreads,
                     "m|mode=s" => \$mode,
                     "c|contiglen=i" => \$mincontiglen,
                     "a=s" => \$assembler,
                     "p=s" => \$project,
                     "s|samples=s" => \$equivfile,
                     "f|seq=s" => \$rawfastq, 
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nopfam" => \$nopfam,   
		     "nobins" => \$nobins,   
		     "e|evalue=f" => \$evalue,   
		     "minidentity=f" => \$miniden,   
		     "spades_options=s" => \$spadesoptions,
		     "megahit_options=s" => \$megahitoptions,
		    );

#-- Set some default values

if(!$numthreads) { $numthreads=12; }
if(!$mincontiglen) { $mincontiglen=1200; }
if(!$assembler) { $assembler="megahit"; }
if(!$evalue) { $evalue=1e-03; }
if(!$miniden) { $miniden=50; }
if(!$nocog) { $nocog=0; }
if(!$nokegg) { $nokegg=0; }
if(!$nopfam) { $nopfam=0; }
if(!$nobins) { $nobins=0; }

#-- Check if we have all the needed options

my $helptext="\nUsage: squeezeM.pl -m <mode> -p <projectname> -s <equivfile> -f <raw fastq dir> <options>\n\n Arguments:\n -m: Mode (sequential, coassembly, merged) (REQUIRED)\n -p: Project name (REQUIRED in coassembly and merged modes)\n -s|-samples: Samples file (REQUIRED)\n -f|-seq: Fastq read files' directory (REQUIRED)\n -t: Number of threads (Default:$numthreads)\n -a: assembler [megahit,spades] (Default:$assembler)\n -c|-contiglen: Minimum length of contigs (Default:$mincontiglen)\n --nocog: Skip COG assignment (Default: no)\n --nokegg: Skip KEGG assignment (Default: no)\n --nopfam: Skip Pfam assignment  (Default: no)\n --nobins: Skip binning  (Default: no)\n -e|-evalue: max evalue for diamond run  (Default: 1e-03)\n -miniden: minimum identity perc for diamond run  (Default: 50)\n --megahit_options: Options for megahit assembler\n --spades_options: Options for spades assembler\n"; 

print "\nSqueezeM v$version - (c) J. Tamames, CNB-CSIC, April 2018\n";

if((!$rawfastq) || (!$equivfile) || (!$mode)) { die "$helptext\n"; }
if(($mode!~/sequential/i) && (!$project)) { die "$helptext\n"; }
if(($mode=~/sequential/i) && ($project)) { die "$helptext\nPlease DO NOT specify project name in sequential mode. The name will be read from the samples in $equivfile\n"; }
if($mode!~/sequential|coassembly|merged/i) { die "$helptext\n"; }
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
		print outfile2 "Run started for $thissample, ",scalar localtime,"\n";
		print outfile4 "Project: $project\n";
		print outfile4 "Map file: $equivfile\n";
		print outfile4 "Fastq directory: $rawfastq\n";
		print outfile4 "[",$currtime->pretty,"]: STEP0 -> preparestuff_sequential.pl\n";
		print outfile2 "[",$currtime->pretty,"]: STEP0 -> preparestuff_sequential.pl\n";
	
		print "Now creating directories\n";
		open(infile2,"$scriptdir/squeezeM_conf.pl") || die;
	
		#-- Creation of the new configuration file for this sample
	
		open(outfile5,">$projectdir/squeezeM_conf.pl") || die;
		while(<infile2>) {
			if($_=~/^\$basedir/) { print outfile5 "\$basedir=\"$pwd\";\n"; }
			elsif($_=~/^\$installpath/) { print outfile5 "\$installpath=\"$installpath\";\n"; }
			elsif($_=~/^\$projectname/) { print outfile5 "\$projectname=\"$project\";\n"; }
			elsif($_=~/^\$evalue/) { print outfile5 "\$evalue=$evalue;\n"; }
			elsif($_=~/^\$miniden/) { print outfile5 "\$miniden=$miniden;\n"; }
			elsif($_=~/^\$nocog/) { print outfile5 "\$nocog=$nocog;\n"; }
			elsif($_=~/^\$nokegg/) { print outfile5 "\$nokegg=$nokegg;\n"; }
			elsif($_=~/^\$nopfam/) { print outfile5 "\$nopfam=$nopfam;\n"; }
			elsif($_=~/^\$nobins/) { print outfile5 "\$nobins=$nobins;\n"; }
			else { print outfile5 $_; }
        	}
	 	close infile2; 

 		if($assembler eq "megahit") { $assembler_options=$megahitoptions; } else { $assembler_options=$spadesoptions; }
 		print outfile5 "\n#-- Options\n\n\$numthreads=$numthreads;\n\$mincontiglen=$mincontiglen;\n\$assembler=$assembler;\n";
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
 		print "Now preparing files\n";
 		my($ca1,$ca2)="";
 		foreach my $afiles(sort keys %{ $ident{$thissample} }) {
 			next if($noassembly{$afiles});
  			my $gzfiles=$afiles;
  			#if($gzfiles!~/gz$/) { $gzfiles.=".gz"; }
 			if($ident{$thissample}{$afiles} eq "pair1") { $ca1.="$datapath/raw_fastq/$gzfiles "; $par1files++; } 
			else { $ca2.="$datapath/raw_fastq/$gzfiles "; $par2files++; } 
		}
		if($par1files>1) { 
			my $command="cat $ca1 > $datapath/raw_fastq/par1.fastq.gz"; 
			print "$command\n"; 
			system($command); 
		} 
		else { 
			my $command="cp $ca1 $datapath/raw_fastq/par1.fastq.gz"; 
			print "$command\n"; 
			system($command); 
		}
 		if($par2files>1) { system("cat $ca2 > $datapath/raw_fastq/par2.fastq.gz"); } 
		else { system("cp $ca2 $datapath/raw_fastq/par2.fastq.gz"); }
		
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
	print outfile4 "Run started ",scalar localtime," in $mode mode\n";
	print outfile4 "Command: $commandline\n"; 
	print outfile4 "Project: $project\n";
	print outfile4 "Map file: $equivfile\n";
	print outfile4 "Fastq directory: $rawfastq\n";
	print outfile4 "Options: threads=$numthreads; contiglen=$mincontiglen; assembler=$assembler;\n";
	print outfile4 "[",$currtime->pretty,"]: STEP0 -> preparestuff.pl\n";
     
	print "Now creating directories\n";
	
	#-- Creation of the new configuration file for this sample
	open(infile3,"$scriptdir/squeezeM_conf.pl") || die "Cannot open $scriptdir/squeezeM_conf.pl\n";
	open(outfile6,">$projectdir/squeezeM_conf.pl") || die;

	while(<infile3>) {
		if($_=~/^\$basedir/) { print outfile6 "\$basedir=\"$pwd\";\n"; }
		elsif($_=~/^\$installpath/) { print outfile6 "\$installpath=\"$installpath\";\n"; }
		elsif($_=~/^\$projectname/) { print outfile6 "\$projectname=\"$project\";\n"; }
		elsif($_=~/^\$evalue/) { print outfile6 "\$evalue=$evalue;\n"; }
		elsif($_=~/^\$miniden/) { print outfile6 "\$miniden=$miniden;\n"; }
		elsif($_=~/^\$nocog/) { print outfile6 "\$nocog=$nocog;\n"; }
		elsif($_=~/^\$nokegg/) { print outfile6 "\$nokegg=$nokegg;\n"; }
		elsif($_=~/^\$nopfam/) { print outfile6 "\$nopfam=$nopfam;\n"; }
		elsif($_=~/^\$nobins/) { print outfile6 "\$nobins=$nobins;\n"; }
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
	print outfile3 "Samples:$numsamples\nMode:coassembly\n0\n";
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
		if($par2files>1) { system("cat $ca2 > $par2name"); } else { system("ln -s $ca2 $par2name"); }
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
		print outfile3 "1\trun_assembly.pl ($assembler)\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1 -> run_assembly.pl ($assembler)\n";
		print "[",$currtime->pretty,"]: STEP1 -> run_assembly.pl ($assembler)\n";
		system("perl $scriptdir/run_assembly.pl $project ");
		if((-e $contigsfna) && (!(-z $contigsfna))) {} else { die "Assembly didn't provide results\n"; }
	}

		#-- In merged mode. Includes merging assemblies

	elsif($mode=~/merged/) {
		print outfile3 "1\trun_assembly_merged.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1 -> run_assembly_merged.pl\n";
		print "[",$currtime->pretty,"]: STEP1 -> run_assembly_merged.pl\n";
		system("perl $scriptdir/run_assembly_merged.pl $project");
	
			#-- Merging individual assemblies
 
		print outfile3 "1.5\tmerge_assemblies.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP1.5 -> merge_assemblies.pl\n";
		print "[",$currtime->pretty,"]: STEP1.5 -> merge_assemblies.pl\n";
		system("perl $scriptdir/merge_assemblies.pl $project");
	}
	
		#-- In sequential mode. 
     
	elsif($mode=~/sequential/) {
 		print outfile3 "1\trun_assembly.pl\n";
 		$currtime=timediff();
 		print outfile4 "[",$currtime->pretty,"]: STEP1 -> run_assembly.pl\n";
 		print "[",$currtime->pretty,"]: STEP1 -> run_assembly.pl\n";
 		system("perl $scriptdir/run_assembly.pl $project");
	}		
			
    #-------------------------------- STEP2: Run RNA prediction

	if($rpoint<=2) {
		print outfile3 "2\trun_barrnap.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP2 -> run_barrnap.pl\n";
		print "[",$currtime->pretty,"]: STEP2 -> run_barrnap.pl\n";
		system("perl $scriptdir/run_barrnap.pl $project");
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if($rpoint<=3) {
 		print outfile3 "3\trun_prodigal.pl\n";
 		$currtime=timediff();
 		print outfile4 "[",$currtime->pretty,"]: STEP3 -> run_prodigal.pl\n";
 		print "[",$currtime->pretty,"]: STEP3 -> run_prodigal.pl\n";
 		system("perl $scriptdir/run_prodigal.pl $project");
		 if((-e $aafile) && (!(-z $aafile))) {} else { die "Prodigal didn't provide results\n"; }
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if($rpoint<=4) {
		print outfile3 "4\trundiamond.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP4 -> rundiamond.pl\n";
		print "[",$currtime->pretty,"]: STEP4 -> rundiamond.pl\n";
		system("perl $scriptdir/rundiamond.pl $project");
		if((-e $taxdiamond)  && (!(-z $taxdiamond))) {} else { die "Diamond didn't provide results\n"; }
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if($rpoint<=5) {
		if(!$nopfam) {
			print outfile3 "5\trunhmmer.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP5 -> runhmmer.pl\n";
			print "[",$currtime->pretty,"]: STEP5 -> runhmmer.pl\n";
			system("perl $scriptdir/run_hmmer.pl $project");
			if((-e $pfamhmmer)  && (!(-z $pfamhmmer))) {} else { die "hmmsearch didn't provide results\n"; }
		}
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if($rpoint<=6) {
		print outfile3 "6\tlca.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP6 -> lca.pl\n";
		print "[",$currtime->pretty,"]: STEP6 -> lca.pl\n";
		system("perl $scriptdir/lca.pl $project");
		if((-e $fun3tax)  && (!(-z $fun3tax))) {} else { die "lca didn't provide results\n"; }
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if($rpoint<=7) {
		if((!$nocog) || (!$nokegg) || (!$nopfam)) {
		print outfile3 "7\tfun3assign.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP7 -> fun3assign.pl\n";
		print "[",$currtime->pretty,"]: STEP7 -> fun3assign.pl\n";
		system("perl $scriptdir/fun3assign.pl $project");
		if((-e $fun3cog)  && (!(-z $fun3cog))) {} else { die "fun3assign didn't provide results\n"; }
		}
	}
			
    #-------------------------------- STEP8: Taxonomic annotation for the contigs (consensus of gene annotations)

	if($rpoint<=8) {
		print outfile3 "8\tsummarycontigs3.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP8 -> summarycontigs3.pl\n";
		print "[",$currtime->pretty,"]: STEP8 -> summarycontigs3.pl\n";
		system("perl $scriptdir/summarycontigs3.pl $project");
		if((-e $alllog)  && (!(-z $alllog))) {} else { die "summarycontigs didn't provide results\n"; }
	}
			
    #-------------------------------- STEP9: Mapping of reads onto contigs for abundance calculations
	
	if($rpoint<=9) {
		print outfile3 "9\tmapbamsamples.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP9 -> mapbamsamples.pl\n";
		print "[",$currtime->pretty,"]: STEP9 -> mapbamsamples.pl\n";
		system("perl $scriptdir/mapbamsamples.pl $project");
		if((-e $rpkmfile)  && (!(-z $rpkmfile))) {} else { die "mapbamsamples didn't provide results\n"; }
	}
			
    #-------------------------------- STEP10: Count of taxa abundances
	
	if($rpoint<=10) {
		print outfile3 "10\tmcount.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP9 -> mcount.pl\n";
		print "[",$currtime->pretty,"]: STEP9 -> mcount.pl\n";
		system("perl $scriptdir/mcount.pl $project");
		if((-e $mcountfile)  && (!(-z $mcountfile))) {} else { die "mcount didn't provide results\n"; }
	}
			
    #-------------------------------- STEP11: Count of function abundances
	
	if($rpoint<=11) {
		print outfile3 "11\tfuncover.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP11 -> funcover.pl\n";
		print "[",$currtime->pretty,"]: STEP11 -> funcover.pl\n";
		system("perl $scriptdir/funcover.pl $project");
		# if((-e $funcoverfile)  && (!(-z $funcoverfile))) {} else { die "funcover didn't provide results\n"; }
	}
			
    #-------------------------------- STEP12: Generation of the gene table
		
	if($rpoint<=12) {
		print outfile3 "12\tmergeannot2.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP12 -> mergeannot2.pl\n";
		print "[",$currtime->pretty,"]: STEP12 -> mergeannot2.pl\n";
		system("perl $scriptdir/mergeannot2.pl $project");
		if((-e $mergedfile) && (!(-z $mergedfile))) {} else { die "mergeannot didn't provide results\n"; }
	}
			
    #-------------------------------- STEP13: Running Maxbin (only for merged or coassembly modes)		
	
	if(($mode!~/sequential/i) && ($numsamples>1) && (!$nobins)) {	       
		if($rpoint<=13) {
			print outfile3 "13\tbin_maxbin.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP13 -> bin_maxbin.pl\n";
			print "[",$currtime->pretty,"]: STEP13 -> bin_maxbin.pl\n";
			system("perl $scriptdir/bin_maxbin.pl $project >> $tempdir/$project.log");
			my $dirbin=$bindirs{maxbin};
			my $firstfile="$dirbin/maxbin.001.fasta";
			# if((-e $firstfile) && (!(-z $firstfile))) {} else { die "bin_maxbin didn't provide results\n"; }
		}
			
    #-------------------------------- STEP14: Running Metabat (only for merged or coassembly modes)		
	
		if($rpoint<=14) {
			print outfile3 "14\tbin_metabat2.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP14 -> bin_metabat2.pl\n";
			print "[",$currtime->pretty,"]: STEP14 -> bin_metabat2.pl\n";
			system("perl $scriptdir/bin_metabat2.pl $project >> $tempdir/$project.log");
			my $dirbin=$bindirs{metabat2};
			my $firstfile="$dirbin/metabat2.1.fa";
			# if((-e $firstfile) && (!(-z $firstfile))) {} else { die "bin_metabat2 didn't provide results\n"; }
		}
			
    #-------------------------------- STEP15: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if($rpoint<=15) {
			print outfile3 "15\taddtax2.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP15 -> addtax2.pl\n";
			print "[",$currtime->pretty,"]: STEP15 -> addtax2.pl\n";
			system("perl $scriptdir/addtax2.pl $project >> $tempdir/$project.log");
			if((-e $bintax)  && (!(-z $bintax))) {} else { die "addtax2 didn't provide results\n"; }
		}
			
    #-------------------------------- STEP16: Checking of bins for completeness and contamination (checkM)		
	
		if($rpoint<=16) {
			print outfile3 "16\tcheckM_batch.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP16 -> checkM_batch.pl\n";
			print "[",$currtime->pretty,"]: STEP16 -> checkM_batch.pl\n";
			system("perl $scriptdir/checkM_batch.pl $project >> $tempdir/$project.log");
			# if((-e $checkmfile)  && (!(-z $checkmfile))) {} else { die "checkM_batch didn't provide results\n"; }
		}
			
    #-------------------------------- STEP17: Make bin table		
	
		if($rpoint<=17) {
			print outfile3 "17\tgetbins.pl\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP17 -> getbins.pl\n";
			print "[",$currtime->pretty,"]: STEP17 -> getbins.pl\n";
			system("perl $scriptdir/getbins.pl $project");
			if((-e $bincov)  && (!(-z $bincov))) {} else { die "getbins didn't provide results\n"; }
		}
	}

    #-------------------------------- STEP18: Make contig table		

	if($rpoint<=18) {
		print outfile3 "18\tgetcontigs.pl\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP18 -> getcontigs.pl\n";
		print "[",$currtime->pretty,"]: STEP18 -> getcontigs.pl\n";
		system("perl $scriptdir/getcontigs.pl $project");
		# if((-e $contigsinbins)  && (!(-z $contigsinbins))) {} else { die "getcontigs didn't provide results\n"; }
	}

    #-------------------------------- END OF PIPELINE		

	print outfile3 "END\n";
	$currtime=timediff();
	print outfile4 "[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
	print "[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
}



