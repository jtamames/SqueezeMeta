#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use Term::ANSIColor qw(:constants);
use lib ".";
use strict;

my $start_run = time();

my $verbose=0;    #-- Reports an explanation msg for each of the steps

###scriptdir patch v2, Fernando Puente-SÃ¡nchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

our $scriptdir;
if(-l __FILE__)
	{
	my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $scriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
	{
	$scriptdir = abs_path(dirname(__FILE__));
	}
our $installpath = abs_path("$scriptdir/..");
###
$scriptdir="$installpath/scripts";

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

our $pwd=cwd();

our($nodiamond,$fastnr,$binners,$nocog,$nokegg,$nopfam,$singletons,$euknofilter,$opt_db,$nobins,$onlybins,$nomaxbin,$nometabat,$empty,$verbose,$lowmem,$minion,$consensus,$doublepass,$force_overwrite)="0";
our($numsamples,$numthreads,$canumem,$mode,$mincontiglen,$contigid,$assembler,$extassembly,$extbins,$mapper,$projectdir,$userdir,$mapping_options,$projectname,$project,$equivfile,$rawfastq,$blocksize,$evalue,$miniden,$assembler_options,$cleaning,$cleaningoptions,$ver,$hel,$methodsfile,$test,$norename,$restart,$rpoint);
our($binresultsdir,$databasepath,$extdatapath,$newtaxdb,$softdir,$datapath,$resultpath,$extpath,$tempdir,$interdir,$mappingfile,$protclust,$extdatapath,$contigsfna,$gff_file_blastx,$contigslen,$mcountfile,$checkmfile,$rnafile,$gff_file,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$mapcountfile,$mappingstat,$contigcov,$contigtable,$mergedfile,$bintax,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bwa_soft,$minimap2_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft,$canu_soft,$trimmomatic_soft,$dastool_soft,$taxbinmode,$gtdbtk,$gtdbtk_data_path,$gtdbtkfile);
our(%bindirs,%dasdir,%binscripts,%assemblers);

#-- Load the database path from SqueezeMeta_conf.pl so we can use it to check that the databases are in place
eval(`grep "^\\\$databasepath" $scriptdir/SqueezeMeta_conf.pl`);

#-- Define help text

my $helpshort = <<END_MESSAGE;
Usage: SqueezeMeta.pl -m <mode> -p <project name> -s <samples file> -f <sequence dir> [options]

END_MESSAGE

my $helptext = <<END_MESSAGE;
Usage: SqueezeMeta.pl -m <mode> -p <project name> -s <samples file> -f <sequence dir> [options]

Arguments:

 Mandatory parameters:
   -m <mode>: Mode (sequential, coassembly, merged, seqmerge) (REQUIRED)
   -s|-samples <samples file>: Samples file (REQUIRED)
   -f|-seq <sequence dir>: fastq/fasta read files' directory (REQUIRED)
   -p <project name>: Project name (REQUIRED in coassembly and merged modes)
  
 Restarting
   --restart: Restarts the given project where it stopped (project must be speciefied with -p option) (will NOT overwite previous results, unless --force-overwrite is also provided)
   -step <step number>: In combination with --restart, restarts the project starting in the given step number (combine with --force_overwrite to regenerate results)
   --force_overwrite: Do not check for previous results, and overwrite existing ones
   
 Filtering: 
   --cleaning: Filters with Trimmomatic (Default: No)
   -cleaning_options [options]: Options for Trimmomatic (Default:LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30)
   
 Assembly: 
   -a: assembler <megahit, spades, rnaspades, spades-base, canu, flye> (Default: megahit)
   -assembly_options [options]: Extra options to be passed when calling the mapper
   -c|-contiglen <size>: Minimum length of contigs (Default: 200)
   -extassembly <file>: External assembly, path to a fasta file with contigs (overrides all assembly steps).
   --sg|--singletons: Add unassembled reads to the contig file, as if they were contigs  
   -contigid <string>: Nomenclature for contigs (Default: assembler´s name)
   --norename: Don't rename contigs (Use at your own risk, characters like '_' in contig names will make it crash)
   
 Mapping: 
   -map: mapping software <bowtie, bwa, minimap2-ont, minimap2-pb, minimap2-sr> (Default: bowtie) 
   -mapping_options [options]: Extra options to be passed when calling the mapper

 ONT support: 
   --minion: Run on MinION reads (assembler: canu; mapper: minimap2-ont; consensus: 20) (Default: no)

 Annotation:  
   -db <file>: Specify a new taxonomic database
   --nodiamond: Check if Diamond results are already in place, and just in that case skips the Diamond run (Default: no)
   --nocog: Skip COG assignment (Default: no)
   --nokegg: Skip KEGG assignment (Default: no)
   --nopfam: Skip Pfam assignment  (Default: no)
   --fastnr: Run DIAMOND in --fast mode for taxonomic assignment (Default: no)
   --euk: Drop identity filters for eukaryotic annotation (Default: no)
   -consensus <value>: Minimum percentage of genes for a taxon needed for contig consensus (Default: 50)
   --D|--doublepass: First pass looking for genes using gene prediction, second pass using Diamond BlastX  (Default: no)
   -extdb <database file>: List of user-provided databases
   -b|-block-size <block size>: block size for diamond against the nr database (Default: calculate automatically)
   
 Binning:
   --nobins: Skip all binning  (Default: no). Overrides -binners
   --onlybins: Run only assembly, binning and bin statistics (including GTDB-Tk if requested) (Default: no)
   -binners: Comma-separated list with the binning programs to be used (available: maxbin, metabat2, concoct)  (Default: concoct,metabat2)
   -taxbinmode <s,c,s+c,c+s>: Source of taxonomy annotation of bins (s: SqueezeMeta; c: CheckM; s+c: SqueezeMeta+CheckM;  c+s: CheckM+SqueezeMeta; (Default: s)
   --gtdbtk: Run GTDB-Tk to classify the bins. Requires a working GTDB-Tk installation in available in your environment
   -gtdbtk_data_path: Path to the GTDB database, by default it is assumed to be present in /path/to/SqueezeMeta/db/gtdb
   -extbins: Path to a directory containing external genomes/bins provided by the user. There must be one file per genome/bin, containing each contigs in the fasta format. This overrides the assembly and binning steps

 
 Performance:
   -t <threads>: Number of threads (Default: 12)
   -canumem <mem>: memory for canu in Gb (Default: 32)
   --lowmem: run on less than 16Gb of memory (Default:no)

 Other:
   -test <step>: Running in test mode, stops AFTER the given step number
   --empty: Creates a empty directory structure and conf files, does not run the pipeline
   
 Information:
   -v: Version number  
   -h: This help 
     
END_MESSAGE

#-- Handle variables from command line

my $result = GetOptions ("t=i" => \$numthreads,
                     "lowmem" => \$lowmem,
		     "canumem=i" => \$canumem,
                     "m|mode=s" => \$mode,
                     "c|contiglen=i" => \$mincontiglen,
                     "contigid=s" => \$contigid,  
	             "a=s" => \$assembler,
                     "map=s" => \$mapper,
                     "p=s" => \$projectdir,
                     "s|samples=s" => \$equivfile,
                     "extassembly=s" => \$extassembly,
                     "f|seq=s" => \$rawfastq,
		     "nodiamond" => \$nodiamond,
                     "sg|singletons" => \$singletons,
		     "db=s" => \$newtaxdb,
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nopfam" => \$nopfam,
		     "fastnr" => \$fastnr,  
		     "euk" => \$euknofilter,
		     "protclust" => \$protclust,
		     "extdb=s" => \$opt_db, 
		     "nobins" => \$nobins,
		     "onlybins" => \$onlybins,
		     "binners=s" => \$binners, 
		     "taxbinmode=s" => \$taxbinmode,
		     "gtdbtk" => \$gtdbtk,
		     "gtdbtk_data_path=s" => \$gtdbtk_data_path,
		     "extbins=s" => \$extbins,
		     "D|doublepass" => \$doublepass, 
		     "b|block_size=i" => \$blocksize,
		     "e|evalue=f" => \$evalue,   
		     "minidentity=f" => \$miniden,   
		     "assembly_options=s" => \$assembler_options,
		     "norename" => \$norename,
		     "restart" => \$restart,
		     "step=i" => \$rpoint,
		     "force_overwrite" => \$force_overwrite,
		     "cleaning" => \$cleaning,
		     "cleaning_options=s" => \$cleaningoptions,
		     "mapping_options=s" => \$mapping_options,
                     "minion" => \$minion,
		     "test=i" => \$test,
		     "empty" => \$empty,
		     "verbose" => \$verbose,
		     "v" => \$ver,
		     "h" => \$hel
		    );

#-- Set some default values

if(!$numthreads) { $numthreads=12; }
if(!$canumem) { $canumem="NF"; }
if(!$mincontiglen) { $mincontiglen=200; }
if(!$assembler) { $assembler="megahit"; }
if(!$mapper) { $mapper="bowtie"; }
if(!$blocksize) { $blocksize="NF"; }
if(!$nodiamond) { $nodiamond=0; }
if(!$singletons) { $singletons=0; }
if(!$nocog) { $nocog=0; }
if(!$nokegg) { $nokegg=0; }
if(!$nopfam) { $nopfam=0; }
if(!$fastnr) { $fastnr=0; }
if(!$euknofilter) { $euknofilter=0; }
if(!$doublepass) { $doublepass=0; }
if(!$nobins) { $nobins=0; }
if(!$onlybins) { $onlybins = 0; }
if(!$gtdbtk) { $gtdbtk=0; }
if(!$gtdbtk_data_path) { eval(`grep "^\\\$gtdbtk_data_path" $scriptdir/SqueezeMeta_conf.pl`); }
if(!$binners) { $binners="concoct,metabat2"; }
if(!$taxbinmode) { $taxbinmode="s"; }
if(!$nomaxbin) { $nomaxbin=0; }
if(!$nometabat) { $nometabat=0; }
if(!$norename) { $norename=0; }
if(!$force_overwrite) { $force_overwrite=0; }
if(!$restart) { $restart=0; }
if(!$cleaningoptions) { $cleaningoptions="LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30"; }
if(!$cleaning) { $cleaning=0; $cleaningoptions=""; } 
if($consensus) { $consensus/=100; }

$mode=~tr/A-Z/a-z/;
if($opt_db) { $opt_db = abs_path($opt_db); }
if($newtaxdb) { 
	$newtaxdb = abs_path($newtaxdb);
	if($newtaxdb!~/\.dmnd$/) { $newtaxdb.=".dmnd"; }
	}

#-- Override settings if running on lowmem or MinION mode.
if($lowmem) { $blocksize=3; $canumem=15; }

if($minion) { $assembler="canu"; $mapper="minimap2-ont"; }


#-- Check if we have all the needed options
my($dietext,$finaltrace);
if($ver) { print "$version\n"; exit; }
if($hel) { die "$helptext\n"; } 

#-- If we are restarting, just load the SqueezeMeta_conf.pl (all needed information is there)

if($restart) {
	if(!$projectdir) { $dietext.="MISSING ARGUMENT: -p: Project name to restart\n"; }
	do "$projectdir/SqueezeMeta_conf.pl";
	$equivfile="$projectdir/data/00.$projectname.samples";
	$rawfastq=$userdir;
	}
else {
	if(!$rawfastq)  { $dietext.="MISSING ARGUMENT: -f|-seq: Fastq read files' directory\n"; }
	if(!$equivfile) { $dietext.="MISSING ARGUMENT: -s|-samples: Samples file\n"; }
	if(!$mode)      { $dietext.="MISSING ARGUMENT: -m: Run mode (sequential, coassembly, merged)\n"; }
	if(($mode!~/sequential$/i) && (!$projectdir)) { $dietext.="MISSING ARGUMENT: -p: Project name\n"; }
	if(($mode=~/sequential$/i) && ($projectdir))  { $dietext.="Please DO NOT specify project name in sequential mode. The name will be read from the samples in the samples file $equivfile\n"; }
	if($mode!~/sequential|coassembly|merged|seqmerge/i) { $dietext.="UNRECOGNIZED mode $mode (valid ones are sequential, coassembly, merged or seqmerge\n"; }
	if($mapper!~/bowtie|bwa|minimap2-ont|minimap2-pb|minimap2-sr/i) { $dietext.="UNRECOGNIZED mapper $mapper (valid ones are bowtie, bwa, minimap2-ont, minimap2-pb or minimap2-sr\n"; }
	# if($assembler!~/megahit|spades|rnaspades|canu|flye/i) { $dietext.="UNRECOGNIZED assembler $assembler (valid ones are megahit, spades, canu or flye)\n"; }
	if($newtaxdb) { if(-e "$newtaxdb.dmnd") {}  else { $dietext.="New taxonomy database specified in $newtaxdb not found\n"; } }
	if($extassembly && $extbins) { $dietext.="-extassembly and -extbins can not be provided at the same time\n"; }
        if($nobins and $onlybins)    { $dietext.="--nobins --onlybins can not be provided at the same time\n"; }
	if($rawfastq=~/^\//) {} else { $rawfastq=abs_path($rawfastq); }
	if($gtdbtk) {
		if(! -e "$gtdbtk_data_path/fastani/genome_paths.tsv") {
			$dietext.="--gtdbtk was provided but we can't find the GTDB-Tk database at $gtdbtk_data_path. Please provide the right path to the database through the -gtdbtk_data_path argument\n"; 
			$dietext.="Note that the GTDB-Tk database is not provided as part of SqueezeMeta, and needs to be downloaded separetely\n";
			$dietext.="More info can be found at https://ecogenomics.github.io/GTDBTk/installing/index.html\n";
			}
		}
	 if($dietext) { print BOLD "$helpshort"; print RESET; print RED; print "$dietext"; print RESET;  exit; } 
	}

#------------------------------------- CHECKING FILES AND START RUN -----------------------------------------------

$projectdir = abs_path($projectdir);
$projectname = (split '/', $projectdir)[-1];
my $syslogfile="$projectdir/syslog";
if (($mode!~/sequential$/i) && (-d $projectdir) && (!$restart)) { print RED; print "Project name $projectdir already exists. Please remove it or change the project name\n"; print RESET; die; } 
elsif(!$restart && $mode ne "sequential") { system("mkdir $projectdir"); }

my(%allsamples,%ident,%noassembly,%pairsample);
my($sample,$file,$iden,$mapreq);
tie %allsamples,"Tie::IxHash";


	#-- Check that everything is correct in the samples file

open(infile1,$equivfile) || die "Can't open samples file (-s) in $equivfile. Please check that it is the correct file\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	$_=~s/\r//g;
	my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
	$pairsample{$sample}.="$iden;";
	if($_=~/ /) { print RED; print "Please do not use blank spaces in the samples file\n"; print RESET;  exit; }
	if(($iden ne "pair1") && ($iden ne "pair2")) { print RED; print "Samples file, line $_: file label must be \"pair1\" or \"pair2\". For single reads, use \"pair1\"\n"; print RESET;  exit; }
	if((!$sample) || (!$file) || (!$iden)) { print RED; print "Bad format in samples file $equivfile. Missing fields\n"; print RESET;  exit; }
	if(-e "$rawfastq/$file") {} elsif(!$empty) { print RED; print "Can't find sample file $rawfastq/$file for sample $sample in the samples file. Please check\n"; print RESET;  exit; }
	$allsamples{$sample}{$file}=1;
	$ident{$sample}{$file}=$iden;
}
close infile1;
foreach my $chsam(keys %pairsample) { 
	if($pairsample{$chsam}!~/pair1/) { print RED; print "Sample $chsam has not pair1 in the samples file. Please check\n"; print RESET;  exit; }
	}

my $currtime=timediff();
print BOLD "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sánchez, Frontiers in Microbiology 9, 3349 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;
if($test) { print GREEN "Running in test mode! I will stop after step $test\n\n"; print RESET; }
print "Run started ",scalar localtime," in $mode mode\n";

my @nmg=keys %allsamples;
$numsamples=$#nmg+1;
print "$numsamples metagenomes found: @nmg";
print "\n\n";


#------------------------------------- WRITE CONF AND START PIPELINE (v1.6) -----------------------------------------------



my $pdir;

my %conf=('version',$version,'mode',$mode,'installpath',$installpath,'projectname',$projectname,'userdir',$rawfastq,
  'blocksize',$blocksize,'nodiamond',$nodiamond,'singletons',$singletons,'nocog',$nocog,'nokegg',$nokegg,
  'nopfam',$nopfam,'fastnr',$fastnr,'euknofilter',$euknofilter,'doublepass',$doublepass,'nobins',$nobins,'onlybins',$onlybins,'binners',$binners,'gtdbtk',$gtdbtk,'gtdbtk_data_path', $gtdbtk_data_path,
  'norename',$norename,'mapper',$mapper,'mapping_options',$mapping_options,'cleaning',$cleaning,
  'cleaningoptions',$cleaningoptions,'consensus',$consensus,'numthreads',$numthreads,'mincontiglen',$mincontiglen,
  'assembler',$assembler,'canumem',$canumem,'contigid',$contigid,'assembler_options',$assembler_options,
  'extassembly',$extassembly,'extbins',$extbins, 'opt_db',$opt_db,'samples',$equivfile,'commandline',$commandline,
  'miniden',$miniden, 'evalue',$evalue,'taxbinmode',$taxbinmode,'overwrite',$force_overwrite,'newtaxdb',$newtaxdb);


if($mode!~/sequential/) {   #-- FOR ALL COASSEMBLY AND MERGED MODES
		
	#-- Creation of the new configuration file, syslog, directories
		  
	if(!$restart) { 
		writeconf($projectdir,$scriptdir,%conf); 
		if($cleaning) { cleaning($projectdir,$scriptdir,"",%conf); }
		}
	  		

	if(!$empty) {
 	
		#-- CALL TO THE STANDARD PIPELINE
		if($cleaning) { cleaning($projectdir,$scriptdir,"",%conf); }
		pipeline();
		}
		
		
	else { die "  Directory structure and conf files created. Exiting\n"; }  #-- If --empty invoked
	close outfile4;  #-- Closing log file for the sample
	close outfile3;	 #-- Closing progress file for the sample
	}

	#-- Sequential mode
		
else {      #-- FOR SEQUENTIAL MODE
	
	my $rootdir=$projectdir;
		
	foreach my $thissample(keys %allsamples) { 
	
		print "--- SAMPLE $thissample ---\n";
	
	
		#-- Creation of the new configuration file, syslog, and directories
	
	  
		if(!$restart) {
	                $projectdir="$rootdir/$thissample";
	                $conf{'projectname'}=$thissample;
			writeconf($projectdir,$scriptdir,%conf); 
			if($cleaning) { cleaning($projectdir,$scriptdir,$thissample,%conf); }
			} 

		
		if(!$empty) {
 	
			#-- CALL TO THE STANDARD PIPELINE
			
			pipeline();
			
			}
		else { die "  Directory structure and conf files created. Exiting\n"; }  #-- If --empty invoked
		close outfile4;  #-- Closing log file for the sample
		close outfile3;	  #-- Closing progress file for the sample
	 		
				
		}	


}                        #------ END



#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}

#---------------------------------------- PIPELINE --------------------------------

sub pipeline {

	if(-e "$tempdir/$projectname.log") { system("rm $tempdir/$projectname.log"); }
	my $DAS_Tool_empty=0;


    #-------------------------------- STEP1: Run assembly

	if(($rpoint<=1) && ((!$test) || ($test>=1))) {
		if(!$assemblers{$assembler}) { 
			my $assemblerlist=join(", ",keys %assemblers);
			$dietext.="UNRECOGNIZED assembler $assembler (valid ones are $assemblerlist)\n"; 
		}
		if(($assembler=~/flye/i) && ($mode=~/merge/i)) { $dietext.="Invalid combination of mode and assembler\n (We are sorry for this, the low number of contigs provided by Flye prevents minimus2 needed in $mode mode to work correctly\n Please use coassembly, or a different assembler)\n"; }
		if($dietext) { print BOLD "$helpshort"; print RESET; print RED; print "$dietext"; print RESET;  exit; }


		my $scriptname="01.run_all_assemblies.pl";
                my $wsize=checksize($contigsfna);
		my $wsize2=checksize($contigslen);
                if(($wsize>=2) && ($wsize2>=1) && (!$force_overwrite)) { print "Contig file $contigsfna already found, skipping step 1\n"; }
		else {		
			print outfile3 "1\t$scriptname ($assembler)\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n"; 
			print BLUE "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)	{ error_out(1,$scriptname); }
			my $wc=qx(wc -l $contigsfna);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($singletons) {
				my $scriptname="01.remap.pl";
               			print outfile3 "1\t$scriptname\n";
                		$currtime=timediff();
                		print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname\n";
                		print BLUE "[",$currtime->pretty,"]: STEP1 ->  ADDING SINGLETONS: $scriptname ($assembler)\n"; print RESET;
                		if($verbose) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
                		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                		if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n"; print RESET; die; }
                		my $wsize=checksize($contigsfna);
                		if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname ($assembler). File $contigsfna is empty!\n"; print RESET; die; }
			}
   	        my $wsize=checksize($contigsfna);
		if($wsize<1)	{ error_out(1,$scriptname,$contigsfna); }
		}

	close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP2: Run RNA prediction

	if(($rpoint<=2) && (!$onlybins) && ((!$test) || ($test>=2))) {
		my $masked="$interdir/02.$projectname.maskedrna.fasta";
                my $wsize=checksize($masked);
                if(($wsize>=2) && (!$force_overwrite)) { print "RNA gff file $masked already found, skipping step 2\n"; }
		else {		
			if($verbose) { print " At this point, we already have contigs\n"; }
			my $scriptname="02.rnas.pl";
			print outfile3 "2\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP2 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP2 -> RNA PREDICTION: $scriptname\n"; print RESET;
			if($verbose) { print " (This will run barrnap and Aragorn for predicting putative RNAs in the contigs. This is done before predicting protein-coding genes for avoiding predicting these where there is a RNA)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(2,$scriptname); }
               		my $wsize=checksize($masked);
			if($wsize<2)    { error_out(2,$scriptname,$masked); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}	
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if(($rpoint<=3) && (!$onlybins) && ((!$test) || ($test>=3))) { 
                my $wsize=checksize($aafile);
                if(($wsize>=2) && (!$force_overwrite))  { print "Aminoacid file $aafile already found, skipping step 3\n"; }
		else {		
			my $scriptname="03.run_prodigal.pl";
 			print outfile3 "3\t$scriptname\n";
 			$currtime=timediff();
 			print outfile4 "\n[",$currtime->pretty,"]: STEP3 -> $scriptname\n";
 			print BLUE "[",$currtime->pretty,"]: STEP3 -> ORF PREDICTION: $scriptname\n"; print RESET;
			if($verbose) { print " (This will predict putative protein-coding genes in the contigs, using Prodigal)\n"; }
 			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(3,$scriptname); }
			my $wsize=checksize($aafile);
			if($wsize<2)    { error_out(3,$scriptname,$aafile); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if(($rpoint<=4) && (!$onlybins) && ((!$test) || ($test>=4))) {
                my $wsize=checksize($taxdiamond);
                if(($wsize>=1) && (!$force_overwrite)) { print "Diamond file $taxdiamond already found, skipping step 4\n"; }
		else {		
			if($verbose) { print " At this point, we already have ORFs\n"; }
			my $scriptname="04.rundiamond.pl";
			print outfile3 "4\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP4 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP4 -> HOMOLOGY SEARCHES: $scriptname\n"; print RESET;
			if($verbose) { print " (This will take all ORFs and run homology searches against functional and taxonomic databases)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(4,$scriptname); }
			my $wsize=checksize($taxdiamond);
			if($wsize<1)    { error_out(4,$scriptname,$taxdiamond); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if(($rpoint<=5) && (!$onlybins) && ((!$test) || ($test>=5))) {
		if(!$nopfam) {
            		my $wsize=checksize($pfamhmmer);
             		if(($wsize>=1) && (!$force_overwrite)) { print "Pfam file $pfamhmmer already found, skipping step 5\n"; }	
			else {	
				my $scriptname="05.run_hmmer.pl";
				print outfile3 "5\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "\n[",$currtime->pretty,"]: STEP5 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP5 -> HMMER/PFAM: $scriptname\n"; print RESET;
				if($verbose) { print " (This will take all ORFs and run HMMER searches against the PFAM database, for increased sensitivity in predicting protein families. This step is likely to be slow)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0)        { error_out(5,$scriptname); }			
				my $wsize=checksize($pfamhmmer);
				if($wsize<1)    { error_out(5,$scriptname,$pfamhmmer); }
				}
		}
	outfile4->autoflush;		
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if(($rpoint<=6) && (!$onlybins) && ((!$test) || ($test>=6))) {
		my $lcaresult="$fun3tax.wranks";
            	my $wsize=checksize($lcaresult);
             	if(($wsize>=1) && (!$force_overwrite)) { print "LCA file $lcaresult already found, skipping step 6\n"; }
		else {		
			my $scriptname="06.lca.pl";
			print outfile3 "6\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP6 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP6 -> TAXONOMIC ASSIGNMENT: $scriptname\n"; print RESET;
			if($verbose) { print " (This will use our last common ancestor (LCA) algorithm to try to annotate the taxonomic origin of each ORF, from the homologues found in the previous step)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(6,$scriptname); }
			my $wsize=checksize($lcaresult);
			if($wsize<1)    { error_out(6,$scriptname,$lcaresult); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}		
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if(($rpoint<=7) && (!$onlybins) && ((!$test) || ($test>=7))) {
		# Set wsize to 1 for a method if we were instructed not to use it,
		#  so that the step counts as completed even if those results are not present.
		my($wsizeCOG,$wsizeKEGG,$wsizePFAM,$wsizeOPTDB);
		if($nocog)  { $wsizeCOG  = 1; } else { $wsizeCOG  = checksize($fun3cog ); }
		if($nokegg) { $wsizeKEGG = 1; } else { $wsizeKEGG = checksize($fun3kegg); }
		if($nopfam) { $wsizePFAM = 1; } else { $wsizePFAM = checksize($fun3pfam); }
		$wsizeOPTDB = 1; # this will be 0 if any of the EXTDB files has no results or is missing
		if($opt_db) {
			open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n";
			while(<infile0>) {
				my($dbname,$extdb,$dblist)=split(/\t/,$_);
				my $dbdname="$resultpath/07.$projectname.fun3.$dbname";
				my $wsize=checksize($dbdname);
				if($wsize<1) { $wsizeOPTDB=0; }
				}
			close infile0;
			}	
             	if(($wsizeCOG>=1) && ($wsizeKEGG>=1) && ($wsizePFAM>=1) && ($wsizeOPTDB>=1) && (!$force_overwrite)) { print "Functional assignments already found, skipping step 7\n"; }
		else {
			my $scriptname="07.fun3assign.pl";
			if((!$nocog) || (!$nokegg) || (!$nopfam) || ($opt_db)) {
				print outfile3 "7\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "\n[",$currtime->pretty,"]: STEP7 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP7 -> FUNCTIONAL ASSIGNMENT: $scriptname\n"; print RESET;
				if($verbose) { print " (This will use fun3 algorithm to annotate putative functions for each ORF, from the homologues found in step 5)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir 0 $force_overwrite");
				if($ecode!=0)   {  error_out(7,$scriptname); }
				if($nocog)  { $wsizeCOG  = 1; } else { $wsizeCOG  = checksize($fun3cog ); }
				if($nokegg) { $wsizeKEGG = 1; } else { $wsizeKEGG = checksize($fun3kegg); }
				if($nopfam) { $wsizePFAM = 1; } else { $wsizePFAM = checksize($fun3pfam); }
				$wsizeOPTDB = 1; # this will be 0 if any of the EXTDB files has no results or is missing
				if($opt_db) {
					open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
					while(<infile0>) {
						my($dbname,$extdb,$dblist)=split(/\t/,$_);
						my $dbdname="$resultpath/07.$projectname.fun3.$dbname";
						my $wsize=checksize($dbdname);
						if($wsize<1) { $wsizeOPTDB=0; }
						}
					close infile0;
					}
				if(($wsizeCOG<1) && ($wsizeKEGG<1) && ($wsizePFAM<1) && ($wsizeOPTDB<1)) { error_out(7,$scriptname,"$fun3cog, $fun3kegg and $fun3pfam"); }
			}
		}
	close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP8: Blastx on the unannotated parts of the contigs
	
	if(($rpoint<=8) && (!$onlybins) && ((!$test) || ($test>=8))) {
		if($doublepass) {
			my $wsize=checksize($gff_file_blastx);
             		if(($wsize>=1) && (!$force_overwrite)) { print "Blastx file $gff_file_blastx already found, skipping step 8\n"; }
			else {		
				my $scriptname="08.blastx.pl";
				# print " DOUBLEPASS: Now starting blastx analysis\n";
				print outfile3 "8\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "\n[",$currtime->pretty,"]: STEP8 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP8 -> DOUBLEPASS, Blastx analysis: $scriptname\n"; print RESET;
				if($verbose) { print " (This will do many things: it will mask the parts of the contigs where an ORF has been found, and will run blastx in the remaining gaps, to identify possible genes missed in gene prediction. This is intended to be useful when dealing with eukaryotic or viral sequences, for which gene prediction is less accurate)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir $force_overwrite");
				if($ecode!=0)  { error_out(8,$scriptname); }
				my $wsize=checksize($gff_file_blastx);
				if($wsize<1)   { error_out(8,$scriptname,$gff_file_blastx); }
				}
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
		

    #-------------------------------- STEP9: Taxonomic annotation for the contigs (consensus of gene annotations)


	if(($rpoint<=9) && (!$onlybins) && ((!$test) || ($test>=9))) {
		my $wsize=checksize($alllog);
             	if(($wsize>=1) && (!$force_overwrite)) { print "Contig tax file $alllog already found, skipping step 9\n"; }
		else {		
			my $scriptname="09.summarycontigs3.pl";
			print outfile3 "9\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP9 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP9 -> CONTIG TAX ASSIGNMENT: $scriptname\n"; print RESET;
			if($verbose) { print " (This will produce a consensus taxonomic annotation for each contig, according to the annotations of their constituent ORFs)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(9,$scriptname); }
			my $wsize=checksize($alllog);
			if($wsize<1)         { error_out(9,$scriptname,$alllog); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
			
    #-------------------------------- STEP10: Mapping of reads onto contigs for abundance calculations
	
	if(($rpoint<=10) && ((!$test) || ($test>=10))) {
		my $ns;
		if($mode eq "sequential") { $ns = 1; } else { $ns = $numsamples; }
		my $wsize = checksize($mappingstat);
		if(($wsize == $ns) && (!$force_overwrite)) { print "Mapping file $mappingstat already found, skipping step 10\n"; }
		else {	
			my $scriptname="10.mapsamples.pl";
			print outfile3 "10\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP10 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP10 -> MAPPING READS: $scriptname\n"; print RESET;
			if($verbose) { print " (This will map reads back to the contigs using $mapper and count how many map to each ORF, to estimate their abundances)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir $force_overwrite");
			if($ecode!=0)           { error_out(10,$scriptname); }
			my $wsize = checksize($mappingstat);
			if($wsize!=$ns) { error_out(10,$scriptname,$mappingstat); }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
			
    #-------------------------------- STEP11: Count of taxa abundances
	
	if(($rpoint<=11) && (!$onlybins) && ((!$test) || ($test>=11))) {
		my $wsize=checksize($mcountfile);
             	if(($wsize>=2) && (!$force_overwrite)) { print "Abundance file $mcountfile already found, skipping step 11\n"; }	
		else {	
			my $scriptname="11.mcount.pl";
			print outfile3 "11\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP11 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP11 -> COUNTING TAX ABUNDANCES: $scriptname\n"; print RESET;
			if($verbose) { print " (This will count the abundances of each taxon, to infer the composition of the community)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(11,$scriptname); }
			my $wsize=checksize($mcountfile);
			if($wsize<2)         { error_out(11,$scriptname,$mcountfile); }
			close(outfile4); open(outfile4,">>$syslogfile");
		}		
	}
			
    #-------------------------------- STEP12: Count of function abundances
	
	if(($rpoint<=12) && (!$onlybins) && ((!$test) || ($test>=12))) {
		# Set wsize to 1 for a method if we were instructed not to use it,                                                          
                #  so that the step counts as completed even if those results are not present.
		my $cogfuncover = "$resultpath/12.$projectname.cog.funcover";
		my $keggfuncover= "$resultpath/12.$projectname.kegg.funcover";
		my($wsizeCOG,$wsizeKEGG);
		if($nocog)  { $wsizeCOG  = 1; } else { $wsizeCOG  = checksize($cogfuncover ); }
		if($nokegg) { $wsizeKEGG = 1; } else { $wsizeKEGG = checksize($keggfuncover); }
             	if(($wsizeCOG>=1) && ($wsizeKEGG>=1) && (!$force_overwrite)) { print "Function abundance files already found, skipping step 12\n"; }
		else {	
			my $scriptname="12.funcover.pl";
			if((!$nocog) || (!$nokegg) || ($opt_db)) {
			print outfile3 "12\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP12 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP12 -> COUNTING FUNCTION ABUNDANCES: $scriptname\n"; print RESET;
			if($verbose) { print " (This will count the abundance of each function, to produce a functional profile of the community)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)     { error_out(12,$scriptname); }
			if($nocog)  { $wsizeCOG  = 1; } else { $wsizeCOG  = checksize($cogfuncover ); }
			if($nokegg) { $wsizeKEGG = 1; } else { $wsizeKEGG = checksize($keggfuncover); }
			if(($wsizeCOG<1) && ($wsizeKEGG<1)) { error_out(12,$scriptname,"$cogfuncover and/or $keggfuncover"); }
			}
			close(outfile4); open(outfile4,">>$syslogfile");
		}
	}
			
    #-------------------------------- STEP13: Generation of the gene table
		
	if(($rpoint<=13) && (!$onlybins) && ((!$test) || ($test>=13))) {
		my $wsize = checksize($mergedfile);
             	if(($wsize>=2) && (!$force_overwrite)) { print "ORF table $mergedfile already found, skipping step 13\n"; }
		else {		
			my $scriptname="13.mergeannot2.pl";
			print outfile3 "13\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP13 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP13 -> CREATING GENE TABLE: $scriptname\n"; print RESET;
			if($verbose) { print " (This will create the gene table by merging all the information compiled in previous steps)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(13,$scriptname); }
			my $wsize=checksize($mergedfile);
			if($wsize<2)         { error_out(13,$scriptname,$mergedfile); }
			if($verbose) { print " (Now we already have a GENE TABLE)\n"; }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}
			
    #-------------------------------- STEP14: Running binning methods 		
	
	 if(!$nobins) {	    
	 	my $hayresults; 
		if(($rpoint<=14) && ((!$test) || ($test>=14)) && (!$extbins)) {
			if($verbose) { print " (Now we will start creating bins for separating individual organisms in the community)\n"; }
			my @binner=split(/\,/,$binners);
			foreach my $tbinner(@binner) { #-- Checking for results for all the specified binners
				my @binfiles;
				my $wsize=0;
				my $firstfile;
				my $dirbin="$interdir/binners/$tbinner";
				opendir(indir1,$dirbin);
				@binfiles=grep(/fasta$|fa$/,readdir indir1);
				closedir indir1;
				$firstfile="$dirbin/$binfiles[0]";
				$wsize=checksize($firstfile);
				if($wsize>=2) { $hayresults=1; last; }
			}
           	 	if(($hayresults) && (!$force_overwrite)) { print "Binning results for $binners already found, skipping step 14\n"; }
			else {
				my $scriptname="14.runbinning.pl";
				print outfile3 "14\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP14 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP14 -> BINNING: $scriptname\n"; print RESET;
				if($verbose) { print " (This will use binning programs for creating a set of bins)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
				if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
				}
		}
		
			
 
    #-------------------------------- STEP15: DAS Tool merging of binning results	
	
		if(($rpoint<=15) && ((!$test) || ($test>=15)) && (!$extbins)) {
			opendir(indir2,$binresultsdir);
			my @binfiles=grep(/fa/,readdir indir2);
			closedir indir2;
			my $wsize;
			if(scalar @binfiles) {
				my $firstfile="$binresultsdir/$binfiles[0]";
				$wsize=checksize($firstfile);
				}
			else { $wsize = 0; }
            	 	if(($wsize>=2) && (!$force_overwrite)) { print "DASTool results in $binresultsdir already found, skipping step 15\n"; }
			else {		
				my $scriptname="15.dastool.pl";
				print outfile3 "15\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP15 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP15 -> DAS_TOOL MERGING: $scriptname\n"; print RESET;
				if($verbose) { print " (This will use DASTool for creating a consensus between the sets of bins created in previous steps)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
				if($ecode!=0){ print RED; print "ERROR in STEP15-> $scriptname\n"; print RESET; }
				opendir(indir2,$binresultsdir) || warn "Can't open $binresultsdir directory, no DAStool results\n";
				my @binfiles=grep(/fa/,readdir indir2);
				closedir indir2;
				my $firstfile="$binresultsdir/$binfiles[0]";
				my ($wsize,$rest);
				my $wsize=checksize($firstfile);
				if($wsize<2) {
					print RED; print "WARNING: File $firstfile is empty!. DAStool did not generate results\n"; print RESET;
					$finaltrace.="DAS Tool abnormal termination: file $firstfile is empty. There are NO BINS!\n";
					$DAS_Tool_empty = 1;
					if($verbose) { print " (This will use DASTool for creating a consensus between the sets of bins created in previous steps)\n"; 	}
					}
				close(outfile4); open(outfile4,">>$syslogfile");
				}
		}
			
    #-------------------------------- STEP16: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if(($rpoint<=16) && (!$onlybins) && ((!$test) || ($test>=16))) {
			if(!$DAS_Tool_empty){
				opendir(indir3,$binresultsdir);
				my @binfiles=grep(/fa$/,readdir indir3);
				closedir indir3;
				my $all_ok = 1;
				foreach(@binfiles)
					{
					if(!checksize("$binresultsdir/$_.tax")) { $all_ok = 0; last; }
					}
            	 		if(($all_ok) && (!$force_overwrite))  { print "Bin annotation in $binresultsdir already found, skipping step 16\n"; }
				else {		
					my $scriptname="16.addtax2.pl";
					print outfile3 "16\t$scriptname\n";
					$currtime=timediff();
					print outfile4 "[",$currtime->pretty,"]: STEP16 -> $scriptname\n";
					print BLUE "[",$currtime->pretty,"]: STEP16 -> BIN TAX ASSIGNMENT: $scriptname\n"; print RESET;
					if($verbose) { print " (This will produce taxonomic assignments for each bin as the consensus of the annotations of each of their contigs)\n"; }
					my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
					if($ecode!=0){ print RED; print "Stopping in STEP16 -> $scriptname\n"; print RESET; die; }
					opendir(indir3,$binresultsdir);
					my @binfiles=grep(/fa$/,readdir indir3);
					closedir indir3;
					my $all_ok = 1;
					foreach(@binfiles)
						{
						if(!checksize("$binresultsdir/$_.tax")) { $all_ok = 0; last; }
						}
					if(!$all_ok) { print RED; print "Stopping in STEP16 -> $scriptname. Some bins were not annotated!\n"; print RESET; die; }
					}
				}
			else{ print RED; print "Skipping BIN TAX ASSIGNMENT: DAS_Tool did not predict bins.\n"; print RESET; }
			close(outfile4); open(outfile4,">>$syslogfile");
		}

			
    #-------------------------------- STEP17: Checking of bins for completeness and contamination (checkM)		
		my $new17=0;
		if(($rpoint<=17) && ((!$test) || ($test>=17))) {
			if(!$DAS_Tool_empty){
				my($filetocheck, $minlines);
				if($gtdbtk) { $filetocheck = $gtdbtkfile; $minlines = 2; } else { $filetocheck = $checkmfile; $minlines = 4; }
				my $wsize=checksize($filetocheck);
            	 		if(($wsize>=$minlines) && (!$force_overwrite))  { print "Results in $filetocheck already found, skipping step 17\n"; }
				else {
					$new17=1;	
					my $scriptname="17.checkbins.pl";
					print outfile3 "17\t$scriptname\n";
					$currtime=timediff();
					print outfile4 "[",$currtime->pretty,"]: STEP17 -> $scriptname\n";
					print BLUE "[",$currtime->pretty,"]: STEP17 -> CHECKING BINS: $scriptname\n"; print RESET;
					if($verbose) { print " (This step will use checkM for estimating completeness and contamination for each bin)\n"; }
					my $ecode = system("perl $scriptdir/$scriptname $projectdir");
					if($ecode!=0) { print RED; print "Stopping in STEP17 -> $scriptname\n"; print RESET; die; }
					my $binmethod="DAS";
					my $wsize=checksize($checkmfile);
					if($wsize<4) {
						print RED; print "Can't find $checkmfile\nStopping in STEP17 -> $scriptname\n"; print RESET; die; }
					if($gtdbtk) {
						my $wsize=checksize($checkmfile);
						if($wsize<2) {
							print RED; print "Can't find $gtdbtkfile\nStopping in STEP17 -> $scriptname\n"; print RESET; die; }
						}
					}
				}
			else { print RED; print"Skipping CHECKM: DAS_Tool did not predict bins.\n"; print RESET; }
		}

			
    #-------------------------------- STEP18: Make bin table		
	
		if(($rpoint<=18) && ((!$test) || ($test>=18))) {
			if(!$DAS_Tool_empty){
				my $wsize=checksize($bintable);
            	 		if(($wsize>=2) && (!$force_overwrite) && (!$new17)) { print "Bin table in $bintable already found, skipping step 18\n"; }
				else {		
					my $scriptname="18.getbins.pl";
					print outfile3 "18\t$scriptname\n";
					$currtime=timediff();
					print outfile4 "[",$currtime->pretty,"]: STEP18 -> $scriptname\n";
					print BLUE "[",$currtime->pretty,"]: STEP18 -> CREATING BIN TABLE: $scriptname\n"; print RESET;
					if($verbose) { print " (Now we will compile all previous information for producing a bin table)\n"; }
					my $ecode = system("perl $scriptdir/$scriptname $projectdir");
					if($ecode!=0){ print RED; print "Stopping in STEP18 -> $scriptname\n"; print RESET; die; }
					my $wsize=checksize($bintable);
					if($wsize<2) { print RED; print "Stopping in STEP18 -> $scriptname. File $bintable is empty!\n"; print RESET; die; }
					if($verbose) { print " (Now we have a BIN TABLE)\n"; }
					}
				}
			else { print RED; print "Skipping BIN TABLE CREATION: (You already know: DAS_Tool did not predict bins.)\n"; print RESET; }
			close(outfile4); open(outfile4,">>$syslogfile");
	 }
	}

    #-------------------------------- STEP19: Make contig table		

	if(($rpoint<=19) && (!$onlybins) && ((!$test) || ($test>=19))) {
		my $wsize=checksize($contigtable);
            	if(($wsize>=2) && (!$force_overwrite)) { print "Contig table in $contigtable already found, skipping step 19\n"; }
		else {		
			my $scriptname="19.getcontigs.pl";
			print outfile3 "19\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP19 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP19 -> CREATING CONTIG TABLE: $scriptname\n"; print RESET;
			if($verbose) { print " (In this step we will gather all information on contigs for producing a contig table. We do it now because we wanted to know the correspondence between contigs and bins)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { print RED; print "Stopping in STEP19 -> $scriptname\n"; print RESET; die; }
			my $wsize=checksize($contigtable);
			if($wsize<2)         { print RED; print "Stopping in STEP19 -> $scriptname. File $contigtable is empty!\n"; print RESET; die; }
			if($verbose) { print " (Now we have a CONTIG TABLE)\n"; }
			close(outfile4); open(outfile4,">>$syslogfile");
			}
	}

    #-------------------------------- STEP20: Pathways in bins          

	if(!$nobins) {	       
		if(($rpoint<=20) && (!$onlybins) && ((!$test) || ($test>=20))) {
			if((!$DAS_Tool_empty) && (!$nokegg)) {
				my $minpathfile="$resultpath/20.$projectname.kegg.pathways";
				my $wsize=checksize($minpathfile);
            			if(($wsize>2) && (!$force_overwrite))  { print "Pathways file $minpathfile already found, skipping step 20\n"; }
				else {		
					my $scriptname="20.minpath.pl";
					print outfile3 "20\t$scriptname\n";
					$currtime=timediff();
					print outfile4 "[",$currtime->pretty,"]: STEP20 -> $scriptname\n";
					print BLUE "[",$currtime->pretty,"]: STEP20 -> CREATING TABLE OF PATHWAYS IN BINS: $scriptname\n"; print RESET;
					if($verbose) { print " (In this step we will use MinPath and the gene content information for each bin to predict the possible presence of metabolic pathways in it)\n"; }
					my $ecode = system("perl $scriptdir/$scriptname $projectdir");
					if($ecode!=0){ print RED; print "Stopping in STEP20 -> $scriptname\n"; print RESET; die; }
					my $wsize=checksize($minpathfile);
					if($wsize<3) { print RED; print "Stopping in STEP20 -> $scriptname. File $minpathfile is empty!\n"; print RESET; die; }
					}
				}
			else{ print("Skipping MINPATH: DAS_Tool did not predict bins.\n") ; }
		}
		close(outfile4); open(outfile4,">>$syslogfile");
	}


    #-------------------------------- STEP21: Make stats		

	if($rpoint<=21) {
		my $statfile="$resultpath/21.$projectname.stats";
		my $wsize=checksize($statfile);
            	if(($wsize>2) && (!$force_overwrite)) { print "Statistics in $statfile already found, skipping step 21\n"; }
		else {		
			my $scriptname="21.stats.pl";
			print outfile3 "21\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP21 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP21 -> MAKING FINAL STATISTICS: $scriptname\n"; print RESET;
			if($verbose) { print " (Finally, we will produce final statistics on the run gathering information on genes, contigs and bins\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { print RED; print "Stopping in STEP21 -> $scriptname\n"; print RESET; die; }
			my $wsize=checksize($statfile);
			if($wsize<10)        { print RED; print "Stopping in STEP21 -> $scriptname. File $statfile is empty!\n"; print RESET; die; }
			}
	}

    #-------------------------------- END OF PIPELINE		

	print outfile3 "END\n";
	$currtime=timediff();
	print "\nDeleting temporary files in $tempdir\n";
	print outfile4 "\nDeleting temporary files in $tempdir\n";
	# system("rm -r $tempdir/*");
	if(-e "$datapath/megahit/final.contigs.fa") { system("rm -r $datapath/megahit/intermediate_contigs; rm $datapath/megahit/final.contigs.fa"); } 
	print outfile4 "\n[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
	print BLUE "[",$currtime->pretty,"]: FINISHED -> Have fun!\n"; print RESET;
	if($finaltrace) { print "\nWARNINGS:\n$finaltrace\n"; }
	print "For citation purposes, you can find a summary of methods in the file $methodsfile\n\n";
	if(!$onlybins) {
		print "You can analyze your results using the SQMtools R library (see https://github.com/jtamames/SqueezeMeta/wiki/Using-R-to-analyze-your-SQM-results)\n";
		if($mode eq "sequential") { print "(Please remember that sequential projects must be loaded indepently in SQMtools)\n"; }
	}
}


sub error_out {			#-- Catch abnormal terminations of pipeline steps
	outfile4->autoflush;		
	my($step,$scriptname,$file)=@_;
	if($file) { print RED; print "Stopping in STEP$step -> $scriptname. File $file is empty!\n"; print RESET; print outfile4 "Stopping in STEP$step -> $scriptname. File $file is empty!\n"; }
	else { print RED; print "Stopping in STEP$step -> $scriptname. Program finished abnormally\n"; print RESET; print outfile4 "Stopping in STEP$step -> $scriptname. Program finished abnormally\n"; }
	print outfile4 "_____________\n\nSystem information:\n";
	system("uname -a >> $projectdir/syslog");
	print outfile4 "_____________\n\nTree for the project:\n";
	system("tree -h -D $projectdir >> $projectdir/syslog");
	close outfile4;
	print RED; print "\n  If you don't know what went wrong or want further advice, please look for similar issues in https://github.com/jtamames/SqueezeMeta/issues\n";
	print "  Feel free to open a new issue if you don't find the answer there. Please add a brief description of the problem and upload the $projectdir/syslog file (zip it first)\n";
	print RESET;
	die;
}


sub writeconf {			#-- Create directories and files, write the SqueeeMeta_conf file
	
	my $projectdir=shift;
	my $scriptdir=shift;
	my %conf=@_;
	
		#-- Create directories for sequential mode (in other modes, the only project directory was already created)
	
	if ($mode=~/sequential$/i) {
		if((-d $projectdir) && (!$restart)) { print RED; print "Project name $projectdir already exists. Please remove it or change the project name\n"; print RESET; die; } 
		else { system("mkdir $projectdir"); }
		}

		#-- Write SqueezeMeta_conf.pl
	
	open(outfile5,">$projectdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't write in directory $projectdir. Out of space?\n"; print RESET; die; };
	open(infile2,"$scriptdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't open conf file $scriptdir/SqueezeMeta_conf.pl\n"; print RESET; die; };

	print outfile5 "\$version = \"$conf{version}\";\n";
	print outfile5 "\$mode = \"$conf{mode}\";\n";
	print outfile5 "\$date = \"",scalar localtime,"\";\n\n";
	print outfile5 "\$installpath = \"$installpath\";\n";
	# print outfile5 "\$userdir = \"$rawfastq\";\n";
	if($cleaning) { print outfile5 "\$userdir       = \"$projectdir/data/raw_fastq\";\n"; }
	else { print outfile5 "\$userdir       = \"$rawfastq\";\n";    }       

	while(<infile2>) {
		chomp;
		next if !$_;		
		if   ($_=~/^\$projectname/)     { print outfile5 "\$projectname = \"$conf{projectname}\";\n";         }
		elsif($_=~/^\$blocksize/)       { print outfile5 "\$blocksize       = $conf{blocksize};\n";           }
		elsif($_=~/^\$nodiamond/)       { print outfile5 "\$nodiamond       = $conf{nodiamond};\n";           }
		elsif($_=~/^\$fastnr/)          { print outfile5 "\$fastnr          = $conf{fastnr};\n";              }
		elsif($_=~/^\$singletons/)      { print outfile5 "\$singletons      = $conf{singletons};\n";          }
		elsif($_=~/^\$nocog/)           { print outfile5 "\$nocog           = $conf{nocog};\n";               }
		elsif($_=~/^\$nokegg/)          { print outfile5 "\$nokegg          = $conf{nokegg};\n";              }
		elsif($_=~/^\$nopfam/)          { print outfile5 "\$nopfam          = $conf{nopfam};\n";              }
		elsif($_=~/^\$euknofilter/)     { print outfile5 "\$euknofilter     = $conf{euknofilter};\n";         }
		elsif($_=~/^\$doublepass/)      { print outfile5 "\$doublepass      = $conf{doublepass};\n";          }
		elsif($_=~/^\$nobins/)          { print outfile5 "\$nobins          = $conf{nobins};\n";              }
		elsif($_=~/^\$onlybins/)        { print outfile5 "\$onlybins        = $conf{onlybins};\n";            }
		elsif($_=~/^\$gtdbtk[= ]/)      { print outfile5 "\$gtdbtk          = $conf{gtdbtk};\n";              }
		elsif($_=~/^\$binners/)         { print outfile5 "\$binners         = \"$conf{binners}\";\n";         }
		elsif($_=~/^\$norename/)        { print outfile5 "\$norename        = $conf{norename};\n";            }
		elsif($_=~/^\$preserve/)	{ print outfile5 "\$preserve        = $conf{preserve};\n";            }		
		elsif($_=~/^\$mapper/)          { print outfile5 "\$mapper          = \"$conf{mapper}\";\n";          }
		elsif($_=~/^\$mapping_options/) { print outfile5 "\$mapping_options = \"$conf{mapping_options}\";\n"; }
		elsif($_=~/^\$cleaning\b/)      { print outfile5 "\$cleaning        = $conf{cleaning};\n";            }
		elsif($_=~/^\$cleaningoptions/) { print outfile5 "\$cleaningoptions = \"$conf{cleaningoptions}\";\n"; }
		elsif($_=~/^\$gtdbtk_data_path/) {
			if($conf{gtdbtk_data_path}) { print outfile5 "\$gtdbtk_data_path = \"$conf{gtdbtk_data_path}\";\n"; }
			else { print outfile5 "$_\n"; }
			}
		else { print outfile5 "$_\n"; }
		if($consensus) { print outfile5 "\$consensus=$conf{consensus};\n"; }
		elsif($minion) { print outfile5 "\$consensus=0.2;\n"; }
	}
	close infile2; 

	print outfile5 "\n#-- Options\n\n";
	print outfile5 "\$numthreads         = $conf{numthreads};\n";
	print outfile5 "\$mincontiglen       = $conf{mincontiglen};\n";
	print outfile5 "\$assembler          = \"$conf{assembler}\";\n";
	print outfile5 "\$canumem            = $conf{canumem};\n";
	if($contigid)          { print outfile5 "\$contigid           = \"$conf{contigid}\";\n";          }
	if($assembler_options) { print outfile5 "\$assembler_options  = \"$conf{assembler_options}\";\n"; }
	if($extassembly)       { print outfile5 "\$extassembly        = \"$conf{extassembly}\";\n";       }
	if($extbins)           { print outfile5 "\$extbins            = \"$conf{extbins}\";\n";           }
	if($opt_db)            { print outfile5 "\$opt_db             = \"$conf{opt_db}\";\n";            }
	if($newtaxdb)          { print outfile5 "\$newtaxdb           = \"$conf{newtaxdb}\";\n";          }
	if($taxbinmode)        { print outfile5 "\$taxbinmode         = \"$taxbinmode\";\n";              }
	close outfile5;
	
	#--  Write progress and syslog
	
	open(outfile3,">$projectdir/progress") or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };  #-- An index indicating where are we and which parts of the method finished already. For the global process
	open(outfile4,">$projectdir/syslog")  or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };	#-- A log file for the global proccess
	$currtime=timediff();
	print outfile4 "Run started ",scalar localtime," in $mode mode\n";
	my $params = join(" ", @ARGV);
	# print outfile2 "$0 $params\n";
	print outfile4 "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-SÃ¡nchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n";
	print outfile4 "Run started for $projectname, ",scalar localtime,"\n";
	print outfile4 "Project: $conf{projectname}\n";
	print outfile4 "Map file: $conf{samples}\n";
	print outfile4 "Fastq directory: $conf{userdir}\n";
	print outfile4 "Command: $conf{commandline}\n"; 
	print outfile4 "[",$currtime->pretty,"]: STEP0 -> SqueezeMeta.pl\n";
	if(!$nocog) { print outfile4 " COGS;"; }
	if(!$nokegg) { print outfile4 " KEGG;"; }
	if(!$nopfam) { print outfile4 " PFAM;"; }
	if($opt_db) { print outfile4 " EXT_DB: $opt_db;"; }
	if($euknofilter) { print outfile4 " EUKNOFILTER;"; } 
	if($doublepass) { print outfile4 " DOUBLEPASS;"; }
	if($lowmem) { print outfile4 " LOW MEMORY MODE;"; }
	print outfile4 "\n";
	print outfile2 "[",$currtime->pretty,"]: STEP0 -> SqueezeMeta.pl\n";
	
	print "Now creating directories\n";
	open(infile2,"$scriptdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't open conf file $scriptdir/SqueezeMeta_conf.pl\n"; print RESET; die; };

	#-- Creation of directories
 
	print "Reading configuration from $projectdir/SqueezeMeta_conf.pl\n";
	do "$projectdir/SqueezeMeta_conf.pl";
	if(!$restart) {
		system ("mkdir $datapath");
 		system ("mkdir $resultpath");
		system ("mkdir $resultpath/bins");
 		system ("mkdir $tempdir");
 		system ("mkdir $datapath/raw_fastq"); 
 		system ("mkdir $extpath"); 
		system ("mkdir $interdir");
		system ("mkdir $interdir/binners");
		}
	open(outmet,">$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
	print outmet "Analysis done with SqueezeMeta v$version (Tamames & Puente-Sanchez 2019, Frontiers in Microbiology 9, 3349)\n";
	close outmet;
	open(outcreator,">$projectdir/creator.txt");
	print outcreator "SqueezeMeta v$version\n";
	close outcreator;
	if(!$nobins) {
		my $validbinners=join(",",keys %binscripts);
		my @binner=split(/\,/,$binners);
		foreach my $tbinner(@binner) { 
			if(!$binscripts{$tbinner}) { print RED; print "UNRECOGNIZED binner $tbinner (valid ones are $validbinners)\n"; print RESET; die; }
			}
		}

        #-- Adjusting parameters.pl file for custom identity and evalue

        # system("cp $equivfile $mappingfile");
	my $miniden=$conf{miniden};
	my $evalue=$conf{evalue};
        if((!$miniden) && (!$evalue)) { system("cp $scriptdir/parameters.pl $projectdir"); }
        else {
                open(outpar,">$projectdir/parameters.pl") || die "Cannot create new parameter file in $projectdir/parameters.pl\n";
                open(inpar,"$scriptdir/parameters.pl") || die "Cannot open parameter file in $scriptdir/parameters.pl\n";
                while(<inpar>) {
                        if($miniden && ($_=~/^\$miniden.*?\=(\d+)/)) {
                                $_=~s/$1/$miniden/;
                                print outpar $_;
                                }
                        elsif($evalue && ($_=~/^\$evalue.*?\=([^;]+)/)) {
                                $_=~s/$1/$evalue/;
                                print outpar $_;
                                }
                        else { print outpar $_; }
                        }
                close inpar;
                close outpar;
                }
		
		#--Creation of samples file

	open(infile0,$conf{samples}) || die;       
        open(outfile0,">$mappingfile") || die;
        while(<infile0>) {
                $_=~s/\r//g;
		if($mode ne "sequential") { print outfile0 $_; }
		else {							#-- In sequential mode just put the current sample in the new samples file
			my @k=split(/\t/,$_);
			if($k[0] eq $conf{projectname}) { print outfile0 $_; }
			}
                }
        close outfile0;
        close infile0;


	}


sub cleaning  {

	my $projectdir=shift;
	my $scriptdir=shift;
	my $thissample=shift;
	my %conf=@_;
	my $newuserdir="$projectdir/data/raw_fastq";
	 			
	my($par1files,$par2files)=0;
	my($par1name,$par2name);
	my %prepsamples;

	#-- NOW FILTERING READS
	
	print "  Running trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20) for quality filtering\n  Parameters: $cleaningoptions\n";
	my $trimmomatic_command;
	open(infile4,$equivfile) or do { print RED; print "Can't open samples file (-s) in $equivfile. Please check if that is the correct file, it is present tin that location, and you have reading permissions\n"; print RESET; die; };
	if($cleaning) {
		while(<infile4>) {
 			chomp;
 			next if(!$_ || ($_=~/^\#/));
			$_=~s/\r//g;			#-- Deleting \r in samples file for windows compatibility
			my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
			if($thissample and ($sample ne $thissample)) { next; } #-- If we wanted only one sample, skip the rest
			#-- Store all the input files coming from the same sample and "pair" (pair1 and pair2) into an array
			#-- This will fail if the input files from the same sample have a different ordering for pair1 and pair2
			#--   in the samples file, but this is very unlikely
			if(!exists($prepsamples{$sample}{$iden})) { $prepsamples{$sample}{$iden} = [$file];}
			else { push @{ $prepsamples{$sample}{$iden}} , $file; }
			}
		close infile4;
		foreach my $ts(sort keys %prepsamples) { 
			print "  Working with sample $ts\n";
			if(exists($prepsamples{$ts}{pair2}) and (scalar(@{ $prepsamples{$ts}{pair1} }) != scalar(@{ $prepsamples{$ts}{pair2} }))) {
				die("Different number of pair1 and pair2 files for sample $ts\n");
				}
			my @pair1files = @{ $prepsamples{$ts}{pair1} };
			my ($par1name, $par2name, $trimmedpar1name, $trimmedpar2name);
			for my $i (0 ..  $#pair1files) {
				#-- Go through all the input files for this sample
				$par1name        = "$rawfastq/".$prepsamples{$ts}{pair1}[$i];
				$trimmedpar1name = "$newuserdir/".$prepsamples{$ts}{pair1}[$i];
				$par2name        = "";
				$trimmedpar2name = "";
				if(exists($prepsamples{$ts}{pair2})) {
					$par2name        = "$rawfastq/".$prepsamples{$ts}{pair2}[$i];
					$trimmedpar2name = "$newuserdir/".$prepsamples{$ts}{pair2}[$i];
					}
				if(-e $par2name) { $trimmomatic_command="$trimmomatic_soft PE -threads $numthreads -phred33 $par1name $par2name $trimmedpar1name $trimmedpar1name.removed $trimmedpar2name $trimmedpar2name.removed $cleaningoptions > /dev/null 2>&1"; }
				else { $trimmomatic_command="$trimmomatic_soft SE -threads $numthreads -phred33 $par1name $trimmedpar1name $cleaningoptions > /dev/null 2>&1"; }
				print outfile4 "Running trimmomatic: $trimmomatic_command";
				my $ecode = system $trimmomatic_command;
				if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
				}
			}
		print outmet "Quality filtering was done using Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20)\n";
		}
	}


sub checksize {
	my $tfile=shift;
	my $wsize;
	if(-e $tfile) {
		$wsize=qx(grep -cv "^#" $tfile); # this excludes comments!
		}
	else { $wsize=0; }
	return $wsize;
	}
