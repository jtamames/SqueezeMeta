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

my $longtrace=0;    #-- Reports an explanation msg for each of the steps

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
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

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

our $pwd=cwd();

our($nodiamond,$binners,$nocog,$nokegg,$nopfam,$singletons,$euknofilter,$opt_db,$nobins,$nomaxbin,$nometabat,$empty,$lowmem,$minion,,$consensus,$doublepass)="0";
our($numsamples,$numthreads,$canumem,$mode,$mincontiglen,$contigid,$assembler,$extassembly,$mapper,$projectdir,$mapping_options,$projectname,$project,$equivfile,$rawfastq,$blocksize,$evalue,$miniden,$assembler_options,$cleaning,$cleaningoptions,$ver,$hel,$methodsfile,$test,$rename);
our($binresultsdir,$databasepath,$extdatapath,$softdir,$datapath,$resultpath,$extpath,$tempdir,$interdir,$mappingfile,$contigsfna,$gff_file_blastx,$contigslen,$mcountfile,$checkmfile,$rnafile,$gff_file,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$mapcountfile,$contigcov,$contigtable,$mergedfile,$bintax,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bwa_soft,$minimap2_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft,$canu_soft,$trimmomatic_soft,$dastool_soft);
our(%bindirs,%dasdir,%binscripts);  

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
   
 Filtering: 
   --cleaning: Filters with Trimmomatic (Default: No)
   -cleaning_options [options]: Options for Trimmomatic (Default:LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30)
   
 Assembly: 
   -a: assembler <megahit,spades,rnaspades,canu, flye> (Default: megahit)
   -assembly_options [options]: Extra options to be passed when calling the mapper
   -c|-contiglen <size>: Minimum length of contigs (Default: 200)
   -extassembly <file>: External assembly, file containing a fasta file of contigs (overrides all assembly steps).
   --sg|--singletons: Add unassembled reads to the contig file, as if they were contigs  
   -contigid <string>: Nomenclature for contigs (Default: assembler�s name)
   --norename: Don't rename contigs (Use at your own risk, characters like '_' in contig names will make it crash)
   
 Mapping: 
   -map: mapping software <bowtie, bwa, minimap2-ont, minimap2-pb, minimap2-sr> (Default: bowtie) 
   -mapping_options [options]: Extra options to be passed when calling the mapper

 Annotation:  
   --nodiamond: Check if Diamond results are already in place, and just in that case skips the Diamond run (Default: no)
   --nocog: Skip COG assignment (Default: no)
   --nokegg: Skip KEGG assignment (Default: no)
   --nopfam: Skip Pfam assignment  (Default: no)
   --euk: Drop identity filters for eukaryotic annotation  (Default: no)
   -consensus <value>: Minimum percentage of genes for a taxon needed for contig consensus (Default: 50)
   --D|--doublepass: First pass looking for genes using gene prediction, second pass using BlastX  (Default: no)
   -extdb <database file>: List of user-provided databases
   -b|-block-size <block size>: block size for diamond against the nr database (Default: 8)
   
 Binning:
   --nobins: Skip all binning  (Default: no). Overrides --binners 
   --binners: Specify binning programs to be used (available: maxbin, metabat)  (Default: maxbin,metabat)
   
 Performance:
   -t <threads>: Number of threads (Default: 12)
   -canumem <mem>: memory for canu in Gb (Default: 32)
   --lowmem: run on less than 16Gb of memory (Default:no)

 Other:
   --minion: Run on MinION reads (assembler: canu; mapper: minimap2-ont; consensus: 20) (Default: no)
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
		     "nocog" => \$nocog,   
		     "nokegg" => \$nokegg,   
		     "nopfam" => \$nopfam,  
		     "euk" => \$euknofilter,
		     "extdb=s" => \$opt_db, 
		     "nobins" => \$nobins,   
		     "binners=s" => \$binners,   
		     "D|doublepass" => \$doublepass, 
		     "b|block_size=i" => \$blocksize,
		     "e|evalue=f" => \$evalue,   
		     "minidentity=f" => \$miniden,   
		     "assembly_options=s" => \$assembler_options,
		     "norename" => \$rename,
		     "cleaning" => \$cleaning,
		     "cleaning_options=s" => \$cleaningoptions,
		     "mapping_options=s" => \$mapping_options,
                     "minion" => \$minion,
		     "test=i" => \$test,
		     "empty" => \$empty,
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
if(!$euknofilter) { $euknofilter=0; }
if(!$doublepass) { $doublepass=0; }
if(!$nobins) { $nobins=0; }
if(!$binners) { $binners="maxbin,metabat2"; }
if(!$nomaxbin) { $nomaxbin=0; }
if(!$nometabat) { $nometabat=0; }
if(!$rename) { $rename=0; }
if(!$cleaningoptions) { $cleaningoptions="LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30"; }
if(!$cleaning) { $cleaning=0; $cleaningoptions=""; } 
if($consensus) { $consensus/=100; }

$mode=~tr/A-Z/a-z/;
if($opt_db) { $opt_db = abs_path($opt_db); }

#-- Override settings if running on lowmem or MinION mode.
if($lowmem) { $blocksize=3; $canumem=15; }

if($minion) { $assembler="canu"; $mapper="minimap2-ont"; }

#-- Check if we have all the needed options


my($dietext,$finaltrace);
if($ver) { print "$version\n"; exit; }
if($hel) { die "$helptext\n"; } 
if(!$rawfastq) { $dietext.="MISSING ARGUMENT: -f|-seq: Fastq read files' directory\n"; }
if(!$equivfile) { $dietext.="MISSING ARGUMENT: -s|-samples: Samples file\n"; }
if(!$mode) { $dietext.="MISSING ARGUMENT: -m: Run mode (sequential, coassembly, merged)\n"; }
if(($mode!~/sequential$/i) && (!$projectdir)) { $dietext.="MISSING ARGUMENT: -p: Project name\n"; }
if(($mode=~/sequential$/i) && ($projectdir)) { $dietext.="Please DO NOT specify project name in sequential mode. The name will be read from the samples in the samples file $equivfile\n"; }
if($mode!~/sequential|coassembly|merged|seqmerge/i) { $dietext.="UNRECOGNIZED mode $mode (valid ones are sequential, coassembly, merged or seqmerge\n"; }
if($mapper!~/bowtie|bwa|minimap2-ont|minimap2-pb|minimap2-sr/i) { $dietext.="UNRECOGNIZED mapper $mapper (valid ones are bowtie, bwa, minimap2-ont, minimap2-pb or minimap2-sr\n"; }
if($assembler!~/megahit|spades|rnaspades|canu|flye/i) { $dietext.="UNRECOGNIZED assembler $assembler (valid ones are megahit, spades, canu or flye)\n"; }
if(($assembler=~/flye/i) && ($mode=~/merge/i)) { $dietext.="Invalid combination of mode and assembler\n (We are sorry for this, the low number of contigs provided by Flye prevents minimus2 needed in $mode mode to work correctly\n Please use coassembly, or a different assembler)\n"; }
if($rawfastq=~/^\//) {} else { $rawfastq=abs_path($rawfastq); }
if($dietext) { print BOLD "$helpshort"; print RESET; print RED; print "$dietext"; print RESET;  exit; }

$projectdir = abs_path($projectdir);
$projectname = (split '/', $projectdir)[-1];
my $syslogfile="$projectdir/syslog";

	#-- Check that everything is correct in the samples file

my %pairsample;
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
}
close infile1;
foreach my $chsam(keys %pairsample) { 
	if($pairsample{$chsam}!~/pair1/) { print RED; print "Sample $chsam has not pair1 in the samples file. Please check\n"; print RESET;  exit; }
	}

my $currtime=timediff();
print BOLD "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 9, 3349 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;
if($test) { print GREEN "Running in test mode! I will stop after step $test\n\n"; print RESET; }
print "Run started ",scalar localtime," in $mode mode\n";


#--------------------------------------------------------------------------------------------------
#----------------------------------- SEQUENTIAL MODE ----------------------------------------------

if($mode=~/sequential/i) { 

	my(%allsamples,%ident,%noassembly);
	my($sample,$file,$iden,$mapreq);
	tie %allsamples,"Tie::IxHash";

	#-- Reading the sample file given by the -s option, to locate the sample files

	print "Reading samples from $equivfile\n";
	open(infile1,$equivfile) or do { print RED; print "Can't open samples file (-s) in $equivfile. Please check if that is the correct file, it is present in that location, and you have reading permissions\n"; print RESET;  die; };
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		$_=~s/\r//g;
		($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		$allsamples{$sample}{$file}=1;
		$ident{$sample}{$file}=$iden;
		if($mapreq && ($mapreq=~/noassembly/i)) { $noassembly{$file}=1; }   #-- Files marked for no assembly (but they will be mapped)
	}
	close infile1;

	my @nmg=keys %allsamples;
	$numsamples=$#nmg+1;
	print "$numsamples metagenomes found: @nmg";
	print "\n\n";

	open(outfile2,">$pwd/syslog") || do { print RED; print "Can't write in directory $pwd\n"; print RESET; die; }; 		 #-- A log file for the global proccess
	print outfile2 "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n";
	print outfile2 "Run started ",scalar localtime," in SEQUENTIAL mode (it will proccess all metagenomes sequentially)\n";
	print outfile2 "Command: $commandline\n"; 
	print outfile2 "Options: threads=$numthreads; contiglen=$mincontiglen; assembler=$assembler; sample file=$equivfile; raw fastq=$rawfastq\n";

	#-- Now we start processing each sample individually

	foreach my $thissample(keys %allsamples) {

		#-- We start creating directories, progress and log files

		$projectname=$thissample;
		my $projectdir="$pwd/$thissample";
		if (-d $projectdir) { print RED; print "Project name $projectdir already exists. Please remove it or change the project name\n"; print RESET; die; } else { system("mkdir $projectdir"); }
		print "Working with $thissample\n";
	
		open(outfile3,">$projectdir/progress") or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };  #-- An index indicating where are we and which parts of the method finished already. For the global process
		open(outfile4,">$projectdir/syslog")  or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; }; 	#-- A log file for the global proccess
		$currtime=timediff();
		print outfile4 "Run started ",scalar localtime," in SEQUENTIAL mode (it will proccess all metagenomes sequentially)\n";
		print "Run started ",scalar localtime," in SEQUENTIAL mode\n";
		my $params = join(" ", @ARGV);
		# print outfile2 "$0 $params\n";
		print outfile4 "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n";
		print outfile4 "Run started for $thissample, ",scalar localtime,"\n";
		print outfile4 "Project: $projectname\n";
		print outfile4 "Map file: $equivfile\n";
		print outfile4 "Fastq directory: $rawfastq\n";
		print outfile4 "Command: $commandline\n"; 
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
	
		#-- Creation of the new configuration file for this sample
	
		open(outfile5,">$projectdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't write in directory $projectdir. Out of space?\n"; print RESET; die; };

		print outfile5 "\$version = \"$version\";\n";
		print outfile5 "\$mode = \"$mode\";\n";
		print outfile5 "\$date = \"",scalar localtime,"\";\n\n";
		print outfile5 "\$installpath = \"$installpath\";\n";

		while(<infile2>) {
			chomp;
			next if !$_;
			if   ($_=~/^\$projectname/)     { print outfile5 "\$projectname = \"$projectname\";\n";         }
			elsif($_=~/^\$blocksize/)       { print outfile5 "\$blocksize       = $blocksize;\n";           }
			elsif($_=~/^\$nodiamond/)       { print outfile5 "\$nodiamond       = $nodiamond;\n";           }
			elsif($_=~/^\$singletons/)      { print outfile5 "\$singletons      = $singletons;\n";          }
			elsif($_=~/^\$nocog/)           { print outfile5 "\$nocog           = $nocog;\n";               }
			elsif($_=~/^\$nokegg/)          { print outfile5 "\$nokegg          = $nokegg;\n";              }
			elsif($_=~/^\$nopfam/)          { print outfile5 "\$nopfam          = $nopfam;\n";              }
			elsif($_=~/^\$euknofilter/)     { print outfile5 "\$euknofilter     = $euknofilter;\n";         }
			elsif($_=~/^\$doublepass/)      { print outfile5 "\$doublepass      = $doublepass;\n";          }
			elsif($_=~/^\$nobins/)          { print outfile5 "\$nobins          = $nobins;\n";              }
			elsif($_=~/^\$binners/)         { print outfile5 "\$binners         = \"$binners\";\n";         }
			elsif($_=~/^\$rename/)          { print outfile5 "\$rename          = $rename;\n";           }
			elsif($_=~/^\$mapper/)          { print outfile5 "\$mapper          = \"$mapper\";\n";          }
			elsif($_=~/^\$mapping_options/) { print outfile5 "\$mapping_options = \"$mapping_options\";\n"; }
			elsif($_=~/^\$cleaning\b/)      { print outfile5 "\$cleaning        = $cleaning;\n";            }
			elsif($_=~/^\$cleaningoptions/) { print outfile5 "\$cleaningoptions = \"$cleaningoptions\";\n"; }
			else { print outfile5 "$_\n"; }
			if($consensus) { print outfile5 "\$consensus=$consensus;\n"; }
			elsif($minion) { print outfile5 "\$consensus=0.2;\n"; }
		}
	 	close infile2; 

		print outfile5 "\n#-- Options\n\n";
		print outfile5 "\$numthreads         = $numthreads;\n";
		print outfile5 "\$mincontiglen       = $mincontiglen;\n";
		print outfile5 "\$assembler          = \"$assembler\";\n";
		print outfile5 "\$canumem            = $canumem;\n";
		if($contigid) { print outfile5 "\$contigid            = \"$contigid\";\n"; }
		if($assembler_options) { print outfile5 "\$assembler_options  = \"$assembler_options\";\n"; }
		if($extassembly)       { print outfile5 "\$extassembly        = \"$extassembly\";\n";       }
		if($opt_db)            { print outfile5 "\$opt_db             = \"$opt_db\";\n";            }
		close outfile5;

		#-- Creation of directories
 
		print "Reading configuration from $projectdir/SqueezeMeta_conf.pl\n";
		do "$projectdir/SqueezeMeta_conf.pl";
		system ("mkdir $datapath");
 		system ("mkdir $resultpath");
 		system ("mkdir $tempdir");
 		system ("mkdir $datapath/raw_fastq"); 
 		system ("mkdir $extpath"); 
		system ("mkdir $interdir");
		system ("mkdir $interdir/binners");
		open(outmet,">$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
		print outmet "Analysis done with SqueezeMeta v$version (Tamames & Puente-Sanchez 2019, Frontiers in Microbiology 9, 3349)\n";
		close outmet;
		if(!$nobins) {
			my $validbinners=join(",",keys %binscripts);
			my @binner=split(/\,/,$binners);
			foreach my $tbinner(@binner) { 
				if(!$binscripts{$tbinner}) { print RED; print "UNRECOGNIZED binner $tbinner (valid ones are $validbinners)\n"; print RESET; die; }
				}
			}

	
		#-- Linkage of files to put them into our data directories
	
		print "Now linking read files\n";
 		foreach my $file(sort keys %{ $allsamples{$thissample} }) {
 			 if(-e "$rawfastq/$file") { 
  				my $tufile="$datapath/raw_fastq/$file";
 				# system("cp $rawfastq/$file $tufile");
 				 system("ln -s $rawfastq/$file $tufile");
 				# if($tufile!~/\.gz$/) { system("gzip $datapath/raw_fastq/$file"); }
			}
 			 else { print RED; print "Can't find read file $file (Sample $sample). Please check if it exists\n"; print RESET; die; }

	open(infile0,$equivfile) || die;	#-- Deleting \r in samples file for windows compatibility
	open(outfile0,">$mappingfile") || die;
	while(<infile0>) {
		$_=~s/\r//g;
		($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		if($sample eq $thissample) { print outfile0 $_; } # Generate a samples file containing only fastqs from this sample.
		}
	close outfile0;
	close infile0;

	#-- Adjusting parameters.pl file for custom identity and evalue

	# system("cp $equivfile $mappingfile");
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
	
		}

		#-- Preparing the files for the assembly, merging all files corresponding to each pair
		if(!$empty) {		
 			my($par1files,$par2files)=0;
			my($par1name,$par2name);
			my %prepsamples;
 			
			#-- NOW FILTERING READS
			
			my $trimmomatic_command;
			open(infile4,$equivfile) or do { print RED; print "Can't open samples file (-s) in $equivfile. Please check if that is the correct file, it is present tin that location, and you have reading permissions\n"; print RESET; die; };
			if($cleaning) {
				while(<infile4>) {
 					chomp;
 					next if(!$_ || ($_=~/^\#/));
					$_=~s/\r//g;			#-- Deleting \r in samples file for windows compatibility
					my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
					if($sample eq $thissample) { $prepsamples{$sample}{$iden}=$file; }
					}
				close infile4;
				foreach my $ts(sort keys %prepsamples) {
					my $par1name="$datapath/raw_fastq/".$prepsamples{$ts}{pair1};
					my $par2name="$datapath/raw_fastq/".$prepsamples{$ts}{pair2};
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
						print "  Running trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20) for quality filtering\n  Parameters: $cleaningoptions\n";
						print outfile4 "Running trimmomatic: $trimmomatic_command";
						my $ecode = system $trimmomatic_command;
						if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
						print outmet "Quality filtering was done using Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20)\n";
						}
					}
				}

			
			print "Now preparing files\n";
 			my($ca1,$ca2)="";
 			foreach my $afiles(sort keys %{ $ident{$thissample} }) {
 				next if($noassembly{$afiles});
  				my $gzfiles=$afiles;
  				#if($gzfiles!~/gz$/) { $gzfiles.=".gz"; }
 				if($ident{$thissample}{$afiles} eq "pair1") { $ca1.="$datapath/raw_fastq/$gzfiles "; $par1files++; } 
				else { $ca2.="$datapath/raw_fastq/$gzfiles "; $par2files++; } 
				if($gzfiles=~/gz$/) { $par1name="$datapath/raw_fastq/par1.fastq.gz"; $par2name="$datapath/raw_fastq/par2.fastq.gz"; }  # Fixed bug 30/10/2018 JT
				else { $par1name="$datapath/raw_fastq/par1.fastq"; $par2name="$datapath/raw_fastq/par2.fastq"; }
			}
			if(!$par1files) { print RED; print "There must be at least one 'pair1' sequence file in your samples file $mappingfile, and there is none!\n"; print RESET; die;  }
			if($par1files>1) { 
				my $command="cat $ca1 > $par1name"; 
				#print "$command\n"; 
				system($command); 
			} 
			else {
				#my $command="cp $ca1 $par1name"; 
				my $command="ln -s $ca1 $par1name";
				#print "$command\n";
				system($command);
 
			}
 			if($par2files>1) { system("cat $ca2 > $par2name"); } 
			elsif ($par2files==1) { system("ln -s $ca2 $par2name"); }    #-- Support for single reads
	
			#else { system("cp $ca2 $par2name"); }
			
		
			
			#-- CALL TO THE STANDARD PIPELINE
		
			pipeline(); 
		}
		else { print "Directory structure and conf files created for $thissample\n"; }
		
		
 		close outfile4;		#-- Closing log file for the sample
 		close outfile3;		#-- Closing progress file for the sample
							       
	}    #-- End of this sample, go for the next

}            #-- END


#-----------------------------------------------------------------------------------------------------
#----------------------------------- COASSEMBLY AND MERGED MODES -------------------------------------

else {

	if (-d $projectdir) {
		print RED;
		print "Project name $projectdir already exists. Please remove it or change project name\n";
		print RESET; 
		die;
	}else{
		my $ecode = system("mkdir $projectdir");
		if($ecode!=0){
			print RED;
			print "Can't create project directory at $projectdir\n";
			print RESET;
			die;
		}
	}
		
	#-- We start creating directories, progress and log files

	open(outfile3,">$projectdir/progress") or do { print RED; print "Can't write in $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; }; #-- Un indice que indica en que punto estamos (que procedimientos han terminado)
	open(outfile4,">$syslogfile") or do { print RED; print "Can't write in $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };
	my $params = join(" ", @ARGV);
	print outfile4 "\nSqueezeMeta v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n";
	print outfile4 "Run started ",scalar localtime," in $mode mode\n";
	print outfile4 "Command: $commandline\n"; 
	print outfile4 "Project: $projectname\n";
	print outfile4 "Map file: $equivfile\n";
	print outfile4 "Fastq directory: $rawfastq\n";
	print outfile4 "Options: threads=$numthreads; contiglen=$mincontiglen; assembler=$assembler; mapper=$mapper;";
	if(!$nocog) { print outfile4 " COGS;"; }
	if(!$nokegg) { print outfile4 " KEGG;"; }
	if(!$nopfam) { print outfile4 " PFAM;"; }
	if($opt_db) { print outfile4 " EXT_DB: $opt_db;"; }
	if($euknofilter) { print outfile4 " EUKNOFILTER;"; } 
	if($doublepass) { print outfile4 " DOUBLEPASS;"; }
	if($lowmem) { print outfile4 " LOW MEMORY;"; }
	print outfile4 "\n";
	print outfile4 "[",$currtime->pretty,"]: STEP0 -> SqueezeMeta.pl\n";

	print "Now creating directories\n";
	
	#-- Creation of the new configuration file for this sample
	open(infile3,"$scriptdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't open $scriptdir/SqueezeMeta_conf.pl\n"; print RESET; die; };
	open(outfile6,">$projectdir/SqueezeMeta_conf.pl") or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };

	print outfile6 "\$version = \"$version\";\n";
	print outfile6 "\$mode = \"$mode\";\n";
	print outfile6 "\$date = \"",scalar localtime,"\";\n\n";
	print outfile6 "\$installpath = \"$installpath\";\n";
	while(<infile3>) {
		if   ($_=~/^\$projectname/)               { print outfile6 "\$projectname = \"$projectname\";\n";                     }
		elsif($_=~/^\$blocksize/)                 { print outfile6 "\$blocksize       = $blocksize;\n";                       }
		elsif($_=~/^\$singletons/)                { print outfile6 "\$singletons      = $singletons;\n";                      }
		elsif($_=~/^\$nocog/)                     { print outfile6 "\$nocog           = $nocog;\n";                           }
	        elsif($_=~/^\$nodiamond/)                 { print outfile6 "\$nodiamond       = $nodiamond;\n";                       }
		elsif($_=~/^\$nokegg/)                    { print outfile6 "\$nokegg          = $nokegg;\n";                          }
		elsif($_=~/^\$nopfam/)                    { print outfile6 "\$nopfam          = $nopfam;\n";                          }
		elsif($_=~/^\$euknofilter/)               { print outfile6 "\$euknofilter     = $euknofilter;\n";                     }
		elsif($_=~/^\$nobins/)                    { print outfile6 "\$nobins          = $nobins;\n";                          }
		elsif($_=~/^\$binners/)                   { print outfile6 "\$binners         = \"$binners\";\n";                     }
		elsif($_=~/^\$nometabat/)                 { print outfile6 "\$nometabat       = $nometabat;\n";                       }
		elsif($_=~/^\$rename/)      		  { print outfile5 "\$rename          = $rename;\n";                          }
		elsif($_=~/^\$doublepass/)                { print outfile6 "\$doublepass      = $doublepass;\n";                      }
		elsif($_=~/^\$mapper/)                    { print outfile6 "\$mapper          = \"$mapper\";\n";                      }
		elsif($_=~/^\$mapping_options/)		  { print outfile6 "\$mapping_options = \"$mapping_options\";\n";             }
		elsif($_=~/^\$cleaning\b/)                { print outfile6 "\$cleaning        = $cleaning;\n";                        }
		elsif($_=~/^\$cleaningoptions/)           { print outfile6 "\$cleaningoptions = \"$cleaningoptions\";\n";             }
		elsif($_=~/^\%bindirs/) { print outfile6 "\%bindirs = (\"metabat2\",\"\$resultpath/metabat2\",\"maxbin\",\"\$resultpath/maxbin\");\n"; }
		else { print outfile6 $_; }
		if($consensus) { print outfile6 "\$consensus=$consensus;\n"; }
		elsif($minion) { print outfile6 "\$consensus=0.2;\n"; }
	 }
	close infile3;

	print outfile6 "\n#-- Options\n\n";
	print outfile6 "\$numthreads         = $numthreads;\n";
	print outfile6 "\$mincontiglen       = $mincontiglen;\n";
        if($contigid) { print outfile5 "\$contigid            = \"$contigid\";\n"; }
	print outfile6 "\$assembler          = \"$assembler\";\n";
	print outfile6 "\$canumem            = $canumem;\n";
	if($assembler_options) { print outfile6 "\$assembler_options  = \"$assembler_options\";\n"; }
	if($extassembly)       { print outfile6 "\$extassembly        = \"$extassembly\";\n";       }
	if($opt_db)            { print outfile6 "\$opt_db             = \"$opt_db\";\n";            }
	close outfile6;

	print "Reading configuration from $projectdir/SqueezeMeta_conf.pl\n";
	do "$projectdir/SqueezeMeta_conf.pl" or do { print RED; print "Can't write in directory $projectdir. Wrong permissions, or out of space?\n"; print RESET; die; };
	print "Reading samples from $mappingfile\n";

	open(outmet,">$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
	print outmet "Analysis done with SqueezeMeta v$version (Tamames & Puente-Sanchez 2019, Frontiers in Microbiology 9, 3349)\n";
	close outmet;
	if(!$nobins) {
		my $validbinners=join(",",keys %binscripts);
		my @binner=split(/\,/,$binners);
		foreach my $tbinner(@binner) { 
			if(!$binscripts{$tbinner}) { print RED; print "UNRECOGNIZED binner $tbinner (valid ones are $validbinners)\n"; print RESET; die; }
			}
		}

	#-- Creation of directories
	
	system ("mkdir $datapath");
	system ("mkdir $resultpath");
	system ("mkdir $tempdir");
	system ("mkdir $datapath/raw_fastq"); 
 	system ("mkdir $extpath"); 
	system ("mkdir $interdir");
	system ("mkdir $interdir/binners");

	#--Creation of samples file

        open(infile0,$equivfile) || die;        #-- Deleting \r in samples file for windows compatibility
        open(outfile0,">$mappingfile") || die;
        while(<infile0>) {
                $_=~s/\r//g;
                print outfile0 $_;
                }
        close outfile0;
        close infile0;

        #-- Adjusting parameters.pl file for custom identity and evalue

        # system("cp $equivfile $mappingfile");
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

	if(!$empty) {
 
		#-- Preparing the files for the assembly
	   
		moving();
	
		#-- CALL TO THE STANDARD PIPELINE
	
		pipeline();
		}
	else { die "  Directory structure and conf files created. Exiting\n"; }
	close outfile4;  #-- Closing log file for the sample
	close outfile3;	  #-- Closing progress file for the sample

}                        #------ END


#----------------------------PREPARING FILES FOR MERGING AND COASSEMBLY MODES-------------------------------------------

sub moving {
	
	#-- Reading samples from the file specified with -s option
	
	my(%allsamples,%ident,%prepsamples,%noassembly);


	open(infile4,$equivfile) or do { print RED; print "Can't open samples file (-s) in $equivfile. Please check if that is the correct file, it is present tin that location, and you have reading permissions\n"; print RESET; die; };
	while(<infile4>) {
 		chomp;
 		next if(!$_ || ($_=~/^\#/));
		$_=~s/\r//g;			#-- Deleting \r in samples file for windows compatibility
		my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
		$allsamples{$sample}=1;
		$ident{$file}=$iden;
		if(($mapreq) && ($mapreq=~/noassembly/i)) { $noassembly{$file}=1; }    #-- Files flagged for no assembly (but they will be mapped)
		if(-e "$rawfastq/$file") { 
			my $tufile="$datapath/raw_fastq/$file";
			 system("cp $rawfastq/$file $tufile");
			# system("ln -s $rawfastq/$file $tufile");
			print outfile4 "Coppying files: cp $rawfastq/$file $tufile\n";
			# if($tufile!~/\.gz$/) { system("gzip $datapath/raw_fastq/$file"); }
		}
		else { print RED; print "Can't find read file $file (Sample $sample)\n"; print RESET; die; }
	}
	close infile4;

	my $trimmomatic_command;
	open(infile4,$equivfile) or do { print RED; print "Can't open samples file (-s) in $equivfile. Please check if that is the correct file, it is present tin that location, and you have reading permissions\n"; print RESET; die; };
	if($cleaning) {
		while(<infile4>) {
 			chomp;
 			next if(!$_ || ($_=~/^\#/));
			$_=~s/\r//g;			#-- Deleting \r in samples file for windows compatibility
			my ($sample,$file,$iden,$mapreq)=split(/\t/,$_);
			$prepsamples{$sample}{$iden}=$file;
			}
		close infile4;
		foreach my $ts(sort keys %prepsamples) {
			my $par1name="$datapath/raw_fastq/".$prepsamples{$ts}{pair1};
			my $par2name="$datapath/raw_fastq/".$prepsamples{$ts}{pair2};
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
				print "  Running trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20) for quality filtering\n  Parameters: $cleaningoptions\n";
				print outfile4 "Running trimmomatic: $trimmomatic_command";
				my $ecode = system $trimmomatic_command;
				if($ecode!=0) { die "Error running command:    $trimmomatic_command"; }
				print outmet "Quality filtering was done using Trimmomatic (Bolger et al 2014, Bioinformatics 30(15):2114-20)\n";
				}
			}
		}


	my @nmg=keys %allsamples;
	$numsamples=$#nmg+1;
	print outfile3 "Samples:$numsamples\nMode:$mode\n0\n";
	print outfile4 "Samples:$numsamples\nMode:$mode\n";

	if($numsamples==1) { print "$numsamples sample found\n"; }
	else { print "$numsamples samples found: @nmg\n\n"; }


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
  
		if(!$par1files) { print RED; print "There must be at least one 'pair1' sequence file in your samples file $mappingfile, and there is none!\n"; print RESET; die; }
		print outfile4 "Merging files for coassembly: ";
		if($par1files>1) { 
			system("cat $ca1 > $par1name"); 
			my $command="sed -e '/^$/d' $par1name > $par1name.prov; cp $par1name.prov $par1name"; 
			#print "*$command*\n";
			# system($command);
			print outfile4 "cat $ca1 > $par1name; "; 
			} 
		else { system("ln -s $ca1 $par1name"); print outfile4 "ln -s $ca1 $par1name; "; }
		if($par2files>1) { 
			system("cat $ca2 > $par2name"); 
			my $command="sed -e '/^$/d' $par2name > $par2name.prov; cp $par2name.prov $par2name"; 
			# system($command);
			print outfile4 "cat $ca2 > $par2name; "; 
			} 
		elsif($par2files==1) { system("ln -s $ca2 $par2name"); print outfile4 "ln -s $ca2 $par2name; "; }  #-- Support for single reads
		print outfile4 "\n";
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

	if(-e "$tempdir/$projectname.log") { system("rm $tempdir/$projectname.log"); }
	my $rpoint=0;
	my $DAS_Tool_empty=0;


    #-------------------------------- STEP1: Run assembly

		#-- In coassembly mode

	if($mode=~/coassembly/) {
		my $scriptname="01.run_assembly.pl";
		print outfile3 "1\t$scriptname ($assembler)\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n"; 
		print BLUE "[",$currtime->pretty,"]: STEP1 -> RUNNING CO-ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		if($longtrace) { print " (This will take all the metagenomes and assemble them together using $assembler as assembler)\n"; }
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
                	if($longtrace) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
                	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                	if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n"; print RESET; die; }
                	my $wc=qx(wc -l $contigsfna);
                	my($wsize,$rest)=split(/\s+/,$wc);
                	if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname ($assembler). File $contigsfna is empty!\n"; print RESET; die; }
		}

		if($wsize<2)	{ error_out(1,$scriptname,$contigsfna); }
	}

		#-- In merged mode. Includes merging assemblies

	elsif($mode=~/merged|seqmerge/) {
		if(!$extassembly) {
			my $scriptname="01.run_assembly_merged.pl";
			print outfile3 "1\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
			print BLUE "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
			if($longtrace) { print " (This will take all the metagenomes and assemble them independently using $assembler as assembler)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
	                if($ecode!=0)	{ error_out(1,$scriptname); }
			}
	
			#-- Merging individual assemblies 
			#-- We still do it in $extassembly for computing contig lengths and prinseq stuff
 
		my $scriptname;
		if($mode eq "merged") { $scriptname="01.merge_assemblies.pl"; }
		elsif($mode eq "seqmerge") { $scriptname="01.merge_sequential.pl"; }
		print outfile3 "1.5\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP1.5 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP1.5 -> MERGING ASSEMBLIES: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will take all the individual assemblies and merge them in a single, combined assembly)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                if($ecode!=0)   { error_out(1.5,$scriptname); }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($singletons) {
			my $scriptname="01.remap.pl";
               		print outfile3 "1\t$scriptname\n";
                	$currtime=timediff();
                	print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname\n";
                	print BLUE "[",$currtime->pretty,"]: STEP1 ->  ADDING SINGLETONS: $scriptname\n"; print RESET;
                	if($longtrace) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
                	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                	if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname\n"; print RESET; die; }
                	my $wc=qx(wc -l $contigsfna);
                	my($wsize,$rest)=split(/\s+/,$wc);
                	if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname. File $contigsfna is empty!\n"; print RESET; die; }
		}
		if($wsize<2)         { error_out(1.5,$scriptname,$contigsfna); }
	}
	
		#-- In sequential mode. 

	elsif($mode=~/sequential/) {
		my $scriptname="01.run_assembly.pl";
 		print outfile3 "1\t$scriptname\n";
 		$currtime=timediff();
 		print outfile4 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
 		print BLUE "[",$currtime->pretty,"]: STEP1 ->  RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		if($longtrace) { print " (This will take the metagenome and assemble it using $assembler as assembler)\n"; }
 		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(1,$scriptname); }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
                if($wsize<2)    { error_out(1,$scriptname,$contigsfna); }
		if($singletons) {
			my $scriptname="01.remap.pl";
                	print outfile3 "1\t$scriptname\n";
                	$currtime=timediff();
                	print outfile4 "[",$currtime->pretty,"]: STEP1 -> $scriptname\n";
                	print BLUE "[",$currtime->pretty,"]: STEP1 ->  ADDING SINGLETONS: $scriptname ($assembler)\n"; print RESET;
                	if($longtrace) { print " (This will remap reads to contigs and add the unmapped ones as if they were contigs)\n"; }
                	my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                	if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n"; print RESET; die; }
                	my $wc=qx(wc -l $contigsfna);
                	my($wsize,$rest)=split(/\s+/,$wc);		
			if($wsize<2)	{ error_out(1,$scriptname,$contigsfna); }
			}
		}
	close(outfile4); open(outfile4,">>$syslogfile");
			
    #-------------------------------- STEP2: Run RNA prediction

	if(($rpoint<=2) && ((!$test) || ($test>=2))) {
		if($longtrace) { print " At this point, we already have contigs\n"; }
		my $scriptname="02.rnas.pl";
		print outfile3 "2\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP2 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP2 -> RNA PREDICTION: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will run barrnap and Aragorn for predicting putative RNAs in the contigs. This is done before predicting protein-coding genes for avoiding predicting these where there is a RNA)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		my $masked="$interdir/02.$projectname.maskedrna.fasta";
		if($ecode!=0)        { error_out(2,$scriptname); }
		my $wc=qx(wc -l $masked);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(2,$scriptname,$masked); }
		close(outfile4); open(outfile4,">>$syslogfile");	
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if(($rpoint<=3) && ((!$test) || ($test>=3))) { 
		my $scriptname="03.run_prodigal.pl";
 		print outfile3 "3\t$scriptname\n";
 		$currtime=timediff();
 		print outfile4 "\n[",$currtime->pretty,"]: STEP3 -> $scriptname\n";
 		print BLUE "[",$currtime->pretty,"]: STEP3 -> ORF PREDICTION: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will predict putative protein-coding genes in the contigs, using Prodigal)\n"; }
 		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(3,$scriptname); }
		my $wc=qx(wc -l $aafile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(3,$scriptname,$aafile); }
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if(($rpoint<=4) && ((!$test) || ($test>=4))) {
		if($longtrace) { print " At this point, we already have ORFs\n"; }
		my $scriptname="04.rundiamond.pl";
		print outfile3 "4\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP4 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP4 -> HOMOLOGY SEARCHES: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will take all ORFs and run homology searches against functional and taxonomic databases)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(4,$scriptname); }
		my $wc=qx(wc -l $taxdiamond);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<1)    { error_out(4,$scriptname,$taxdiamond); }
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if(($rpoint<=5) && ((!$test) || ($test>=5))) {
		if(!$nopfam) {
			my $scriptname="05.run_hmmer.pl";
			print outfile3 "5\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP5 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP5 -> HMMER/PFAM: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will take all ORFs and run HMMER searches against the PFAM database, for increased sensitivity in predicting protein families. This step is likely to be slow)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(5,$scriptname); }			
			my $wc=qx(wc -l $pfamhmmer);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<4)    { error_out(5,$scriptname,$pfamhmmer); }
		}
	outfile4->autoflush;		
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if(($rpoint<=6) && ((!$test) || ($test>=6))) {
		my $scriptname="06.lca.pl";
		print outfile3 "6\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP6 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP6 -> TAXONOMIC ASSIGNMENT: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will use our last common ancestor (LCA) algorithm to try to annotate the taxonomic origin of each ORF, from the homologues found in the previous step)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(6,$scriptname); }
		my $wc=qx(wc -l "$fun3tax.wranks");
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(6,$scriptname,"$fun3tax.wranks"); }
		close(outfile4); open(outfile4,">>$syslogfile");		
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if(($rpoint<=7) && ((!$test) || ($test>=7))) {
		my $scriptname="07.fun3assign.pl";
		if((!$nocog) || (!$nokegg) || (!$nopfam) || ($opt_db)) {
			print outfile3 "7\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP7 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP7 -> FUNCTIONAL ASSIGNMENT: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will use fun3 algorithm to annotate putative functions for each ORF, from the homologues found in step 5)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)   {  error_out(7,$scriptname); }
			my($wsizeCOG,$wsizeKEGG,$wsizePFAM,$wsizeOPTDB,$rest);
			if(!$nocog) {
				my $wc=qx(wc -l $fun3cog);
				($wsizeCOG,$rest)=split(/\s+/,$wc);
				}
			if(!$nokegg) {
				my $wc=qx(wc -l $fun3kegg);
				($wsizeKEGG,$rest)=split(/\s+/,$wc);
				}
			if(!$nopfam) {
				my $wc=qx(wc -l $fun3pfam);
				($wsizePFAM,$rest)=split(/\s+/,$wc);
				}
			my $optdbsw;
			if($opt_db) {
				open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
				while(<infile0>) {
					my($dbname,$extdb,$dblist)=split(/\t/,$_);
					my $wc=qx(wc -l $resultpath/07.$projectname.fun3.$dbname);
					($wsizeOPTDB,$rest)=split(/\s+/,$wc);
					if($wsizeOPTDB<2) { $optdbsw=$wsizeOPTDB; }
					}
				close infile0;
				}
			if(($wsizeCOG<2) && ($wsizeKEGG<2) && ($wsizePFAM<2) && ($optdbsw<2)) { error_out(7,$scriptname,"$fun3cog, $fun3kegg and $fun3pfam"); }
		}
	close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP8: Blastx on the unannotated parts of the contigs
	
	if(($rpoint<=8) && ((!$test) || ($test>=8))) {
		if($doublepass) {
			my $scriptname="08.blastx.pl";
			# print " DOUBLEPASS: Now starting blastx analysis\n";
			print outfile3 "8\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "\n[",$currtime->pretty,"]: STEP8 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP8 -> DOUBLEPASS, Blastx analysis: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will do many things: it will mask the parts of the contigs where an ORF has been found, and will run blastx in the remaining gaps, to identify possible genes missed in gene prediction. This is intended to be useful when dealing with eukaryotic or viral sequences, for which gene prediction is less accurate)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)  { error_out(8,$scriptname); }
			my $wc=qx(wc -l $gff_file_blastx);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2)         { error_out(8,$scriptname,$gff_file_blastx); }
			}
	close(outfile4); open(outfile4,">>$syslogfile");
	}
		

    #-------------------------------- STEP9: Taxonomic annotation for the contigs (consensus of gene annotations)


	if(($rpoint<=9) && ((!$test) || ($test>=9))) {
		my $scriptname="09.summarycontigs3.pl";
		print outfile3 "9\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP9 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP9 -> CONTIG TAX ASSIGNMENT: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will produce a consensus taxonomic annotation for each contig, according to the annotations of their constituent ORFs)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(9,$scriptname); }
		my $wc=qx(wc -l $alllog);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { error_out(9,$scriptname,$alllog); }
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP10: Mapping of reads onto contigs for abundance calculations
	
	if(($rpoint<=10) && ((!$test) || ($test>=10))) {
		my $scriptname="10.mapsamples.pl";
		print outfile3 "10\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP10 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP10 -> MAPPING READS: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will map reads back to the contigs using $mapper and count how many map to each ORF, to estimate their abundances)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(10,$scriptname); }
		my $wc=qx(wc -l $mapcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { error_out(10,$scriptname,$mapcountfile); }
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP11: Count of taxa abundances
	
	if(($rpoint<=11) && ((!$test) || ($test>=11))) {
		my $scriptname="11.mcount.pl";
		print outfile3 "11\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP11 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP11 -> COUNTING TAX ABUNDANCES: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will count the abundances of each taxon, to infer the composition of the community)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(11,$scriptname); }
		my $wc=qx(wc -l $mcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { error_out(11,$scriptname,$mcountfile); }
		close(outfile4); open(outfile4,">>$syslogfile");		
	}
			
    #-------------------------------- STEP12: Count of function abundances
	
	if(($rpoint<=12) && ((!$test) || ($test>=12))) {
		my $scriptname="12.funcover.pl";
		if((!$nocog) || (!$nokegg) || ($opt_db)) {
		print outfile3 "12\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP12 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP12 -> COUNTING FUNCTION ABUNDANCES: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will count the abundance of each function, to produce a functional profile of the community)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)     { error_out(12,$scriptname); }
		my $cogfuncover="$resultpath/12.$projectname.cog.funcover";
		my $keggfuncover="$resultpath/12.$projectname.kegg.funcover";
		my $wc=qx(wc -l $cogfuncover);
		my($wsizeCOG,$rest)=split(/\s+/,$wc);
		my $wc=qx(wc -l $keggfuncover);
		my($wsizeKEGG,$rest)=split(/\s+/,$wc);
		if(($wsizeCOG<3) && ($wsizeKEGG<3)) { error_out(12,$scriptname,"$cogfuncover and/or $keggfuncover"); }
		}
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP13: Generation of the gene table
		
	if(($rpoint<=13) && ((!$test) || ($test>=13))) {
		my $scriptname="13.mergeannot2.pl";
		print outfile3 "13\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "\n[",$currtime->pretty,"]: STEP13 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP13 -> CREATING GENE TABLE: $scriptname\n"; print RESET;
		if($longtrace) { print " (This will create the gene table by merging all the information compiled in previous steps)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(13,$scriptname); }
		my $wc=qx(wc -l $mergedfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { error_out(13,$scriptname,$mergedfile); }
		if($longtrace) { print " (Now we already have a GENE TABLE)\n"; }
		close(outfile4); open(outfile4,">>$syslogfile");
	}
			
    #-------------------------------- STEP14: Running binning methods 		
	
	 if(!$nobins) {	     
	 	if($longtrace) { print " (Now we will start creating bins for separating individual organisms in the community)\n"; }  
		if(($rpoint<=14) && ((!$test) || ($test>=14))) {
			my $scriptname="14.runbinning.pl";
			print outfile3 "14\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP14 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP14 -> BINNING: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will use binning programs for creating a set of bins)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
		}
		
			
 
    #-------------------------------- STEP15: DAS Tool merging of binning results	
	
		if(($rpoint<=15) && ((!$test) || ($test>=15))) {
			my $scriptname="15.dastool.pl";
			print outfile3 "15\t$scriptname\n";
			$currtime=timediff();
			print outfile4 "[",$currtime->pretty,"]: STEP15 -> $scriptname\n";
			print BLUE "[",$currtime->pretty,"]: STEP15 -> DAS_TOOL MERGING: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will use DASTool for creating a consensus between the sets of bins created in previous steps)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ print RED; print "ERROR in STEP15-> $scriptname\n"; print RESET; }
			my $dirbin=$binresultsdir;
			opendir(indir2,$dirbin) || warn "Can't open $dirbin directory, no DAStool results\n";
			my @binfiles=grep(/fa/,readdir indir2);
			closedir indir2;
			my $firstfile="$dirbin/$binfiles[0]";
			my ($wsize,$rest);
			if(-e $firstfile) {
				my $wc=qx(wc -l $firstfile);
				($wsize,$rest)=split(/\s+/,$wc);
				}
			else { $wsize==0; }
			if($wsize<2) {
				print RED; print "WARNING: File $firstfile is empty!. DAStool did not generate results\n"; print RESET;
				$finaltrace.="DAS Tool abnormal termination: file $firstfile is empty. There are NO BINS!\n";
				$DAS_Tool_empty = 1;
				if($longtrace) { print " (This will use DASTool for creating a consensus between the sets of bins created in previous steps)\n"; 	}
				}
			close(outfile4); open(outfile4,">>$syslogfile");
		}
			
    #-------------------------------- STEP16: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if(($rpoint<=16) && ((!$test) || ($test>=16))) {
			if(!$DAS_Tool_empty){
				my $scriptname="16.addtax2.pl";
				print outfile3 "16\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP16 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP16 -> BIN TAX ASSIGNMENT: $scriptname\n"; print RESET;
				if($longtrace) { print " (This will produce taxonomic assignments for each bin as the consensus of the annotations of each of their contigs)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
				if($ecode!=0){ print RED; print "Stopping in STEP16 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintax);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<1) { print RED; print "Stopping in STEP16 -> $scriptname. File $bintax is empty!\n"; print RESET; die; }
			}
			else{ print RED; print "Skipping BIN TAX ASSIGNMENT: DAS_Tool did not predict bins.\n"; print RESET; }
			close(outfile4); open(outfile4,">>$syslogfile");
		}

			
    #-------------------------------- STEP17: Checking of bins for completeness and contamination (checkM)		
	
		if(($rpoint<=17) && ((!$test) || ($test>=17))) {
			if(!$DAS_Tool_empty){
				my $scriptname="17.checkM_batch.pl";
				print outfile3 "17\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP17 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP17 -> CHECKING BINS: $scriptname\n"; print RESET;
				if($longtrace) { print " (This step will use checkM for estimating completeness and contamination for each bin)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0) { print RED; print "Stopping in STEP17 -> $scriptname\n"; print RESET; die; }
				my $binmethod="DAS";
				my $wc=qx(wc -l $checkmfile);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<4) {
					print RED; print "Can't find $checkmfile\nStopping in STEP18 -> $scriptname\n"; print RESET; die; }
				}
			else { print RED; print"Skipping CHECKM: DAS_Tool did not predict bins.\n"; print RESET; }
		}

			
    #-------------------------------- STEP18: Make bin table		
	
		if(($rpoint<=18) && ((!$test) || ($test>=18))) {
			if(!$DAS_Tool_empty){
				my $scriptname="18.getbins.pl";
				print outfile3 "18\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP18 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP18 -> CREATING BIN TABLE: $scriptname\n"; print RESET;
				if($longtrace) { print " (Now we will compile all previous information for producing a bin table)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP18 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintable);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP18 -> $scriptname. File $bintable is empty!\n"; print RESET; die; }
				if($longtrace) { print " (Now we have a BIN TABLE)\n"; }
			}
			else{ print RED; print "Skipping BIN TABLE CREATION: (You already know: DAS_Tool did not predict bins.)\n"; print RESET; }
		}
		close(outfile4); open(outfile4,">>$syslogfile");
	 }

    #-------------------------------- STEP19: Make contig table		

	if(($rpoint<=19) && ((!$test) || ($test>=19))) {
		my $scriptname="19.getcontigs.pl";
		print outfile3 "19\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP19 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP19 -> CREATING CONTIG TABLE: $scriptname\n"; print RESET;
		if($longtrace) { print " (In this step we will gather all information on contigs for producing a contig table. We do it now because we wanted to know the correspondence between contigs and bins)\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP19 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $contigtable);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { print RED; print "Stopping in STEP19 -> $scriptname. File $contigtable is empty!\n"; print RESET; die; }
		if($longtrace) { print " (Now we have a CONTIG TABLE)\n"; }
		close(outfile4); open(outfile4,">>$syslogfile");
	}

    #-------------------------------- STEP20: Pathways in bins          

	if(!$nobins) {	       
		if(($rpoint<=20) && ((!$test) || ($test>=20))) {
			if((!$DAS_Tool_empty) && (!$nokegg)) {
				my $scriptname="20.minpath.pl";
				print outfile3 "20\t$scriptname\n";
				$currtime=timediff();
				print outfile4 "[",$currtime->pretty,"]: STEP20 -> $scriptname\n";
				print BLUE "[",$currtime->pretty,"]: STEP20 -> CREATING TABLE OF PATHWAYS IN BINS: $scriptname\n"; print RESET;
				if($longtrace) { print " (In this step we will use MinPath and the gene content information for each bin to predict the possible presence of metabolic pathways in it)\n"; }
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP20 -> $scriptname\n"; print RESET; die; }
				my $minpathfile="$resultpath/20.$projectname.kegg.pathways";
				my $wc=qx(wc -l $minpathfile);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP20 -> $scriptname. File $minpathfile is empty!\n"; print RESET; die; }
			}
			else{ print("Skipping MINPATH: DAS_Tool did not predict bins.\n") ; }
		}
		close(outfile4); open(outfile4,">>$syslogfile");
	}


    #-------------------------------- STEP22: Make stats		

	if($rpoint<=21) {
		my $scriptname="21.stats.pl";
		print outfile3 "21\t$scriptname\n";
		$currtime=timediff();
		print outfile4 "[",$currtime->pretty,"]: STEP21 -> $scriptname\n";
		print BLUE "[",$currtime->pretty,"]: STEP21 -> MAKING FINAL STATISTICS: $scriptname\n"; print RESET;
		if($longtrace) { print " (Finally, we will produce final statistics on the run gathering information on genes, contigs and bins\n"; }
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP21 -> $scriptname\n"; print RESET; die; }
		my $statfile="$resultpath/21.$projectname.stats";
		my $wc=qx(wc -l $statfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<10)        { print RED; print "Stopping in STEP21 -> $scriptname. File $statfile is empty!\n"; print RESET; die; }
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
}


sub error_out {
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

