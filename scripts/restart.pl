#!/usr/bin/env perl

#-- Restarts interrupted SqueezeMeta processes


$|=1;

use strict;
use Time::Seconds;
use Getopt::Long;
use Cwd;
use Term::ANSIColor qw(:constants);
use lib ".";

#-- Restarts an interrupted pipeline

my $longtrace=0;    #-- Reports an explanation msg for each of the steps
my $start_run = time();

my($rpoint,$hel); 
my $result = GetOptions ("step=i" => \$rpoint,"h" => \$hel);

my $helptext = <<END_MESSAGE;
Usage: restart.pl [options] project

Options:

      -step: Step of the analysis to restart

END_MESSAGE



if($hel) { die "$helptext\n"; } 

my $projectdir=pop @ARGV;
if(!$projectdir) { print RED; print "Please indicate the project to restart\nUsage: restart.pl  [options] project\n"; print RESET;  die; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { print RED; print "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; print RESET;  die; }

do "$projectdir/SqueezeMeta_conf.pl";
do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,$nomaxbin,$nometabat,$lowmem,$minion);
our($nocog,$nokegg,$nopfam,$nobins,$opt_db);
our($numsamples,$numthreads,$mode,$mincontiglen,$assembler,$extassembly,$equivfile,$rawfastq,$evalue,$miniden,$spadesoptions,$megahitoptions,$assembler_options,$doublepass);
our($scriptdir,$databasepath,$extdatapath,$binresultsdir,$interdir,$softdir,$basedir,$datapath,$resultpath,$tempdir,$mappingfile,$contigsfna,$nomaxbin,$contigslen,$mcountfile,$rnafile,$checkmfile,$gff_file,$gff_file_blastx,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$mapcountfile,$contigcov,$checkmfile,$contigtable,$mergedfile,$bintax,$checkmfile,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft,$canu_soft,$trimmomatic_soft,$dastool_soft);
our(%binscripts); 
our($methodsfile, $projectname,$scriptdir,$databasepath,$extdatapath,$interdir,$softdir,$basedir,$datapath,$resultpath,$tempdir,$mappingfile,$contigsfna,$nomaxbin,$contigslen,$mcountfile,$rnafile,$checkmfile,$gff_file,$gff_file_blastx,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$mapcountfile,$contigcov,$contigtable,$mergedfile,$bintax,$checkmfile,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft,$canu_soft,$trimmomatic_soft,$dastool_soft);
our(%bindirs,%dasdir); 
my($finaltrace);
my $progress="$projectdir/progress";

	#-- Read where the process stopped

my $sflag=$rpoint;
my($numsamples,$mode);
open(infile1,$progress) or do { print RED; print "Can't open progress file file $progress\n"; print RESET;  die; }; 
while(<infile1>) {
	chomp $_;
	next if(!$_);
	if($_=~/^Samples\:(\d+)/) { $numsamples=$1; next; }
	if($_=~/^Mode\:(\w+)/) { $mode=$1; next; }
	my $point=$_;
	if(!$sflag) {	($rpoint,my $rest)=split(/\t/,$point); }
	}
close infile1;

	#-- Create new progress, append to existing syslog

my $currtime;
open(outfile1,">$projectdir/progress") or do { print RED; print "Can't open $projectdir/progress for writing\n"; print RESET;  die; }; 
open(outfile2,">>$projectdir/syslog")  or do { print RED; print "Can't open $projectdir/syslog for writing\n"; print RESET;  die; };
$currtime=timediff();
print outfile2 "Restarting project $projectname, ",scalar localtime,"\n";
print outfile1 "Samples:$numsamples\nMode:$mode\n";
if(-e "$tempdir/$projectname.log") { system("rm $tempdir/$projectname.log"); }



#---------------------------------------- PIPELINE --------------------------------
my $DAS_Tool_empty=0;

	if($rpoint<=1) {


	if(-e "$tempdir/$projectname.log") { system("rm $tempdir/$projectname.log"); }
	my $rpoint=0;

    #-------------------------------- STEP1: Run assembly

		#-- In coassembly mode

	if($mode=~/coassembly/) {
		my $scriptname="01.run_assembly.pl";
		print outfile1 "1\t$scriptname ($assembler)\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING CO-ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir ");
		if($ecode!=0)	{ error_out(1,$scriptname); }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)	{ error_out(1,$scriptname,$contigsfna); }
		}

		#-- In merged mode. Includes merging assemblies

	elsif(($mode=~/merged|seqmerge/) && (!$extassembly)) {
		my $scriptname="01.run_assembly_merged.pl";
		print outfile1 "1\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
	        if($ecode!=0)	{ error_out(1,$scriptname); }
		}
	
		#-- In sequential mode. 

	elsif($mode=~/sequential/) {
 		my $scriptname="01.run_assembly.pl";
 		print outfile1 "1\t$scriptname\n";
 		$currtime=timediff();
 		print outfile2 "\n[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
 		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
                if($ecode!=0)   { error_out(1,$scriptname); }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)	{ error_out(1,$scriptname,$contigsfna); }
		}		
	}   	
		
   #-------------------------------- STEP1.5: Merge assemblies

	if(($mode=~/merged|seqmerge/) && ($rpoint<=1.5)) {
		my $scriptname;
		if($mode=~/merged/) { $scriptname="01.merge_assemblies.pl"; }
		elsif($mode=~/seqmerge/) { $scriptname="01.merge_sequential.pl"; }
		print outfile1 "1.5\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP1.5 -> MERGING ASSEMBLIES: $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP1.5 -> $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
             	if($ecode!=0)   { error_out(1.5,$scriptname); }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)	{ error_out(1.5,$scriptname,$contigsfna); }
	}
	   	
		
    #-------------------------------- STEP2: Run RNA prediction

	if($rpoint<=2) {
		my $scriptname="02.rnas.pl";
		print outfile1 "2\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP2 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP2 -> RNA PREDICTION: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		my $masked="$interdir/02.$projectname.maskedrna.fasta";
		if($ecode!=0)        { error_out(2,$scriptname); }
		my $wc=qx(wc -l $masked);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(2,$scriptname,$masked); }
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if($rpoint<=3) {
		my $scriptname="03.run_prodigal.pl";
 		print outfile1 "3\t$scriptname\n";
 		$currtime=timediff();
 		print outfile2 "\n[",$currtime->pretty,"]: STEP3 -> $scriptname\n";
 		print CYAN "[",$currtime->pretty,"]: STEP3 -> ORF PREDICTION: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(3,$scriptname); }
		my $wc=qx(wc -l $aafile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(3,$scriptname,$aafile); }
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if($rpoint<=4) {
		my $scriptname="04.rundiamond.pl";
		print outfile1 "4\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP4 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP4 -> HOMOLOGY SEARCHES: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(4,$scriptname); }
		my $wc=qx(wc -l $taxdiamond);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<1)    { error_out(4,$scriptname,$taxdiamond); }
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if($rpoint<=5) {
		if(!$nopfam) {
			my $scriptname="05.run_hmmer.pl";
			print outfile1 "5\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "\n[",$currtime->pretty,"]: STEP5 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP5 -> HMMER/PFAM: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)        { error_out(5,$scriptname); }			
			my $wc=qx(wc -l $pfamhmmer);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<4)    { error_out(5,$scriptname,$pfamhmmer); }
		}
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if($rpoint<=6) {
		my $scriptname="06.lca.pl";
		print outfile1 "6\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP6 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP6 -> TAXONOMIC ASSIGNMENT: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(6,$scriptname); }
		my $wc=qx(wc -l "$fun3tax.wranks");
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)    { error_out(6,$scriptname,"$fun3tax.wranks"); }
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if($rpoint<=7) {
		my $scriptname="07.fun3assign.pl";
		if((!$nocog) || (!$nokegg) || (!$nopfam) || ($opt_db)) {
			print outfile1 "7\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "\n[",$currtime->pretty,"]: STEP7 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP7 -> FUNCTIONAL ASSIGNMENT: $scriptname\n"; print RESET;
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
	}
			
    #-------------------------------- STEP8: Blastx on the unannotated parts of the contigs
	
	if($rpoint<=8) {
		if($doublepass) {
			my $scriptname="08.blastx.pl";
			print outfile1 "8\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "\n[",$currtime->pretty,"]: STEP8 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP8 -> DOUBLEPASS, Blastx analysis: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)  { error_out(8,$scriptname); }
			my $wc=qx(wc -l $gff_file_blastx);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2)         { error_out(8,$scriptname,$gff_file_blastx); }
			}
	}
		
    #-------------------------------- STEP9: Taxonomic annotation for the contigs (consensus of gene annotations)


	if($rpoint<=9) {
		my $scriptname="09.summarycontigs3.pl";
		print outfile1 "9\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP9 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP9 -> CONTIG TAX ASSIGNMENT: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(9,$scriptname); }
		my $wc=qx(wc -l $alllog);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { error_out(9,$scriptname,$alllog); }
	}
			
    #-------------------------------- STEP10: Mapping of reads onto contigs for abundance calculations
	
	if($rpoint<=10) {
		my $scriptname="10.mapsamples.pl";
		print outfile1 "10\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP10 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP10 -> MAPPING READS: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(10,$scriptname); }
		my $wc=qx(wc -l $mapcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { error_out(10,$scriptname,$mapcountfile); }
	}
			
    #-------------------------------- STEP11: Count of taxa abundances
	
	if($rpoint<=11) {
		my $scriptname="11.mcount.pl";
		print outfile1 "11\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP11 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP11 -> COUNTING TAX ABUNDANCES: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(11,$scriptname); }
		my $wc=qx(wc -l $mcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { error_out(11,$scriptname,$mcountfile); }
	}
			
    #-------------------------------- STEP12: Count of function abundances
	
	if(($rpoint<=12)) {
		my $scriptname="12.funcover.pl";
		if((!$nocog) || (!$nokegg) || ($opt_db)) {
		print outfile1 "12\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP12 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP12 -> COUNTING FUNCTION ABUNDANCES: $scriptname\n"; print RESET;
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
	}
			
    #-------------------------------- STEP13: Generation of the gene table
		
	if($rpoint<=13) {
		my $scriptname="13.mergeannot2.pl";
		print outfile1 "13\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "\n[",$currtime->pretty,"]: STEP13 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP13 -> CREATING GENE TABLE: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { error_out(13,$scriptname); }
		my $wc=qx(wc -l $mergedfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { error_out(13,$scriptname,$mergedfile); }
	}
			
 			
    #-------------------------------- STEP14: Running binning methods 		
	
	 if(!$nobins) {	     
	 	if($longtrace) { print " (Now we will start creating bins for separating individual organisms in the community)\n"; }  
		if($rpoint<=14) {
			my $scriptname="14.runbinning.pl";
			print outfile1 "14\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP14 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP14 -> BINNING: $scriptname\n"; print RESET;
			if($longtrace) { print " (This will use binning programs for creating a set of bins)\n"; }
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0){ print RED; print "ERROR in STEP14 -> $scriptname\n"; print RESET; }
		}
		
			
 
    #-------------------------------- STEP15: DAS Tool merging of binning results (only for merged or coassembly modes)		
	
		if(($rpoint<=15)) {
			my $scriptname="15.dastool.pl";
			print outfile1 "15\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP15 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP15 -> DAS_TOOL MERGING: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ warn "ERROR in STEP15-> $scriptname\n"; }
			my $dirbin=$binresultsdir;
			opendir(indir2,$dirbin)|| die "Can't open $dirbin directory\n";
			my @binfiles=grep(/fa/,readdir indir2);
			closedir indir2;
			my $firstfile="$dirbin/$binfiles[0]";
			my $wc=qx(wc -l $firstfile);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2) {
				print("WARNING: File $firstfile is empty!. DAStool did not generate results\n");
				$DAS_Tool_empty = 1;
				$finaltrace.="DAS Tool abnormal termination: file $firstfile is empty. There are NO BINS!\n";		
			}
		}
			
    #-------------------------------- STEP16: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if($rpoint<=16) {
			if(!$DAS_Tool_empty){
				my $scriptname="16.addtax2.pl";
				print outfile1 "16\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP16 -> $scriptname\n";
				print CYAN "[",$currtime->pretty,"]: STEP16 -> BIN TAX ASSIGNMENT: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
				if($ecode!=0) { print RED; print "Stopping in STEP16 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintax);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<1) { die "Stopping in STEP16 -> $scriptname. File $bintax is empty!\n"; }
			}
		else{ print RED; print("Skipping BIN TAX ASSIGNMENT: DAS_Tool did not predict bins.\n"); print RESET; }
		}

    #-------------------------------- STEP18: Checking of bins for completeness and contamination (checkM)		
	
		if($rpoint<=17) {
			if(!$DAS_Tool_empty){
				my $scriptname="17.checkM_batch.pl";
				print outfile1 "17\t$scriptname\n";
				$currtime=timediff();
				print CYAN "[",$currtime->pretty,"]: STEP17 -> CHECKING BINS: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0) { print RED; print "Stopping in STEP17 -> $scriptname\n"; print RESET; die; }
					my $binmethod="DAS";
					my $wc=qx(wc -l $checkmfile);
					my($wsize,$rest)=split(/\s+/,$wc);
					if($wsize<4) {
						print RED; print "Can't find $checkmfile\nStopping in STEP18 -> $scriptname\n"; print RESET; die; }
			}
			else{ print RED; print("Skipping CHECKM: DAS_Tool did not predict bins.\n"); print RESET; }
		}

			
    #-------------------------------- STEP18: Make bin table		
	
		if($rpoint<=18) {
			if(!$DAS_Tool_empty){
				my $scriptname="18.getbins.pl";
				print outfile1 "18\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP18 -> $scriptname\n";
				print CYAN "[",$currtime->pretty,"]: STEP18 -> CREATING BIN TABLE: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP18 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintable);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP18 -> $scriptname. File $bintable is empty!\n"; print RESET; die; }
			}
			else{ print RED; print("Skipping BIN TABLE CREATION: (You already know: DAS_Tool did not predict bins.)\n") ; print RESET; }
		}
	}


    #-------------------------------- STEP19: Make contig table		

	if($rpoint<=19) {
		my $scriptname="19.getcontigs.pl";
		print outfile1 "19\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP19 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP19 -> CREATING CONTIG TABLE: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP19 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $contigtable);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { print RED; print "Stopping in STEP19 -> $scriptname. File $contigtable is empty!\n"; print RESET; die; }
	}

    #-------------------------------- STEP20: Pathways in bins          

	if(!$nobins) {	       
  		if($rpoint<=20) {
			if(!$DAS_Tool_empty){
				my $scriptname="20.minpath.pl";
				print outfile1 "20\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP20 -> $scriptname\n";
	   	 		print CYAN "[",$currtime->pretty,"]: STEP20 -> CREATING TABLE OF PATHWAYS IN BINS: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP20 -> $scriptname\n"; print RESET; die; }
				my $minpathfile="$resultpath/20.$projectname.kegg.pathways";
				my $wc=qx(wc -l $minpathfile);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP20 -> $scriptname. File $minpathfile is empty!\n"; print RESET; die; }
	 		}
		else{ print("Skipping MINPATH: DAS_Tool did not predict bins.\n") ; }
		}
	}

    #-------------------------------- STEP22: Make stats		

	if($rpoint<=21) {
		my $scriptname="21.stats.pl";
		print outfile1 "21\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP21 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP21 -> MAKING FINAL STATISTICS: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP21 -> $scriptname\n"; print RESET; die; }
		my $statfile="$resultpath/21.$projectname.stats";
		my $wc=qx(wc -l $statfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<10)        { print RED; print "Stopping in STEP21 -> $scriptname. File $statfile is empty!\n"; print RESET; die; }
	}

    #-------------------------------- END OF PIPELINE		

	print outfile1 "END\n";
	$currtime=timediff();
	print "\nDeleting temporary files in $tempdir\n";
	print outfile2 "\nDeleting temporary files in $tempdir\n";
	system("rm -r $tempdir/*");
	if(-e "$datapath/megahit/final.contigs.fa") { system("rm -r $datapath/megahit/intermediate_contigs; rm $datapath/megahit/final.contigs.fa"); } 
	print outfile2 "\n[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
	print CYAN "[",$currtime->pretty,"]: FINISHED -> Have fun!\n"; print RESET;
	if($finaltrace) { print "\nWARNINGS:\n$finaltrace\n"; }
	print "For citation purposes, you can find a summary of methods in the file $methodsfile\n\n";


#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}

#---------------------------------------- ERROR LOG --------------------------------


sub error_out {
	close outfile2;
	open(outfile2,">>$projectdir/syslog");		#-- This is just to flush syslog before printing errors
	my($step,$scriptname,$file)=@_;
	if($file) { print RED; print "Stopping in STEP$step -> $scriptname. File $file is empty!\n"; print RESET; print outfile2 "Stopping in STEP$step -> $scriptname. File $file is empty!\n"; }
	else { print RED; print "Stopping in STEP$step -> $scriptname. Program finished abnormally\n"; print RESET; print outfile2 "Stopping in STEP$step -> $scriptname. Program finished abnormally\n"; }
	print outfile2 "_____________\n\nSystem information:\n";
	system("uname -a >> $projectdir/syslog");
	print outfile2 "_____________\n\nTree for the project:\n";
	system("tree -h -D $projectdir >> $projectdir/syslog");
	close outfile2;
	print RED; print "\n  If you don't know what went wrong or want further advice, please look for similar issues in https://github.com/jtamames/SqueezeMeta/issues\n";
	print "  Feel free to open a new issue if you don't find the answer there. Please add a brief description of the problem and upload the $projectdir/syslog file (zip it first)\n";
	print RESET;
	die;
}

