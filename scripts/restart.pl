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

my $version="1.2.0";
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
		print outfile2 "[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING CO-ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir ");
		if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n"; print RESET; die; }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname ($assembler). File $contigsfna is empty!"; print RESET; die;  }
		}

		#-- In merged mode. Includes merging assemblies

	elsif(($mode=~/merged|seqmerge/) && (!$extassembly)) {
		my $scriptname="01.run_assembly_merged.pl";
		print outfile1 "1\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n";  print RESET; die; }
		}
	
		#-- In sequential mode. 

	elsif($mode=~/sequential/) {
 		my $scriptname="01.run_assembly.pl";
 		print outfile1 "1\t$scriptname\n";
 		$currtime=timediff();
 		print outfile2 "[",$currtime->pretty,"]: STEP1 -> $scriptname ($assembler)\n";
 		print CYAN "[",$currtime->pretty,"]: STEP1 -> RUNNING ASSEMBLY: $scriptname ($assembler)\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP1 -> $scriptname ($assembler)\n";  print RESET; die; }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP1 -> $scriptname ($assembler). File $contigsfna is empty!\n";  print RESET; die; }
		}		
	}   	
		
   #-------------------------------- STEP1.5: Merge assemblies

	if(($mode=~/merged|seqmerge/) && ($rpoint<=1.5)) {
		my $scriptname;
		if($mode=~/merged/) { $scriptname="01.merge_assemblies.pl"; }
		elsif($mode=~/seqmerge/) { $scriptname="01.merge_sequential.pl"; }
		print outfile1 "1.5\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP1.5 -> MERGING ASSEMBLIES: $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP1.5 -> $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP1.5 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $contigsfna);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP1.5 -> $scriptname. File $contigsfna is empty!"; print RESET; die; }
	}
	   	
		
    #-------------------------------- STEP2: Run RNA prediction

	if($rpoint<=2) {
		my $scriptname="02.rnas.pl";
		print outfile1 "2\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP2 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP2 -> RNA PREDICTION: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		my $masked="$interdir/02.$projectname.maskedrna.fasta";
		if($ecode!=0)        { print RED; print "Stopping in STEP2 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $masked);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP2 -> $scriptname. File $masked is empty!"; print RESET; die; }
	}
			
    #-------------------------------- STEP3: Run gene prediction

	if($rpoint<=3) {
		my $scriptname="03.run_prodigal.pl";
 		print outfile1 "3\t$scriptname\n";
 		$currtime=timediff();
 		print outfile2 "[",$currtime->pretty,"]: STEP3 -> $scriptname\n";
 		print CYAN "[",$currtime->pretty,"]: STEP3 -> ORF PREDICTION: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP3 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $aafile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP3 -> $scriptname. File $aafile is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP4: Run Diamond for taxa and functions

	if($rpoint<=4) {
		my $scriptname="04.rundiamond.pl";
		print outfile1 "4\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP4 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP4 -> HOMOLOGY SEARCHES: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP4 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $taxdiamond);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<1)         { print RED; print "Stopping in STEP4 -> $scriptname. File $taxdiamond is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP5: Run hmmer for PFAM annotation

	if($rpoint<=5) {
		if(!$nopfam) {
			my $scriptname="05.run_hmmer.pl";
			print outfile1 "5\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP5 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP5 -> HMMER/PFAM: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0){ print RED; print "Stopping in STEP5 -> $scriptname\n"; print RESET; die; }
			my $wc=qx(wc -l $pfamhmmer);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<4) { print RED; print "Stopping in STEP5 -> $scriptname. File $pfamhmmer is empty!\n"; print RESET; die; }
		}
	}
			
    #-------------------------------- STEP6: LCA algorithm for taxa annotation

	if($rpoint<=6) {
		my $scriptname="06.lca.pl";
		print outfile1 "6\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP6 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP6 -> TAXONOMIC ASSIGNMENT: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP6 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l "$fun3tax.wranks");
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP6 -> $scriptname. File $fun3tax is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP7: fun3 for COGs, KEGG and PFAM annotation

	if($rpoint<=7) {
		my $scriptname="07.fun3assign.pl";
		if((!$nocog) || (!$nokegg) || (!$nopfam) || ($opt_db)) {
			print outfile1 "7\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP7 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP7 -> FUNCTIONAL ASSIGNMENT: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)   { print RED; print "Stopping in STEP7 -> $scriptname\n"; print RESET; die; }
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
			if(($wsizeCOG<2) && ($wsizeKEGG<2) && ($wsizePFAM<2) && ($optdbsw<2)) {
				print RED; print "Stopping in STEP7 -> $scriptname. Files $fun3cog, $fun3kegg and $fun3pfam are empty!\n"; print RESET; die; }
		}
	}
			
    #-------------------------------- STEP8: Blastx on the unannotated parts of the contigs
	
	if($rpoint<=8) {
		if($doublepass) {
			my $scriptname="08.blastx.pl";
			print outfile1 "8\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP8 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP8 -> DOUBLEPASS, Blastx analysis: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir");
			if($ecode!=0)  { print RED; print "Stopping in STEP8 -> $scriptname\n"; print RESET; die; }
			my $wc=qx(wc -l $gff_file_blastx);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2)         { print RED; print "Stopping in STEP8 -> $scriptname. File $gff_file_blastx is empty!\n"; print RESET; die; }
			}
	}
		
    #-------------------------------- STEP9: Taxonomic annotation for the contigs (consensus of gene annotations)


	if($rpoint<=9) {
		my $scriptname="09.summarycontigs3.pl";
		print outfile1 "9\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP9 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP9 -> CONTIG TAX ASSIGNMENT: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP9 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $alllog);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP9 -> $scriptname. File $alllog is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP10: Mapping of reads onto contigs for abundance calculations
	
	if($rpoint<=10) {
		my $scriptname="10.mapsamples.pl";
		print outfile1 "10\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP10 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP10 -> MAPPING READS: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP10 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $mapcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { print RED; print "Stopping in STEP10 -> $scriptname. File $mapcountfile is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP11: Count of taxa abundances
	
	if($rpoint<=11) {
		my $scriptname="11.mcount.pl";
		print outfile1 "11\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP11 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP11 -> COUNTING TAX ABUNDANCES: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP11 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $mcountfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<2)         { print RED; print "Stopping in STEP11 -> $scriptname. File $mcountfile is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP12: Count of function abundances
	
	if(($rpoint<=12)) {
		my $scriptname="12.funcover.pl";
		if((!$nocog) || (!$nokegg) || (!$nopfam)) {
		print outfile1 "12\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP12 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP12 -> COUNTING FUNCTION ABUNDANCES: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)     { print RED; print "Stopping in STEP12 -> $scriptname\n"; print RESET; die; }
		my $cogfuncover="$resultpath/12.$projectname.cog.funcover";
		my $keggfuncover="$resultpath/12.$projectname.kegg.funcover";
		my $wc=qx(wc -l $cogfuncover);
		my($wsizeCOG,$rest)=split(/\s+/,$wc);
		my $wc=qx(wc -l $keggfuncover);
		my($wsizeKEGG,$rest)=split(/\s+/,$wc);
		if(($wsizeCOG<3) && ($wsizeKEGG<3)) {
			print RED; print "Stopping in STEP12 -> $scriptname. Files $cogfuncover and/or $keggfuncover are empty!\n"; print RESET; die; }
			}
	}
			
    #-------------------------------- STEP13: Generation of the gene table
		
	if($rpoint<=13) {
		my $scriptname="13.mergeannot2.pl";
		print outfile1 "13\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP13 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP13 -> CREATING GENE TABLE: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP13 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $mergedfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { print RED; print "Stopping in STEP13 -> $scriptname. File $mergedfile is empty!\n"; print RESET; die; }
	}
			
    #-------------------------------- STEP14: Running Maxbin (only for merged or coassembly modes)		
	
	if(!$nobins) {	       
		if(($rpoint<=14) && (!$nomaxbin)) {
			my $scriptname="14.bin_maxbin.pl";
			print outfile1 "14\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP14 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP14 -> MAXBIN BINNING: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ warn "ERROR in STEP14 -> $scriptname\n"; }
			my $dirbin=$bindirs{maxbin};
			opendir(indir1,$dirbin) || die "Can't open $dirbin directory\n";
			my @binfiles=grep(/maxbin.*fasta/,readdir indir1);
                        closedir indir1;
			my $firstfile="$dirbin/$binfiles[0]";
			my $wc=qx(wc -l $firstfile);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2) { warn "WARNING in STEP14 -> $scriptname. File $firstfile is empty, no MaxBin results!\n";  $finaltrace.="WARNING in STEP15: No Maxbin results!\n"; }
		}
			
    #-------------------------------- STEP15: Running Metabat (only for merged or coassembly modes)		
	
		if(($rpoint<=15) && (!$nometabat)) {
			my $scriptname="15.bin_metabat2.pl";
			print outfile1 "15\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP15 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP15 -> METABAT BINNING: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ warn "ERROR in STEP15 -> $scriptname\n"; }
			my $dirbin=$bindirs{metabat2};
			opendir(indir2,$dirbin) || print "Can't open $dirbin directory\n";
			my @binfiles=grep(/fa/,readdir indir2);
			closedir indir2;
			my $firstfile="$dirbin/$binfiles[0]";
			my $wc=qx(wc -l $firstfile);
			my($wsize,$rest)=split(/\s+/,$wc);
			if($wsize<2) { warn "WARNING in STEP15 -> $scriptname. File $firstfile is empty, no Metabat2 results!\n"; $finaltrace.="WARNING in STEP15: No Metabat2 results!\n"; }
		}
 
    #-------------------------------- STEP16: DAS Tool merging of binning results (only for merged or coassembly modes)		
	
		if(($rpoint<=16)) {
			my $scriptname="16.dastool.pl";
			print outfile1 "16\t$scriptname\n";
			$currtime=timediff();
			print outfile2 "[",$currtime->pretty,"]: STEP16 -> $scriptname\n";
			print CYAN "[",$currtime->pretty,"]: STEP16 -> DAS_TOOL MERGING: $scriptname\n"; print RESET;
			my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
			if($ecode!=0){ warn "ERROR in STEP16-> $scriptname\n"; }
			my $dirbin=$dasdir{DASTool};
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
			
    #-------------------------------- STEP17: Taxonomic annotation for the bins (consensus of contig annotations)		
	
		if($rpoint<=17) {
			if(!$DAS_Tool_empty){
				my $scriptname="17.addtax2.pl";
				print outfile1 "17\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP17 -> $scriptname\n";
				print CYAN "[",$currtime->pretty,"]: STEP17 -> BIN TAX ASSIGNMENT: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir >> $tempdir/$projectname.log");
				if($ecode!=0) { print RED; print "Stopping in STEP17 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintax);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<1) { die "Stopping in STEP17 -> $scriptname. File $bintax is empty!\n"; }
			}
		else{ print("Skipping BIN TAX ASSIGNMENT: DAS_Tool did not predict bins.\n"); }
		}

    #-------------------------------- STEP18: Checking of bins for completeness and contamination (checkM)		
	
		if($rpoint<=18) {
			if(!$DAS_Tool_empty){
				my $scriptname="18.checkM_batch.pl";
				print outfile1 "18\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP18 -> $scriptname\n";
				print CYAN "[",$currtime->pretty,"]: STEP18 -> CHECKING BINS: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0) { print RED; print "Stopping in STEP18 -> $scriptname\n"; print RESET; die; }
				foreach my $binmethod(keys %dasdir) {
					$checkmfile="$interdir/18.$projectname.$binmethod.checkM";
					my $wc=qx(wc -l $checkmfile);
					my($wsize,$rest)=split(/\s+/,$wc);
					if($wsize<4) {
						print RED; print "Can't find $checkmfile\nStopping in STEP18 -> $scriptname\n"; print RESET; die; }
				}
			}
			else{ print("Skipping CHECKM: DAS_Tool did not predict bins.\n"); }
		}

			
    #-------------------------------- STEP19: Make bin table		
	
		if($rpoint<=19) {
			if(!$DAS_Tool_empty){
				my $scriptname="19.getbins.pl";
				print outfile1 "19\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP19 -> $scriptname\n";
				print CYAN "[",$currtime->pretty,"]: STEP19 -> CREATING BIN TABLE: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP19 -> $scriptname\n"; print RESET; die; }
				my $wc=qx(wc -l $bintable);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP19 -> $scriptname. File $bintable is empty!\n"; print RESET; die; }
			}
			else{ print("Skipping BIN TABLE CREATION: DAS_Tool did not predict bins.\n") ; }
		}
	}


    #-------------------------------- STEP20: Make contig table		

	if($rpoint<=20) {
		my $scriptname="20.getcontigs.pl";
		print outfile1 "20\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP20 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP20 -> CREATING CONTIG TABLE: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP20 -> $scriptname\n"; print RESET; die; }
		my $wc=qx(wc -l $contigtable);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<3)         { print RED; print "Stopping in STEP20 -> $scriptname. File $contigtable is empty!\n"; print RESET; die; }
	}

    #-------------------------------- STEP21: Pathways in bins          

	if(!$nobins) {	       
  		if($rpoint<=21) {
			if(!$DAS_Tool_empty){
				my $scriptname="21.minpath.pl";
				print outfile1 "21\t$scriptname\n";
				$currtime=timediff();
				print outfile2 "[",$currtime->pretty,"]: STEP21 -> $scriptname\n";
	   	 		print CYAN "[",$currtime->pretty,"]: STEP21 -> CREATING TABLE OF PATHWAYS IN BINS: $scriptname\n"; print RESET;
				my $ecode = system("perl $scriptdir/$scriptname $projectdir");
				if($ecode!=0){ print RED; print "Stopping in STEP21 -> $scriptname\n"; print RESET; die; }
				my $minpathfile="$resultpath/21.$projectname.kegg.pathways";
				my $wc=qx(wc -l $minpathfile);
				my($wsize,$rest)=split(/\s+/,$wc);
				if($wsize<3) { print RED; print "Stopping in STEP21 -> $scriptname. File $minpathfile is empty!\n"; print RESET; die; }
	 		}
		else{ print("Skipping MINPATH: DAS_Tool did not predict bins.\n") ; }
		}
	}

    #-------------------------------- STEP21: Make stats		

	if($rpoint<=22) {
		my $scriptname="22.stats.pl";
		print outfile1 "22\t$scriptname\n";
		$currtime=timediff();
		print outfile2 "[",$currtime->pretty,"]: STEP22 -> $scriptname\n";
		print CYAN "[",$currtime->pretty,"]: STEP22 -> MAKING FINAL STATISTICS: $scriptname\n"; print RESET;
		my $ecode = system("perl $scriptdir/$scriptname $projectdir");
		if($ecode!=0)        { print RED; print "Stopping in STEP22 -> $scriptname\n"; print RESET; die; }
		my $statfile="$resultpath/22.$projectname.stats";
		my $wc=qx(wc -l $statfile);
		my($wsize,$rest)=split(/\s+/,$wc);
		if($wsize<10)        { print RED; print "Stopping in STEP22 -> $scriptname. File $statfile is empty!\n"; print RESET; die; }
	}

    #-------------------------------- END OF PIPELINE		

	print outfile1 "END\n";
	$currtime=timediff();
	print "\nDeleting temporary files in $tempdir\n";
	print outfile2 "\nDeleting temporary files in $tempdir\n";
	system("rm -r $tempdir/*");
	if(-e "$datapath/megahit/final.contigs.fa") { system("rm -r $datapath/megahit/intermediate_contigs; rm $datapath/megahit/final.contigs.fa"); } 
	print outfile2 "[",$currtime->pretty,"]: FINISHED -> Have fun!\n";
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

