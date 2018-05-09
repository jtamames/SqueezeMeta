#!/usr/bin/perl

# v0.1.0 07/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
# v0.1.1 07(05/2018 Added dynamic path detection. FPS.

#-- Restarts interrupted squeezeM processes


$|=1;

use strict;
use Time::Seconds;
use Cwd;

#-- Restarts an interrupted pipeline

my $start_run = time();

###scriptdir patch, Fernando Puente-Sánchez, 07-V-2018
use File::Basename;
our $scriptdir = dirname(__FILE__);
our $installpath = "$scriptdir/..";
###

my $pwd=cwd();	

my $project=$ARGV[0];				#-- THIS MUST POINT TO THE INTERRUPTED PROCESS		
if(!$project) { die "Please indicate the project to restart\nUsage: restart.pl <project>\n"; }

my $progress="$pwd/$project/progress";

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($datapath,$assembler,$outassembly,);
our($nocog,$nokegg,$nopfam,$nobins);
our($numsamples,$numthreads,$mode,$mincontiglen,$assembler,$project,$equivfile,$rawfastq,$evalue,$miniden,$spadesoptions,$megahitoptions,$assembler_options);
our($databasepath,$extdatapath,$softdir,$basedir,$datapath,$resultpath,$tempdir,$mappingfile,$contigsfna,$contigslen,$mcountfile,$rnafile,$gff_file,$aafile,$ntfile,$daafile,$taxdiamond,$cogdiamond,$keggdiamond,$pfamhmmer,$fun3tax,$fun3kegg,$fun3cog,$fun3pfam,$allorfs,$alllog,$rpkmfile,$coveragefile,$contigcov,$contigtable,$mergedfile,$bintax,$checkmfile,$bincov,$bintable,$contigsinbins,$coglist,$kegglist,$pfamlist,$taxlist,$nr_db,$cog_db,$kegg_db,$lca_db,$bowtieref,$pfam_db,$metabat_soft,$maxbin_soft,$spades_soft,$barrnap_soft,$bowtie2_build_soft,$bowtie2_x_soft,$bedtools_soft,$diamond_soft,$hmmer_soft,$megahit_soft,$prinseq_soft,$prodigal_soft,$cdhit_soft,$toamos_soft,$minimus2_soft);
our %bindirs;  

	#-- Read where the process stopped

my($numsamples,$mode,$rpoint);
open(infile1,$progress) || die; 
while(<infile1>) {
	chomp $_;
	next if(!$_);
	if($_=~/^Samples\:(\d+)/) { $numsamples=$1; next; }
	if($_=~/^Mode\:(\d+)/) { $mode=$1; next; }
	my $point=$_;
	($rpoint,my $rest)=split(/\t/,$point);
	}
close infile1;

	#-- Create new progress, append to existing syslog

my $currtime;
open(outfile1,">$pwd/$project/progress") || die;  
open(outfile2,">>$pwd/$project/syslog") || die;
$currtime=timediff();
print outfile2 "Restarting project $project, ",scalar localtime,"\n";
print outfile1 "Samples:$numsamples\nMode:$mode\n";
system("rm $tempdir/$project.log");



#---------------------------------------- PIPELINE --------------------------------

	if($rpoint<=1) {


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



#---------------------------------------- TIME CALCULATIONS --------------------------------

sub timediff {
	my $end_run = time();
	my $run_time = $end_run - $start_run;
	my $timesp = Time::Seconds->new( $run_time );
	return $timesp;
}

