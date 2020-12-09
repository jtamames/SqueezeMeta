#!/usr/bin/perl

#-- Part of the SqueezeMeta distribution. (c) J Tamames Aug 2020

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";
use Term::ANSIColor qw(:constants);
use Getopt::Long;

my $pwd=cwd();
use File::Basename;
use Cwd 'abs_path';

our $utilsdir;
if(-l __FILE__)
	{
	my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
	{
	$utilsdir = abs_path(dirname(__FILE__));
	}
our $installpath = abs_path("$utilsdir/..");
if(-s "$installpath/scripts/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $installpath/scriptsls \n"; }
do "$installpath/scripts/SqueezeMeta_conf.pl";

our($scriptdir,$databasepath,$cdhit_soft,$rdpclassifier_soft, $mothur_soft, $mothur_r, $mothur_t);

my $numthreads=4;
my $seqidtres=0.98;
my $cover_rank=4;		#-- Number of order of the OTU for which calculate coverage
my $covertarget=5;		#-- Required coverage
my $classifier="mothur";		#-- Can be rdp or mothur
my $outputdir="cover";
my ($fileseqs,$hel);
my $datafile="$installpath/data/sizeandcopies.txt";
my $parentsfile="$databasepath/LCA_tax/parents.txt";
my($ecode);


my $helptext = <<END_MESSAGE;
Usage: cover.pl -i <input file> [options]

Arguments:

 Mandatory parameters:
   -i <input file>: Fasta file containing 16S sequences (REQUIRED)
   
 Options  
   -t: Number of threads (Default: $numthreads)
   -idcluster: Identity threshold for collapsing OTUs (Default: $seqidtres)
   -c|-coverage: Target coverage (Default: $covertarget)
   -r|-rank: Rank of target OTU (Default: $cover_rank)
             (Default values imply looking for $covertarget x coverage for the $cover_rank th most abundant OTU)
   -cl|-classifier: Classifier to use (RDP or Mothur) (Default: $classifier)
   -d|-dir: Output directory (Default: $outputdir)	     
   -h|-help: This help
     
END_MESSAGE
   
   

my $result = GetOptions ("t=i" => \$numthreads,
			 "i=s" => \$fileseqs,
			 "idcluster=f" => \$seqidtres,
			 "c|coverage=f" => \$covertarget,
			 "r|rank=i" => \$cover_rank,
			 "cl|classifier=s" => \$classifier,
			 "d|dir=s" => \$outputdir,
			 "h" => \$hel);

if($hel) { die "\n$helptext\n"; } 
if(!$fileseqs) { print BOLD "\nMissing input file\n\n"; print RESET; print "$helptext\n"; die; }
$classifier=~tr/A-Z/a-z/;
if(($classifier ne "rdp") && ($classifier ne "mothur")) { print BOLD "\nUnknown classifier $classifier\n\n"; print RESET; print "$helptext\n"; die; }
if(-e $outputdir) {} else { system("mkdir $outputdir"); }

print BOLD "\nCOVER is part of SqueezeMeta - (c) J. Tamames, F. Puente-SÃ¡nchez CNB-CSIC, Madrid, SPAIN\n\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 9, 3349 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349; Tamames et al, Environ Microbiol Rep. 4:335-41 (2012). doi: 10.1111/j.1758-2229.2012.00338.x\n\n"; print RESET;
print "Looking for coverage $covertarget","x for $cover_rank","th OTU, classifying with $classifier\n\n";


my $outfile="$outputdir/cover.out";
open(out,">$outfile") || die;
print out "#-- Created by $0, ",scalar localtime,", for $fileseqs. \$seqidtres=$seqidtres, $classifier classifier. Looking for $covertarget x coverage for the $cover_rank th most abundant OTU\n";
print out "#-- Citation: Tamames & Puente-Sanchez, Frontiers in Microbiology 9, 3349 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349; Tamames et al, Environ Microbiol Rep. 4(3):335-41 (2012). doi: 10.1111/j.1758-2229.2012.00338.x\n\n"; print RESET;

my @ranks = ( "superkingdom", "phylum", "class", "order", "family", "genus", "species" );
my %cyano2 = (
	"Bacteria;Cyanobacteria;Cyanobacteria;Family III;GpIII",
	"Bacteria;Cyanobacteria;Prochlorales;Prochlorotrichaceae;Prochlorothrix",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family II;GpIIa",
	"Bacteria;Cyanobacteria;Prochlorales;Prochlorococcaceae;Prochlorococcus",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family IV;GpIV",
	"Bacteria;Cyanobacteria;Prochlorales;Prochlorotrichaceae;Prochlorothrix",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family IX;GpIX",
	"Bacteria;Cyanobacteria;Pleurocapsales;Dermocarpella",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family I;GpI",
	"Bacteria;Cyanobacteria;Nostocales;Nostocaceae;Nodularia",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family VII;GpVII",
	"Bacteria;Cyanobacteria;Oscillatoriales;Halospirulina",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family V;GpV",
	"Bacteria;Cyanobacteria;Nostocales;Nostocaceae;Anabaena",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family XIII;GpXIII",
	"Bacteria;Cyanobacteria;Oscillatoriales;Arthrospira",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family XII;GpXII",
	"Bacteria;Cyanobacteria;Oscillatoriales;Oscillatoria",
	"Bacteria;Cyanobacteria;Cyanobacteria;Family XI;GpXI",
	"Bacteria;Cyanobacteria;Chroococcales;Microcystis"
);

my $singletons;

#-- Reading taxonomy

my %parents=('Bacteria','superkingdom:Bacteria','Archaea','superkingdom:Archaea','Eukaryota','superkingdom:Eukaryota');
open(infile3,$parentsfile) || die "Can't open parents file $parentsfile\n";
while(<infile3>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	$tax=~s/ \<.*//g;
	$parents{$tax}=$par;
	}
close infile3;

#-- Reading data for sizes and 16S copies

my(%datos,%abund,%taxa,%corrabund,%copies,%sizes,%toget,%pi);
open(in,$datafile) || die "Cannot open data file $datafile\n";
while(<in>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @d=split(/\t/,$_);
	my($rank,$tax)=split(/\:/,$d[0]);
	$datos{$rank}{$tax}{size}=$d[1]*1e06;
        $datos{$rank}{$tax}{copies}=$d[2];
	}
close in;
    
#-- Running cd-hit for making OTUs

my $filename=$fileseqs;
my @filen=split(/\//,$filename);
$filename=$filen[$#filen];
my $outputcdhit="$outputdir/$filename.cdhit";
if(-e $outputcdhit) { print "  cd-hit results found in $outputcdhit, skipping run\n"; }
else {
	print BOLD "Running cd-hit-est"; print RESET; print " (Schmieder et al 2011, Bioinformatics 27(6):863-4)\n";
	my $command_cdhit = "$cdhit_soft -T $numthreads -c $seqidtres -M 0 -r 1 -l 100 -d 1000 -i $fileseqs -o $outputcdhit > $outputdir/cdhit.log";  
	# print "$command_cdhit\n";
	$ecode = system $command_cdhit;
	if($ecode!=0) { die "Error running command:    $command_cdhit"; }
	print "  Results created in $outputcdhit\n";
	}
	
my($rep,$incl,$totalseqs,$numotus,$singletons);
open(incl,"$outputcdhit.clstr") || die;
while(<incl>) {
	chomp;
	next if !$_;
	if($_=~/^\>/) { 
		$numotus++;
		if($rep) { $abund{$rep}=$incl; }
		($rep,$incl)="";
	}
	elsif($_=~/\>(.*)?\.\.\. \*$/) { $rep=$1; $incl++; $totalseqs++; }
	else { $incl++; $totalseqs++; }
	}
close inc;
if($rep) { $abund{$rep}=$incl; }

foreach my $p(sort keys %abund) {
	if($abund{$p}==1) { $singletons++; }
	}
print "  Number of sequences: $totalseqs\n  Number of OTUs: $numotus\n  Singletons: $singletons\n";
print out "Number of sequences: $totalseqs\nNumber of OTUs: $numotus\nSingletons: $singletons\n";

#-- Unobserved sequences

my $nonobserved = $singletons / $totalseqs ;

printf "  Unobserved fraction: %.2f\n",$nonobserved;
printf out "Unobserved fraction: %.2f\n",$nonobserved;

if($classifier eq "rdp") { rdp(); } 
elsif($classifier eq "mothur") { mothur(); }
else { die "Unknown classifier\n"; }

#-- Adjust abundances by copy numbers

my($totcounts,$uncl_num);
my %totalabund;
foreach my $ttax(sort keys %abund) {
	my $ad=0;
	my $fulltax=$taxa{$ttax};
	my @r=split(/\;/,$fulltax);
	foreach my $ltax(reverse @r) {
		my($lrank,$lname)=split(/\:/,$ltax);
		$totalabund{$lrank}{$lname}+=$abund{$ttax};
		}
	foreach my $ltax(reverse @r) {
		my($lrank,$lname)=split(/\:/,$ltax);
		if($datos{$lrank}{$lname}{copies}) {
			$corrabund{$ttax}=$abund{$ttax}/$datos{$lrank}{$lname}{copies};
			$copies{$ttax}=$datos{$lrank}{$lname}{copies};
			# print "Correcting $ttax by copy number of $ltax ($datos{$lrank}{$lname}{copies}): $abund{$ttax} -> $corrabund{$ttax}\n";
			$totcounts+=$corrabund{$ttax};
			$ad=1;	
			last;
			}
		}
	if(!$ad) { 
		# print "Uncorrected $ttax: $fulltax\n";
		$corrabund{$ttax}=$abund{$ttax};
		$copies{$ttax}=1;
		$totcounts+=$corrabund{$ttax}; 
		}
	}

map{ $corrabund{$_}/=$totcounts; } keys %corrabund;

#-- Output with taxa abundances

open(outt,">$outputdir/abundances.txt") || die;
foreach my $tr(sort keys %totalabund) {
	foreach my $tnam(sort { $totalabund{$tr}{$b}<=>$totalabund{$tr}{$a}; } keys %{ $totalabund{$tr} }) {
		print outt "$tr\t$tnam\t$totalabund{$tr}{$tnam}\n";
		}
	}
close outt;

#-- Calculate sizes and the probability of picking a base from a genome (pi)

my $accumsizes;
foreach my $ttax(sort keys %abund) {
	my $ad=0;
	my $fulltax=$taxa{$ttax};
	my @r=split(/\;/,$fulltax);
	foreach my $ltax(reverse @r) {
		my($lrank,$lname)=split(/\:/,$ltax);
		if($datos{$lrank}{$lname}{size}) { 
			$sizes{$ttax}=$datos{$lrank}{$lname}{size};
			$accumsizes+=($sizes{$ttax}*$corrabund{$ttax});
			$ad=1;
			last;
			}
		}
	if(!$ad) { 
		# print "Unknown size $ttax: $fulltax\n";
		$sizes{$ttax}=4000000;
		$accumsizes+=($sizes{$ttax}*$corrabund{$ttax});
		}
	}
	
foreach my $ttax(sort keys %abund) { $pi{$ttax}=$sizes{$ttax}*$corrabund{$ttax}/$accumsizes; }
	

#-- Calculating requested depth

my($humanread,$humanmag);
my @abundrank=sort { 
	$corrabund{$b}<=>$corrabund{$a} ||  
	$sizes{$b}<=>$sizes{$a};
	} keys %corrabund;
# foreach my $sortrank(@abundrank) { print "$sortrank\t$corrabund{$sortrank}\n"; }
# print "Now estimating sequencing depth for having $covertarget coverage for the $cover_rank","th OTU\n";
my $targetotu=$abundrank[$cover_rank-1];
my $targetabund=$corrabund{$targetotu};
my $targetdepth=int($covertarget*$sizes{$targetotu}/$pi{$targetotu});
# printf "Target: $targetotu, size %d bps, abundance %.3f\n",$sizes{$targetotu},$targetabund;
$humanread=$targetdepth;
if($targetdepth>=1000000000) { $humanread/=1000000000; $humanmag="Gb"; }
elsif($targetdepth>=1000000) { $humanread/=1000000; $humanmag="Mb"; }
else { $humanmag="pb"; }
print BOLD;
printf "\nNeeded %.2f $humanmag, uncorrected\n",$humanread ;
print out "\n\nNeeded $targetdepth bases, uncorrected\n";
my $corrected_depth=int($targetdepth/(1-$nonobserved));
$humanread=$corrected_depth;
if($targetdepth>=1000000000) { $humanread/=1000000000; $humanmag="Gb"; }
elsif($targetdepth>=1000000) { $humanread/=1000000; $humanmag="Mb"; }
else { $humanmag="pb"; }
printf "Correcting by unobserved: %.2f $humanmag\n\n",$humanread;
print out "Correcting by unobserved: $corrected_depth bases\n";
print RESET; 

#-- Calculating coverage for all OTUs

print out "Expected coverages with the proposed sequencing depth:\n";
print out "\nOTU\tSize\tRaw abundance\tRNA copies\tCorrected abundance\tPi\t% Genome sequenced\tCoverage\tTaxon\n";
foreach my $curank(@abundrank) {
	my $cucover=($targetdepth*$pi{$curank})/$sizes{$curank};
	my $genomesequenced = ( 1 - ( exp( -1 * $cucover ) ) ) * 100;
	my $rk=$taxa{$curank};
	my @rk=split(/\;/,$rk);
	my $thistax=$rk[$#rk];
	if(!$thistax) { $thistax="Unclassified"; $uncl_num++; }
	printf out "$curank\t%d\t%.d\t%.1f\t%.3f\t%.4f\t%.2f\t%.1f\t$thistax\n",$sizes{$curank},$abund{$curank},$copies{$curank},$corrabund{$curank},$pi{$curank},$genomesequenced,$cucover;
	}
close out;
print "Finished! Results in $outfile\n";
my $uncratio=$uncl_num/$numotus;
if($uncratio>=0.01) { 
	print RED "WARNING:\n"; print RESET; print "  Some OTUs ($uncl_num out of $numotus) could not be classified\n  For these OTUs, default values of size (4 Mb) and copy number (1) were used\n  This can substantially affect the results, especially if these OTUs are abundant\n";
	if($classifier eq "rdp") { print "  You can try using Mothur as classifier for improving the results\n"; }
	}
	
	
		
sub rdp {

	#-- Running RDP classifier

	my $rdp_outfile="$outputdir/16S.RDP.out";
	my $command="$rdpclassifier_soft classify $outputcdhit -o $rdp_outfile -f filterbyconf";
	print BOLD "Running RDP classifier"; print RESET; print "  (Wang et al 2007, Appl Environ Microbiol 73, 5261-7)\n";
	if(-e $rdp_outfile) { print "  RDP results found in rdp_outfile, skipping run\n"; }
	else {
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		print "  Results created in $rdp_outfile\n";
		}

	#-- Parsing the RDP results

	open(outfile5,">$outputdir/16S.cover.txt") || die "Can't open $outputdir/16S.cover.txt for writing\n";
	print outfile5 "#-- Created by $0, ",scalar localtime,"\n# ORF\tModel\tLast tax\tRank\tFull tax\n";
	open(infile4,$rdp_outfile) || die "Can't open $rdp_outfile\n";
	my @ranks=('superkingdom','phylum','class','order','family','genus','species');
	while(<infile4>) {
		chomp;
		next if !$_;
		my @f=split(/\t/,$_);
		next if(!$f[0]);
		my $lasttax=$f[$#f];
		$lasttax=~s/unclassified\_//;	
		$lasttax=~s/\"//g;
		$lasttax=~s/\_.*//;
		my $phyl=$parents{$lasttax};
		my $rankg;
		for(my $p=1; $p<=$#f; $p++) {
			if($f[$p]!~/unclassified/) { $rankg=$ranks[$p-1]; }
			}
		if((($lasttax=~/^Gp/) || ($lasttax=~/^Subdivision/)) && (!$phyl) && ($_=~/Cyanobacteria/)) { $phyl="superkingdom:Bacteria;phylum:Cyanobacteria"; }
		print outfile5 "$f[0]\t$lasttax\t$rankg\t$phyl\n";
		$taxa{$f[0]}=$phyl;
		}
	close infile4;
	close outfile5;
	}			


sub mothur {

	#-- Running Mothur classifier
	
	my $mothur_outfile="$outputdir/16S.mothur.out";
	my $command="$mothur_soft \"#classify.seqs(fasta=$outputcdhit, taxonomy=$mothur_t, reference=$mothur_r, cutoff=50, processors = $numthreads)\"  > /dev/null 2>&1";
	print BOLD "Running mothur classifier"; print RESET; print " (Schloss et al, Appl Environ Microbiol, 2009. 75(23):7537-41)\n";
	if(-e $mothur_outfile) { print "  Mothur results found in $mothur_outfile, skipping run\n"; }
	else {
		my $ecode = system $command;
		if($ecode!=0) { die "Error running command:    $command"; }
		my $motname=$outputcdhit;
		$motname=~s/\.[^.]+$//;
		$motname.=".nr_v132.wang.taxonomy";
		system("mv $motname $mothur_outfile");
		system("mv mothur.*logfile $outputdir");	
		print "  Results created in $mothur_outfile\n";
		}
	
	#-- Parsing the Mothur results

	open(outfile5,">$outputdir/16S.cover.txt") || die "Can't open $outputdir/16S.cover.txt for writing\n";
	print outfile5 "#-- Created by $0, ",scalar localtime,"\n# ORF\tLast tax\tFull tax\n";
	open(infile4,$mothur_outfile) || die "Can't open $mothur_outfile\n";
	my @ranks=('superkingdom','phylum','class','order','family','genus','species');
	while(<infile4>) {
		chomp;
		next if !$_;
		$_=~s/\;$//;
		my @f=split(/\t/,$_);
		next if(!$f[0]);
		my $oldname=$f[0];
		$f[0]=~s/\_/\:/g;	#-- Mothur things (why do you change sequence names?)
		my @r=split(/\;/,$f[1]);
		my($lasttax,$phyl);
		foreach $lasttax(reverse @r) {
			$lasttax=~s/unclassified\_//;	
			$lasttax=~s/\"//g;
			$lasttax=~s/\_.*//;
			$lasttax=~s/\(.*//;			
			$phyl=$parents{$lasttax};
			last if($phyl);
			}
		my @p=split(/\;/,$phyl);
		my $rankg=$p[$#p];
		if((($lasttax=~/^Gp/) || ($lasttax=~/^Subdivision/)) && (!$phyl) && ($_=~/Cyanobacteria/)) { $phyl="superkingdom:Bacteria;phylum:Cyanobacteria"; }
		print outfile5 "$f[0]\t$rankg\t$phyl\n";
		$taxa{$f[0]}=$phyl; 
		$taxa{$oldname}=$phyl;
		}
	close infile4;
	close outfile5;
	}			
	
