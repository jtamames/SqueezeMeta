#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs barrnap for RNA prediction. Sequentially runs model for bacteria, archaea, eukaryote and mitochondria
#-- Runs aragorn for tRNA/tmRNA prediction (21/11/2019)

use strict;
use Cwd;
use Tie::IxHash;
use lib ".";

my $pwd=cwd();

my $projectdir=$ARGV[0];
if(!$projectdir) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectdir/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectdir. Is the project path ok?"; }
do "$projectdir/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectdir/parameters.pl";

#-- Configuration variables from conf file

our($resultpath,$interdir,$contigsfna,$tempdir,$barrnap_soft,$rdpclassifier_soft,$numthreads,$rnafile,$databasepath,$aragorn_soft,$trnafile,$methodsfile,$syslogfile);

open(outmet,">>$methodsfile") || warn "Cannot open methods file $methodsfile for writing methods and references\n";
open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

my %king;
my $command;
tie %king,"Tie::IxHash";

%king=('bac','Bacteria','arc','Archaea','euk','Eukaryote','mito','Mitochondrial');
my $targetfile="$interdir/02.$project.maskedrna.fasta";
system("cp $contigsfna $targetfile");
if(-e $rnafile) { system("rm $rnafile"); }

#-- Loop for all kingdoms (Bac, Arch, Euk) plus mitochondria, looking for the respective RNAs

print "  Running barrnap (Seeman 2014, Bioinformatics 30, 2068-9) for predicting RNAs: ";
my %rname;
open(outfile4,">$tempdir/16S.fasta") || die "Can't open $tempdir/16S.fasta for writing\n";
foreach my $kingdom(keys %king) {
	my(%rna,%inrna)=();
	my $output="$tempdir/$kingdom.gff";

	#-- Run barrnap

	my $command="$barrnap_soft --quiet --threads $numthreads --kingdom $kingdom --reject 0.1 $targetfile --dbdir $databasepath > $output";
	print outsyslog "Running barrnap for $king{$kingdom}: $command\n";
	print " $king{$kingdom}";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:    $command"; }

	#-- Reformat the output, adding the type of RNA found and the ORF ID

	open(outfile1,">$output.mod") || die "Can't open $output.mod for writing\n";
	open(infile1,$output) || die "Can't open $output\n";
	my $mol;
	while(<infile1>) {
 		chomp;
		if($_=~/^\#/) { print outfile1 "$_\n"; next; }
		my @k=split(/\t/,$_);
		my $idx="$k[3]-$k[4]"; 
		$inrna{$k[0]}++;
		my $newid="$k[0]\_$idx";
		if($k[8]=~/Name\=(.*)?\_/) { $mol=$1; } else  { $mol=""; }
		$rna{$k[0]}{$idx}{molecule}="$mol rRNA";
		$rna{$k[0]}{$idx}{name}=$newid ; 
		$rna{$k[0]}{$idx}{kingdom}=$king{$kingdom} ;
		$rname{$newid}=$king{$kingdom} ;
		$k[8]="ID=$newid;".$k[8];                  # We must add the ID in the gff, for reading abundances later
		my $newline=join("\t",@k);
		print outfile1 "$newline\n";
		}
	close infile1;
	close outfile1;

	#-- Concatenate all RNA files, and mask the contigs for not predicting these RNAs as proteins (in upcoming gene prediction)

	open(outfile2,">>$rnafile") || die "Can't open $rnafile for writing\n";
	open(outfile3,">$tempdir/contigs.prov") || die "Can't open contigs.prov for writing\n";
	open(infile2,$targetfile) || die "Can't open $targetfile\n";
	my($seq,$current)="";
	while(<infile2>) {
		chomp;
		if($_=~/^\>([^ ]+)/) { 
			my $contigname=$1; 
			if($current) {

				#-- Masking with 'N's

				foreach my $rns(sort keys %{ $rna{$current} }) {
				my($init,$end)=split(/\-/,$rns);
				my $longr=($end-$init)+1;
				my $replace=('N' x $longr);
				my $rnaseq=substr($seq,$init-1,$longr,$replace);
				print outfile2 ">$rna{$current}{$rns}{name}\t$rna{$current}{$rns}{molecule};$current|$rns;$rna{$current}{$rns}{kingdom}\n$rnaseq\n";
				if($rna{$current}{$rns}{molecule}=~/16S/) { print outfile4 ">$rna{$current}{$rns}{name}\t$rna{$current}{$rns}{molecule};$current|$rns;$rna{$current}{$rns}{kingdom}\n$rnaseq\n"; } 
				}
			}
			if($current) { print outfile3 ">$current\n$seq\n"; }	
			$seq="";
			$current=$contigname;     
		}
		else { $seq.=$_; }
	}

	print outfile3 ">$current\n$seq\n";	
	close infile2;
	close outfile2;
	close outfile3;
	system("mv $tempdir/contigs.prov $targetfile");
}
print "\n";
close outfile4;
print outmet "RNAs were predicted using Barrnap (Seeman 2014, Bioinformatics 30, 2068-9)\n";
		

#-- Running RDP classifier for 16S sequences

$command="$rdpclassifier_soft classify $tempdir/16S.fasta -o $tempdir/16S.out -f filterbyconf";
print outsyslog "Running RDP classifier: $command\n";
print "  Running RDP classifier (Wang et al 2007, Appl Environ Microbiol 73, 5261-7)\n";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:    $command"; }
print outmet "16S rRNA sequences were taxonomically classified using the RDP classifier (Wang et al 2007, Appl Environ Microbiol 73, 5261-7)\n";

my %parents=('Bacteria','superkingdom:Bacteria','Archea','superkingdom:Archaea','Eukaryota','superkingdom:Eukaryota');
open(infile3,"$databasepath/LCA_tax/parents.txt") || die "Can't open $databasepath/LCA_tax/parents.txt\n";
while(<infile3>) {
	chomp;
	next if !$_;
	my ($tax,$par)=split(/\t/,$_);
	$tax=~s/ \<.*//g;
	$parents{$tax}=$par;
	}
close infile3;

open(outfile5,">$resultpath/02.$project.16S.txt") || die "Can't open $resultpath/02.$project.16S.txt for writing\n";
print outfile5 "#-- Created by $0, ",scalar localtime,"\n# ORF\tModel\tLast tax\tRank\tFull tax\n";
open(infile4,"$tempdir/16S.out") || die "Can't open $tempdir/16S.out\n";
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
	print outfile5 "$f[0]\t$rname{$f[0]}\t$lasttax\t$rankg\t$phyl\n";
	}
close infile4;
close outfile5;

#-- Running Aragorn

print "  Running Aragorn (Laslett & Canback 2004, Nucleic Acids Res 31, 11-16) for tRNA/tmRNA prediction\n";
my $temparagorn="$tempdir/trnas.aragorn";
$command="$aragorn_soft -w $targetfile -o $temparagorn";
print outsyslog "Running Aragorn: $command\n";
system($command);
open(infile5,$temparagorn) || die "Cannot open Aragorn result file $temparagorn\n";
open(outfile6,">$trnafile") || die;
open(outfile7,">$tempdir/trna.gff.mod") || die;
my($incontig);
my %trnas;
while(<infile5>) {
	chomp;
	next if !$_;
	if($_=~/^\>([^ ]+)/) { $incontig=$1; next; }
	if($_=~/tRNA|tmRNA/) {
		my @fields=split(/\s+/,$_);
		my($pos,$rest)=split(/\t/,$fields[2]);
		$pos=~s/\[//;
		$pos=~s/\]//;
		$pos=~s/c//;
		my($pos1,$pos2)=split(/\,/,$pos);
		if($pos1<0) { $pos1=1; }
		if($pos2<0) { $pos2=1; }		
		my $dire="+";
		if($pos1>$pos2) { $pos="$pos2-$pos1"; my $tpo=$pos1; $pos1=$pos2; $pos2=$tpo; $dire="-"; }
		else { $pos="$pos1-$pos2"; }
		my $genname="$incontig\_$pos";
		print outfile6 "$genname\t$fields[1]\n";
		$trnas{$incontig}{$pos}=1;
		print outfile7 "$incontig\taragorn\ttRNA\t$pos1\t$pos2\t-\t$dire\t.\tID=$genname;Name=$fields[1]\n";
		}
	}
close infile5;
print outmet "tRNA/tmRNA sequences were predicted using Aragorn (Laslett & Canback 2004, Nucleic Acids Res 31, 11-16)\n";
close outfile6;
close outfile7;
		
		
	#-- Masking with 'N's
	
open(outfile7,">$tempdir/contigs.prov") || die "Can't open contigs.prov for writing\n";

open(outfile8,">$resultpath/02.$project.trnas.fasta") || die "Can't open trna file for writing\n";
open(infile6,$targetfile) || die "Can't open $targetfile\n";
my($seq,$current)="";
while(<infile6>) {
	chomp;
	if($_=~/^\>([^ ]+)/) { 
		my $contigname=$1; 
		if($current) {

			#-- Masking with 'N's

			foreach my $rns(sort keys %{ $trnas{$current} }) {
			my($init,$end)=split(/\-/,$rns);
			my $longr=($end-$init)+1;
			my $replace=('N' x $longr);
			my $rnaseq=substr($seq,$init-1,$longr,$replace);
			print outfile8 ">$current\_$rns\n$rnaseq\n";
			}
		}
		if($current) { print outfile7 ">$current\n$seq\n"; }	
		$seq="";
		$current=$contigname;     
	}
	else { $seq.=$_; }
}
if($current) {	#-- Last entry
	foreach my $rns(sort keys %{ $trnas{$current} }) {
		my($init,$end)=split(/\-/,$rns);
		my $longr=($end-$init)+1;
		my $replace=('N' x $longr);
		my $rnaseq=substr($seq,$init-1,$longr,$replace);
		print outfile8 ">$current\_$rns\n$rnaseq\n";
		}
	print outfile7 ">$current\n$seq\n"; 	
	}
close infile6;
close outfile7;
close outfile8;
close outmet;

system("mv $tempdir/contigs.prov $targetfile");
		
#-- Creating the RNAs gff file

my $gffout="$tempdir/02.$project.rna.gff";			     		
if(-e $gffout) { system("rm $gffout"); }
$command="cat $tempdir/*gff.mod > $gffout";
print outsyslog "Creating new gff file: $command\n";
system($command);
close outsyslog;
