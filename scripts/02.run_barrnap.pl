#!/usr/bin/perl

#-- Part of squeezeM distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Runs barrnap for RNA prediction. Sequentially runs model for bacteria, archaea, eukaryote and mitochondria

use strict;
use warnings;
use Cwd;
use Tie::IxHash;

my $pwd=cwd();

my $project=$ARGV[0];

do "$project/squeezeM_conf.pl";

#-- Configuration variables from conf file

our($resultpath,$contigsfna,$tempdir,$barrnap_soft,$numthreads,$rnafile);

my %king;
tie %king,"Tie::IxHash";

%king=('bac','Bacteria','arc','Archaea','euk','Eukaryote','mito','Mitochondrial');
my $targetfile="$resultpath/02.$project.maskedrna.fasta";
system("cp $contigsfna $targetfile");
if(-e $rnafile) { system("rm $rnafile"); }

#-- Loop for all kingdoms (Bac, Arch, Euk) plus mitochondria, looking for the respective RNAs

foreach my $kingdom(keys %king) {
	my(%rna,%inrna)=();
	my $output="$tempdir/$kingdom.gff";

	#-- Run barrnap

	my $command="$barrnap_soft --quiet --threads $numthreads --kingdom $kingdom --reject 0.1 $targetfile > $output";
	print "Running barrnap for $king{$kingdom}: $command\n";
	system $command;

	#-- Reformat the output, adding the type of RNA found and the ORF ID

	open(outfile1,">$output.mod") || die;
	open(infile1,$output) || die;
	my $mol;
	while(<infile1>) {
 		chomp;
		if($_=~/^\#/) { print outfile1 "$_\n"; next; }
		my @k=split(/\t/,$_);
		my $idx="$k[3]-$k[4]"; 
		$inrna{$k[0]}++;
		my $newid="$k[0]\_RNA$inrna{$k[0]}";
		if($k[8]=~/Name\=(.*)?\_/) { $mol=$1; } else  { $mol=""; }
		$rna{$k[0]}{$idx}{molecule}="$mol rRNA";
		$rna{$k[0]}{$idx}{name}=$newid ; 
		$rna{$k[0]}{$idx}{kingdom}=$king{$kingdom} ;
		$k[8]="ID=$newid;".$k[8];                  # We must add the ID in the gff, for reading abundances later
		my $newline=join("\t",@k);
		print outfile1 "$newline\n";
		}
	close infile1;
	close outfile1;

	#-- Concatenate all RNA files, and mask the contigs for not predicting these RNAs as proteins (in upcoming gene prediction)

	open(outfile2,">>$rnafile");
	open(outfile3,">contigs.prov");
	open(infile2,$targetfile) || die;
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
	system("mv contigs.prov $targetfile");
}
		
#-- Creating the RNAs gff file

my $gffout="$tempdir/02.$project.rna.gff";			     		
if(-e $gffout) { system("rm $gffout"); }
my $command="cat $tempdir/*gff.mod > $gffout";
system($command);
