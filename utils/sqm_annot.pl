#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

$|=1;

my $commandline=$0 . " ". (join " ", @ARGV);

use Time::Seconds;
use Cwd;
use Getopt::Long;
use Tie::IxHash;
use Linux::MemInfo;
use Term::ANSIColor qw(:constants);
use lib ".";
use strict;
use threads;

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;
use Cwd 'abs_path';

my($symlinkpath,$symlinkdest,$utilsdir);
if(-l __FILE__)
        {
        $symlinkpath = dirname(__FILE__);
        $symlinkdest = readlink(__FILE__);
        $utilsdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $utilsdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$utilsdir/..");
my $scriptdir = "$installpath/scripts";
my $auxdir = "$installpath/lib/SQM_reads";
my $datadir = "$installpath/data/";
###

open(inv,"$installpath/VERSION") || die;
my $version=<inv>;
chomp $version;
close inv;

my($numthreads,$aadir,$project,$samplesfile,$notax,$nocog,$nokegg,$blastmode,$blocksize,$hel,$blockoption,$printversion,$opt_db);

my $start_run = time();
print BOLD "\nSQM_annot v$version - (c) J. Tamames, F. Puente-Sánchez CNB-CSIC, Madrid, SPAIN\n\nThis is part of the SqueezeMeta distribution (https://github.com/jtamames/SqueezeMeta)\nPlease cite: Tamames & Puente-Sanchez, Frontiers in Microbiology 10.3389 (2019). doi: https://doi.org/10.3389/fmicb.2018.03349\n\n"; print RESET;

my $helpshort="Usage: SQM_annot.pl -f <faa directory> -s <samples file> [options]\n";

my $helptext = <<END_MESSAGE;
Usage: SQM_annot.pl -f <faa directory> -s <samples file> [options]

Mandatory parameters:
   -f: Directory for sequence files (REQUIRED)
   -s: Samples file (REQUIRED)

 Options:
   -t: Number of threads (Default: 12)
   -b: Diamond block size
   --notax: Skip taxonomic annotation
   --nocog: Skip COGs annotation
   --nofun: Skip KEGG annotation
   -extdb <database file>: List of user-provided databases
   -version: Print version
   -h: this help
   
END_MESSAGE

my $result = GetOptions ("t=i" => \($numthreads=12),
                     "s=s" => \$samplesfile,
                     "f=s" => \$aadir,
		     "b=i" => \$blocksize,
		     "notax" => \$notax,
		     "nocog" => \$nocog,
		     "extdb=s" => \$opt_db, 
		     "nokegg" => \$nokegg,
		     "v|version" => \$printversion,
		     "h" => \$hel
		    );

if($hel) { print "$helptext\n"; exit; }
if($printversion) { exit; }
if(!$samplesfile) { die "Please specify a file of samples\n"; }
if(!$aadir) { $aadir="."; }
if($blocksize) { $blockoption="-b $blocksize"; }
if($blastmode=~/blastx/) { $blastmode="blastx"; } else { $blastmode="blastp"; }
if(!$notax) { $notax="0"; }
my $coglist="$datadir/coglist.txt";
my $kegglist="$datadir/keggfun2.txt";

my(%files,%type);
tie %files,"Tie::IxHash";
open(in,$samplesfile) || die; 
my $tempsample="temp.samples";  #-- We need a samples file to trick SqueezeMeta
my $anum;
while(-e $tempsample) { 	#-- For preventing overwriting an already existing "temp.samples" file (happened already)
	$anum++;
	$tempsample="temp.samples.$anum";
	}
open(out,">$tempsample") || die;
while(<in>) {
	chomp;
	my @e=split(/\t/,$_);
	$files{$e[0]}=$e[1];
	$type{$e[0]}=$e[2];
	my $fname=$e[1];
	$fname=~s/\.gz$//;
	print out "$e[0]\t$fname\tpair1\n";
	}
close out;

my $edb;
if($opt_db) { $edb="-extdb $opt_db"; }
print "\nNow I will call SqueezeMeta to do my stuff. Please hold on.\n\n";
my $command="perl $scriptdir/SqueezeMeta.pl  -s $tempsample -f $aadir -m sequential $edb $blockoption -t $numthreads --nopfam -c 0 --empty";
my $ecode = system($command);
# if($ecode!=0) { catch_error(); }
my $diamond_command;
foreach my $thissample(keys %files) {
	print "Working with $thissample\n";
	my $aafile=$files{$thissample};
	if($aafile=~/\.gz$/) { system("gunzip $aafile"); $aafile=~s/\.gz$//; }
	$project=$thissample;
	my $typefile=$type{$thissample};
	if($typefile=~/genome/i) {
		my $command="cp $aadir/$aafile SQM/$project/results/01.$project.fasta";
		system($command);
		print "***$command***n";
		$command="perl $scriptdir/02.rnas.pl SQM/$project; perl $scriptdir/03.run_prodigal.pl SQM/$project;";
		my $ecode = system $command;
		if($ecode!=0) { catch_error(); }
		$blastmode="blastp";
		}
	elsif($typefile=~/aa/i) { 
		system("cp $aadir/$aafile SQM/$project/results/03.$project.faa"); 
		$blastmode="blastp";
	}
	elsif($typefile=~/nt/i) { 
		system("cp $aadir/$aafile SQM/$project/results/03.$project.faa"); 
		$blastmode="blastx";
	}
	$diamond_command="perl $scriptdir/04.rundiamond.pl SQM/$project $notax $blastmode";
	$command="perl $scriptdir/04.rundiamond.pl SQM/$project $notax $blastmode";
	my $ecode = system $command;
	if($ecode!=0) { catch_error(); }
	if($typefile=~/nt/i) { system("mv SQM/$project/results/03.$project.faa SQM/$project/results/03.$project.fna"); }

	if(!$notax) { 
		$command="perl $scriptdir/06.lca.pl SQM/$project"; 
		my $ecode = system $command; 
		if($ecode!=0) { catch_error(); }
		}
	if((!$nocog) || (!$nokegg)) { 
		$command="perl $scriptdir/07.fun3assign.pl SQM/$project"; 
		my $ecode = system $command; 
		if($ecode!=0) { catch_error(); }
		}
print "\n";
if(!$notax) { print "Taxonomic assignment stored in SQM/$project/results/06.$project.fun3.tax.wranks\n"; }
if(!$nocog) { print "COG functional assignment stored in SQM/$project/results/07.$project.fun3.cog\n"; }	
if(!$nokegg) { print "KEGG functional assignment stored in SQM/$project/results/07.$project.fun3.kegg\n"; }	

system("rm $tempsample");
funclass();

	}
	
print "Have a nice day!\n";


sub funclass {

	my(%funs,%funannot,%accumc);

	#-- Reading KEGG functions and pathways

	open(infile1,$kegglist) || warn "Missing KEGG equivalence file $kegglist\n";
	while(<infile1>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my @t=split(/\t/,$_);
		$funs{KEGG}{$t[0]}{name}=$t[1];
		$funs{KEGG}{$t[0]}{fun}=$t[2];
		$funs{KEGG}{$t[0]}{path}=$t[3];
		}
	close infile1;

	#-- Reading COG functions and pathways

	open(infile2,$coglist) || warn "Missing COG equivalence file $coglist\n";
	while(<infile2>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my @t=split(/\t/,$_);
		$funs{COG}{$t[0]}{name}=$t[0];
		$funs{COG}{$t[0]}{fun}=$t[1];
		$funs{COG}{$t[0]}{path}=$t[2];
		}
	close infile2;

	#-- Reading results
	
	my %fils=("COG","SQM/$project/results/07.$project.fun3.cog","KEGG","SQM/$project/results/07.$project.fun3.kegg");
	foreach my $tcase(sort keys %fils) {
		my $output=$fils{$tcase};
		open(inc,$output) || die "Cannot open $tcase output in $output\n";
		while(<inc>) {
			chomp;
			next if(!$_ || ($_=~/^\#/));
			my @k=split(/\t/,$_);
			$funannot{$tcase}{$k[1]}.="$k[0];"; 
			$accumc{$tcase}{$k[1]}++;
			}
		close inc;
		}
	foreach my $tcase(sort keys %fils) {
		my $outfile="SQM/$project/results/$tcase.summary";
		open(out,">$outfile") || die "Cannot open output file $outfile\n";
		print out "$tcase\tAbundance\tName\tFuction\tClass/Pathway\tORFs\n";
		foreach my $tfun(sort { $accumc{$tcase}{$b}<=>$accumc{$tcase}{$a}; } keys %{ $accumc{$tcase} }) {
			print out "$tfun\t$accumc{$tcase}{$tfun}\t$funs{$tcase}{$tfun}{name}\t$funs{$tcase}{$tfun}{fun}\t$funs{$tcase}{$tfun}{path}\t$funannot{$tcase}{$tfun}\n";
			}
		close outfile;
		print "$tcase summary created in $outfile\n";
		}
		
	#-- Summary for extdbs
	
	open(in,$opt_db) || die;
	while(<in>) {
		chomp;
		next if !$_;
		my($db,$dbfile,$dblist)=split(/\t/,$_);
	

		#-- Reading functions and pathways

		open(infile2,$dblist) || warn "Missing DB equivalence file $dblist\n";
		while(<infile2>) {
			chomp;
			next if(!$_ || ($_=~/\#/));
			my @t=split(/\t/,$_);
			$funs{$db}{$t[0]}{name}=$t[0];
			$funs{$db}{$t[0]}{fun}=$t[1];
			$funs{$db}{$t[0]}{path}=$t[2];
			}
		close infile2;
	
	
			#-- Reading results
	
			my $output="SQM/$project/results/07.$project.fun3.$db";
			open(inc,$output) || die "Cannot open $db output in $output\n";
			while(<inc>) {
				chomp;
				next if(!$_ || ($_=~/^\#/));
				my @k=split(/\t/,$_);
				$funannot{$db}{$k[1]}.="$k[0];"; 
				$accumc{$db}{$k[1]}++;
				}
			close inc;
			my $outfile="SQM/$project/results/$db.summary";
			open(out,">$outfile") || die "Cannot open output file $outfile\n";
			print out "$db\tAbundance\tName\tFuction\tClass/Pathway\tORFs\n";
			foreach my $tfun(sort { $accumc{$db}{$b}<=>$accumc{$db}{$a}; } keys %{ $accumc{$db} }) {
			print out "$tfun\t$accumc{$db}{$tfun}\t$funs{$db}{$tfun}{name}\t$funs{$db}{$tfun}{fun}\t$funs{$db}{$tfun}{path}\t$funannot{$db}{$tfun}\n";
				}
			close outfile;
			print "$db summary created in $outfile\n";

		}
	close in;									
	}
	

sub catch_error {
	system("rm $tempsample");
	die;
	}
		
