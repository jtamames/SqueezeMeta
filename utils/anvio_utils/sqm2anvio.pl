#!/usr/bin/env perl

# (c) Javier Tamames, CNB-CSIC

#-- This program prepares sqm output for being loaded in anvio

use strict;
use Cwd;
use lib ".";

use File::Basename;
use Cwd 'abs_path';

$|=1;

my $pwd=cwd();
my $projectpath=$ARGV[0];
my $outdir=$ARGV[1];
my $anvio_version=$ARGV[2];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";


if((!$project) or (!$outdir)) { die "Usage: sqm2anvio.pl <project name> <output dir> <anvio version>\n"; }


#-- Configuration variables from conf file

our($installdir,$resultpath,$datapath,$gff_file,$gff_file_blastx,$mergedfile,$contigsfna,$contigsinbins,$datapath);

our $scriptdir = abs_path(dirname(__FILE__));
our $installpath = "$scriptdir/../..";

my $version="1.2";
my $gff;
my(%genstore,%genindex);

my $project_name = ( split(/\//,$project) )[-1];

print "======== sqm2anvio.pl, v$version ========\n";

my $ecode = system("mkdir $outdir 2> /dev/null");
if($ecode!=0) { die("Directory $outdir already exists or can't be created"); }
system("mkdir $outdir/bam");

if(-e $gff_file_blastx) { $gff=$gff_file_blastx; }
elsif(-e $gff_file) { $gff=$gff_file; }
else { die "Cannot open gff file\n"; }

my $genes_out="$project_name\_anvio_genes.txt";
my $equivalence_out="$project_name\_anvio_SQMequivalence.txt";
my $contigs_out="$project_name\_anvio_contigs.txt";
my $functions_out="$project_name\_anvio_functions.txt";
my $taxonomy_out="$project_name\_anvio_taxonomy.txt";
my $bin_out="$project_name\_anvio_bins.txt";

open(infile1,$gff) || die "Cannot open gff file $gff\n";
open(outfile1,">$outdir/$genes_out") || die "Cannot open gen outfile $outdir/$genes_out\n";
open(outfile2,">$outdir/$equivalence_out") || die "Cannot open equivalence outfile $outdir/$equivalence_out\n";
if ($anvio_version >= 7){ print outfile1 "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\n"; }
else { print outfile1 "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n"; }
print outfile2 "anvio_ID\tSQM_ID\n";
my $geneidx;
while(<infile1>) {
	chomp;	
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $contigid=$k[0];
	my $source=$k[1];
	my $type=$k[2];
	my $ntlong=$k[4]-$k[3]+1;
	my $initgen=$k[3]-1;
	my $endgen=$k[4];	#-- Not -1 as it should be. This is weird and makes no sense, but anvio has its own way to index genes and wants it this way
	my $direction=$k[6];
	if($direction eq "+") { $direction="f"; }
	elsif($direction eq "-") { $direction="r"; }
	else { $direction="f"; }
	my($geneid,$partial);
	if($k[8]=~/ID\=([^;]+)/) { 
		$geneid=$1;
		$geneidx++;
		$genindex{$geneid}=$geneidx;
		print outfile2 "$geneidx\t$geneid\n";
		}
	if($k[8]=~/partial\=00/) { $partial=0; } else { $partial=1; } 
	my $call_type;
	if($type eq 'CDS') { $call_type=1; }
	elsif($type=~/RNA/) { $call_type=2; }
	else { $call_type=3; } # Unknown (should not really be happening)
	if ($anvio_version >= 7) { print outfile1 "$genindex{$geneid}\t$contigid\t$initgen\t$endgen\t$direction\t$partial\t$call_type\t$source\t$version\n"; }
	else { print outfile1 "$genindex{$geneid}\t$contigid\t$initgen\t$endgen\t$direction\t$partial\t$source\t$version\n"; }
	}
close outfile1;
close outfile2;
close infile1;


my %orftax;
my $firstline;

# get sqm2tables.py taxonomy
my $ecode = system("$installpath/utils/sqm2tables.py $projectpath $outdir --sqm2anvio");
if($ecode!=0) { die("Error running subprocess: $installpath/utils/sqm2tables.py $project $outdir --sqm2anvio"); }
open(infile2, "$outdir/$project_name.orf.tax.allfilter.tsv") || die "Cannot open taxonomy table $outdir/$project_name.orf.tax.allfilter.tsv";
while(<infile2>) {
	chomp;
	if(!$firstline) { $firstline=$_ ; next;}
	my ($orf, $tax) = split(/\t/, $_, 2);
       	$orftax{$orf} = $tax;
	}

my $header;
my @head;
open(outfile3,">$outdir/$functions_out") || die "Cannot open functions outfile $outdir/$functions_out\n";
print outfile3 "gene_callers_id\tsource\taccession\tfunction\te_value\n";
open(outfile4,">$outdir/$taxonomy_out") || die "Cannot open taxonomy outfile $outdir/$taxonomy_out\n";
print outfile4 "gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n";
open(infile3,$mergedfile)  || die "Cannot open gene table $mergedfile\n";
while(<infile3>) {
	chomp;	
	next if(!$_ || ($_=~/^\#/));
	if(!$header) { $header=$_; @head=split(/\t/,$header); next; }
	my @k=split(/\t/,$_);
	foreach(my $pos=0; $pos<=$#k; $pos++) {
		my $f=$head[$pos];
		$k[$pos]=~s/\*$//;
		if(($f eq "KEGG ID") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tKEGG\t$k[$pos]\t$k[$pos+1]\t0\n"; }
		if(($f eq "KEGGPATH") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tKEGGPATH\t\t$k[$pos]\t0\n"; }
		if(($f eq "COG ID") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tCOG\t$k[$pos]\t$k[$pos+1]\t0\n"; }
                if(($f eq "COGPATH") && $k[$pos]) { print outfile3  "$genindex{$k[0]}\tCOGPATH\t\t$k[$pos]\t0\n"; }
		if(($f eq "PFAM") && $k[$pos]) { 
			my($pfid,$pffun)=split(/\s+/,$k[$pos], 2);
			$pffun=~s/\[|\]//g;
			print outfile3 "$genindex{$k[0]}\tPFAM\t$pfid\t$pffun\t0\n"; 
			}
		}
	print outfile4 "$genindex{$k[0]}\t$orftax{$k[0]}\n";
#		my($kingdom,$phylum,$class,$order,$family,$genus,$species);
#		if(($f eq "TAX ORF") && $k[$pos]) { 
#			if($k[$pos]=~/k_([^;]+)/) { $kingdom=$1; }
#			if($k[$pos]=~/p_([^;]+)/) { $phylum=$1; }
#			if($k[$pos]=~/o_([^;]+)/) { $order=$1; }
#			if($k[$pos]=~/c_([^;]+)/) { $class=$1; }
#			if($k[$pos]=~/f_([^;]+)/) { $family=$1; }
#			if($k[$pos]=~/g_([^;]+)/) { $genus=$1; }
#			if($k[$pos]=~/s_([^;]+)/) { $species=$1; }		
#			print outfile4 "$genindex{$k[0]}\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
#			}
#		}
	}
close infile3;
close outfile3;
close outfile4;

system("rm $outdir/$project_name.orf.tax.allfilter.tsv"); # we don't need this anymore
system("mv $outdir/$project_name.contig.tax.allfilter.tsv $outdir/$project_name\_anvio_contig_taxonomy.txt"); # so Natalia is happy

open(infile4,$contigsfna) || die "Cannot open contig file $contigsfna\n";
open(outfile5,">$outdir/$contigs_out") || die "Cannot open output in $outdir/$contigs_out\n";
while(<infile4>) {
        chomp;
	if($_=~/^\>([^ ]+)/) { print outfile5 ">$1\n"; } else { print outfile5 "$_\n"; }
	}
close infile4;
close outfile5;

if(-e $contigsinbins) {
	open(infile5,$contigsinbins) || warn "Cannot open bin file $contigsinbins\n"; 
	open(outfile6,">$outdir/$bin_out") || die "Cannot open taxonomy outfile $outdir/$bin_out\n";
	while(<infile5>) {
		chomp;
		next if($_=~/^\#/);
		my @k=split(/\t/,$_);
		$k[2]=~s/\./\_/g;
		print outfile6 "$k[0]\t$k[2]\n";
		}
	close infile5;
	close outfile6;
	}
else { print "Cannot find contigs in bins $contigsinbins file. No binning?\n"; }

my $samdir="$datapath/sam";
my $samkeep;
opendir(indir1,$samdir) || die;
my @samfiles=grep(/\.sam$/,readdir indir1);
my $samlist=join(" ",@samfiles);
closedir indir1;
#if($#samfiles>=0) { 
#	print "SAM files found for this run ($samlist)\nDo you want to compress them and include them in the output folder (y/n)? ";
#        while(1) {
#		$samkeep=<STDIN>;
#		chomp $samkeep;
#		if($samkeep eq 'y' or $samkeep eq 'yes'){ $samkeep=1; print "Compressing SAM files to the BAM format\n"; last }
#		elsif($samkeep eq 'n' or $samkeep eq 'no') { $samkeep=0; print "SAM files will be ignored\n"; last }
#               else { print "Only y(es) or n(o) are valid answers\n" }
#		}
#	}

$samkeep=1;
if($samkeep) { 
	foreach my $sam(@samfiles) {
		(my $bam = $sam) =~ s/\.sam$/-RAW.bam/;
                $bam =~ s/$project\.//;
        	my $ecode = system("LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$installpath/lib $installpath/bin/samtools view -b $samdir/$sam > $outdir/bam/$bam");
		if($ecode!=0) { die "Error running samtools"; }
		}
	}

print "Created anvio gene file (import with anvi-gene-contigs-database): $outdir/$genes_out\n";
print "Created anvio contig file (import with anvi-gene-contigs-database): $outdir/$contigs_out\n";
print "Created anvio functions file (import with anvi-import-functions): $outdir/$functions_out\n";
print "Created anvio taxonomy file (import with anvi-import-taxonomy-for-genes): $outdir/$taxonomy_out\n";
if(-e $bin_out) { print "Created anvio bins file (import with anvi-import-collection): $outdir/$bin_out\n"; }
print "Created anvio-SQM equivalence (just for reference, don't import this): $outdir/$equivalence_out\n";
if($samkeep) { print "Added SAM files: $samlist\n and converted them to the BAM format\n"; }

