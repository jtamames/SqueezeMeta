#!/usr/bin/perl

# (c) Javier Tamames, CNB-CSIC

#-- This program prepares sqm output for being loaded in anvio

use strict;
use Cwd;
use lib "."; 

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//;
if(!$project) { die "Usage: sqm2anvio.pl <project name>\n"; }

if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in project $project. Please check project name\n"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($resultpath,$datapath,$gff_file,$gff_file_blastx,$mergedfile,$contigsfna);

my $version="1.0";
my $gff;
my(%genstore,%genindex);

if(-e $gff_file_blastx) { $gff=$gff_file_blastx; }
elsif(-e $gff_file) { $gff=$gff_file; }
else { die "Cannot open gff file\n"; }
my $genes_out="$project\_anvio_genes.txt";
my $equivalence_out="$project\_anvio_SQMequivalence.txt";
my $contigs_out="$project\_anvio_contigs.txt";
my $functions_out="$project\_anvio_functions.txt";
my $taxonomy_out="$project\_anvio_taxonomy.txt";

open(infile1,$gff) || die "Cannot open gff file $gff\n";
open(outfile1,">$genes_out") || die "Cannot open gen outfile $genes_out\n";
open(outfile2,">$equivalence_out") || die "Cannot open equivalence outfile $equivalence_out\n";
print outfile1 "gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tsource\tversion\n";
print outfile2 "anvio_ID\tSQM_ID\n";
my $geneidx;
while(<infile1>) {
	chomp;	
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $contigid=$k[0];
	my $source=$k[1];
	my $ntlong=$k[4]-$k[3]+1;
	my $initgen=$k[3]-1;
	my $endgen=$k[4];	#-- This is incorrect and weird, but anvio has its own way to index genes and wants it this way
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
	print outfile1 "$genindex{$geneid}\t$contigid\t$initgen\t$endgen\t$direction\t$partial\t$source\t$version\n";
	}
close outfile1;
close outfile2;
close infile1;

my $header;
my @head;
open(outfile3,">$functions_out") || die "Cannot open functions outfile $functions_out\n";
print outfile3 "gene_callers_id\tsource\taccesion\tfunction\te_value\n";
open(outfile4,">$taxonomy_out") || die "Cannot open taxonomy outfile $taxonomy_out\n";
print outfile4 "gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n";
open(infile2,$mergedfile)  || die "Cannot open gene table $mergedfile\n"; 
while(<infile2>) {
	chomp;	
	next if(!$_ || ($_=~/^\#/));
	if(!$header) { $header=$_; @head=split(/\t/,$header); next; }
	my @k=split(/\t/,$_);
	foreach(my $pos=0; $pos<=$#k; $pos++) {
		my $f=$head[$pos];
		$k[$pos]=~s/\*$//;
		if(($f eq "KEGG ID") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tKEGG\t$k[$pos]\t$k[$pos+1]\t0\n"; }
		if(($f eq "COG ID") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tCOG\t$k[$pos]\t$k[$pos+1]\t0\n"; }
		if(($f eq "PFAM") && $k[$pos]) { print outfile3 "$genindex{$k[0]}\tPFAM\t$k[$pos]\t\t0\n"; }
		my($kingdom,$phylum,$class,$order,$family,$genus,$species);
		if(($f eq "TAX ORF") && $k[$pos]) { 
			if($k[$pos]=~/k_([^;]+)/) { $kingdom=$1; }
			if($k[$pos]=~/p_([^;]+)/) { $phylum=$1; }
			if($k[$pos]=~/o_([^;]+)/) { $order=$1; }
			if($k[$pos]=~/c_([^;]+)/) { $class=$1; }
			if($k[$pos]=~/f_([^;]+)/) { $family=$1; }
			if($k[$pos]=~/g_([^;]+)/) { $genus=$1; }
			if($k[$pos]=~/s_([^;]+)/) { $species=$1; }		
			print outfile4 "$genindex{$k[0]}\t$kingdom\t$phylum\t$class\t$order\t$family\t$genus\t$species\n";
			}
		}
	}
close infile2;
close outfile3;
close outfile4;

open(infile3,$contigsfna) || die "Cannot open contig file $contigsfna\n";
open(outfile5,">$contigs_out") || die "Cannot open output in $contigs_out\n";
while(<infile3>) {
	if($_=~/^\>([^ ]+)/) { print outfile5 ">$1\n"; } else { print outfile5 $_; }
	}
close infile3;
close outfile5;

system("tar cf anvio\_$project.tar $genes_out $contigs_out $functions_out $taxonomy_out $equivalence_out");
system("gzip anvio\_$project.tar");
system("rm $genes_out $contigs_out $functions_out $taxonomy_out");
 
print "Created anvio gene file (import with anvi-gene-contigs-database): $genes_out\n";
print "Created anvio contig file (import with anvi-gene-contigs-database): $contigs_out\n";
print "Created anvio functions file (import with anvi-import-functions): $functions_out\n";
print "Created anvio taxonomy file (import with anvi-import-taxonomy-for-genes): $taxonomy_out\n";
print "Created anvio-SQM equivalence (just for reference, don't import this): $equivalence_out\n";

print "\n============= TAR file =============\nCreated tar.gz file: anvio\_$project.tar.gz\n\n";	
		
	
