#!/usr/bin/perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Creates gene table putting together all the information from previous steps
#-- Modified 18/01/19 JT for working with new mapcount files

use strict;
use Tie::IxHash;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project path ok?"; }
do "$project/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($datapath,$resultpath,$coglist,$kegglist,$aafile,$ntfile,$gff_file,$rnafile,$fun3tax,$alllog,$nocog,$nokegg,$nopfam,$doublepass,$fun3kegg,$fun3cog,$fun3pfam,$opt_db,$fun3tax_blastx,$fun3kegg_blastx,$fun3cog_blastx,$gff_file_blastx,$fna_blastx,$mapcountfile,$mergedfile,$doublepass);

my $seqsinfile=0;     # Put sequences in the output table (0=no, 1=yes)

my(%orfdata,%contigdata,%cog,%kegg,%datafiles,%mapping);
tie %orfdata,"Tie::IxHash";
tie %mapping,"Tie::IxHash";


	#-- Reading gff for getting the list of ORFs

my $gff;
my %ingff;
if($doublepass) { $gff=$gff_file_blastx; } else { $gff=$gff_file; }
open(infile1,$gff) || die "Missing gff file $gff\n";
print "Reading GFF in $gff\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	if($_=~/ID\=([^;]+)/) { $ingff{$1}=1; }
	}
close infile1;



	#-- Reading data for COGs (names, pathways)

open(infile1,$coglist) || warn "Missing COG equivalence file\n";
print "Reading COG list\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$cog{$t[0]}{fun}=$t[1];
	$cog{$t[0]}{path}=$t[2]; 
            }
close infile1;

	#-- Reading data for KEGGs (names, pathways)

open(infile2,$kegglist) || warn "Missing KEGG equivalence file\n";
print "Reading KEGG list\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$kegg{$t[0]}{name}=$t[1];
	$kegg{$t[0]}{fun}=$t[2];
	$kegg{$t[0]}{path}=$t[3];
	}
close infile2;

	#-- Reading aa sequences 

open(infile3,$aafile) || die "I need the protein sequences from the prediction\n";
print "Reading aa sequences\n";
my($thisorf,$aaseq);
while(<infile3>) {
	chomp;
	if($_=~/^\>([^ ]+)/) {		#-- If we are reading a new ORF, store the data for the last one
		if($aaseq) { 
			$orfdata{$thisorf}{aaseq}=$aaseq; 
			$orfdata{$thisorf}{length}=(length $aaseq)+1; 
			$orfdata{$thisorf}{molecule}="CDS";
			$orfdata{$thisorf}{method}="Prodigal";
			}
		$thisorf=$1;
		$aaseq="";
		}
	else { $aaseq.=$_; }		#-- Otherwise store the sequence of the current	      
	}
close infile3;
if($aaseq) { 
	$orfdata{$thisorf}{aaseq}=$aaseq; 
	$orfdata{$thisorf}{length}=(length $aaseq)+1; 
	$orfdata{$thisorf}{molecule}="CDS";
	$orfdata{$thisorf}{method}="Prodigal";
	}


	#-- Reading RNAs

open(infile4,$rnafile) || warn "I need the RNA sequences from the prediction\n";
print "Reading RNA sequences\n";
my($thisrna,$rnaseq);
while(<infile4>) {
	chomp;
	if($_=~/^\>/) {			#-- If we are reading a new ORF, store the data for the last one
		$_=~s/^\>//;
		my @mt=split(/\t/,$_);
		if($rnaseq) { 
			$orfdata{$thisrna}{ntseq}=$rnaseq;
			$orfdata{$thisrna}{length}=(length $rnaseq)+1;
			$orfdata{$thisrna}{molecule}="RNA";
			$orfdata{$thisrna}{method}="barrnap";
			}
		$thisrna=$mt[0];
		my @l=split(/\s+/,$_,2);
		my @ll=split(/\;/,$l[1]);
		my $rnaname=$ll[0];
		$orfdata{$thisorf}{name}=$rnaname;  
		$rnaseq="";
		}
	else { $rnaseq.=$_; }		#-- Otherwise store the sequence of the current		      
}
close infile4;
if($rnaseq) { 
	$orfdata{$thisrna}{ntseq}=$rnaseq; 
	$orfdata{$thisrna}{length}=(length $rnaseq)+1;
	}

	#-- Reading taxonomic assignments

my $taxfile;
if($doublepass) { $taxfile="$fun3tax_blastx.wranks"; } else { $taxfile="$fun3tax.wranks"; }
open(infile5,$taxfile) || warn "Cannot open allorfs file $fun3tax.wranks\n";
print "Reading ORF information\n";
while(<infile5>) { 
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	my $mdat=$t[1];
	$mdat=~s/\;$//;
	$orfdata{$t[0]}{tax}=$mdat;
	$datafiles{'allorfs'}=1;
}
close infile5;

	#-- Reading nt sequences for calculating GC content

my($ntorf,$ntseq,$gc);
open(infile6,$ntfile) || warn "Cannot open nt file $ntfile\n";
print "Calculating GC content for genes\n";
while(<infile6>) { 
	chomp;
	if($_=~/^\>([^ ]+)/) {			#-- If we are reading a new ORF, store the data for the last one
		if($ntseq) { 
		$gc=gc_count($ntseq,$ntorf);
		$orfdata{$ntorf}{gc}=$gc;
		}
	$ntorf=$1;
	$ntseq="";
	}
	else { $ntseq.=$_; }		#-- Otherwise store the sequence of the current			      
}
close infile6;
if($ntseq) { $gc=gc_count($ntseq); }		#-- Last ORF in the file
$orfdata{$ntorf}{gc}=$gc; 
$datafiles{'gc'}=1;

	#-- Reading nt sequences for calculating GC content, from blastx

if($doublepass) {
	my($ntorf,$ntseq,$gc);
	open(infile6,$fna_blastx) || warn "Cannot open nt file $ntfile\n";
	print "Calculating GC content for blastx genes\n";
	while(<infile6>) { 
		chomp;
		if($_=~/^\>([^ ]+)/) {			#-- If we are reading a new ORF, store the data for the last one
			if($ntseq) { 
			$gc=gc_count($ntseq,$ntorf);
			my @sf=split(/\_/,$ntorf);
			my $ipos=pop @sf;
			my $contname=join("_",@sf);
			my($poinit,$poend)=split(/\-/,$ipos);
			$orfdata{$ntorf}{gc}=$gc;
			$orfdata{$ntorf}{length}=int(($poend-$poinit+1)/3);
			$orfdata{$ntorf}{molecule}="CDS";
			$orfdata{$ntorf}{method}="blastx";
			}
		$ntorf=$1;
		$ntseq="";
		}
		else { $ntseq.=$_; }		#-- Otherwise store the sequence of the current			      
	}
	close infile6;
	if($ntseq) { $gc=gc_count($ntseq); }		#-- Last ORF in the file
	$orfdata{$ntorf}{gc}=$gc; 
	$datafiles{'gc'}=1;
	}

	#-- Reading nt sequences for calculating GC content for RNAs

($ntorf,$ntseq,$gc)="";
open(infile7,$rnafile) || warn "Cannot open RNA file $rnafile\n";
print "Calculating GC content for RNAs\n";
while(<infile7>) { 
	chomp;
	if($_=~/^\>/) {			#-- If we are reading a new ORF, store the data for the last one
		$_=~s/^\>//;
		my @mt=split(/\t/,$_);
		if($ntseq) { 
			$gc=gc_count($ntseq,$ntorf);
			$orfdata{$ntorf}{gc}=$gc; 
		}
	$ntorf=$mt[0];
	$ntseq="";
                      }
	else { $ntseq.=$_; }		#-- Otherwise store the sequence of the current		      
}
close infile7;
if($ntseq) { $gc=gc_count($ntseq); }
$orfdata{$ntorf}{gc}=$gc; 

	#-- Reading taxonomic assignment and disparity for the contigs

open(infile8,$alllog) || warn "Cannot open contiglog file $alllog\n";
print "Reading contig information\n";
while(<infile8>) { 
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$contigdata{$t[0]}{tax}=$t[1]; 
	if($t[3]=~/Disparity\: (.*)/i) { $contigdata{$t[0]}{chimerism}=$1; }
	$datafiles{'alllog'}=1;
}
close infile8;

	#-- Reading KEGG annotations for the ORFs

if(!$nokegg) {
	if($doublepass) { $fun3kegg=$fun3kegg_blastx; }
	open(infile9,$fun3kegg) || warn "Cannot open fun3 KEGG annotation file $fun3kegg\n";
	print "Reading KEGG annotations\n";
	while(<infile9>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($gen,$f,$ko)=split(/\t/,$_);
		if($f) { 
			$orfdata{$gen}{kegg}=$f; 
			$orfdata{$gen}{name}=$kegg{$f}{name};	#-- Name of the gene (symbol), taken from KEGG
		}		
		if($ko) { $orfdata{$gen}{keggaver}=1; }	#-- Best aver must be the same than best hit, we just mark if there is best aver or not 
		$datafiles{'kegg'}=1;
	}
	close infile9;          
}
  
	#-- Reading COG annotations for the ORFs

if(!$nocog) {
	if($doublepass) { $fun3cog=$fun3cog_blastx; }
	open(infile10,$fun3cog) || warn "Cannot open fun3 COG annotation file $fun3cog\n";;
	print "Reading COGs annotations\n";
	while(<infile10>) { 
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($gen,$f,$co)=split(/\t/,$_);
		if($f) { $orfdata{$gen}{cog}=$f; }
		if($co) { $orfdata{$gen}{cogaver}=1; } #-- Best aver must be the same than best hit, we just mark if there is best aver or not
		$datafiles{'megancog'}=1;
	}
	close infile10;            
}
 
	#-- Reading Pfam annotations for the ORFs

if(!$nopfam) {
	open(infile11,$fun3pfam) || warn "Cannot open fun3 Pfam annotation file $fun3pfam\n";;
	print "Reading Pfam annotations\n";
	while(<infile11>) { 
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($gen,$co)=split(/\t/,$_);
		if($co) { $orfdata{$gen}{pfam}=$co; }
		$datafiles{'pfam'}=1;
	}
	close infile11; 
}           			       
  
	#-- Reading RPKM and coverage values for the ORFs in the different samples

open(infile12,$mapcountfile) || warn "Cannot open mapping file $mapcountfile\n";
print "Reading RPKMs and Coverages\n";
while(<infile12>) {
	chomp;
	next if(!$_ || ($_=~/\#/) || ($_=~/^Gen/));
	my($orf,$longg,$rawreads,$rawbases,$rpkm,$coverage,$idfile)=split(/\t/,$_);
	$mapping{$idfile}{$orf}{rpkm}=$rpkm;		#-- RPKM values
	$mapping{$idfile}{$orf}{raw}=$rawreads; 		#-- Raw counts
	$mapping{$idfile}{$orf}{coverage}=$coverage;	#-- Coverage values
	$mapping{$idfile}{$orf}{rawbases}=$rawbases;	#-- Coverage values
	#  print "$idfile*$orf*$fpkm\n"
}
close infile12;	     
  
	#-- CREATING GENE TABLE

print "Creating table\n";
open(outfile1,">$mergedfile") || die "I need an output file\n";

	#-- Headers

print outfile1 "#--Created by $0, ",scalar localtime,"\n";
print outfile1 "ORF\tCONTIG ID\tMOLECULE\tMETHOD\tLENGTH\tGC perc\tGENNAME\tTAX ORF\tKEGG ID\tKEGGFUN\tKEGGPATH\tCOG ID\tCOGFUN\tCOGPATH\tPFAM";
foreach my $cnt(sort keys %mapping) { print outfile1 "\tRPKM $cnt"; }
foreach my $cnt(sort keys %mapping) { print outfile1 "\tCOVERAGE $cnt"; }
foreach my $cnt(sort keys %mapping) { print outfile1 "\tRAW READ COUNT $cnt"; }
foreach my $cnt(sort keys %mapping) { print outfile1 "\tRAW BASE COUNT $cnt"; }
if($seqsinfile) { print outfile1 "\tAASEQ"; }
print outfile1 "\n";

	#-- ORF data
	
		#-- Sorting first by contig ID, then by position in contig

my (@listorfs,@sortedorfs);
foreach my $orf(keys %orfdata) {
	next if(!$ingff{$orf});		#-- Excluding ORFs not in gff table (removed in doublepass by overlapping hits in blastx)
	my @sf=split(/\_/,$orf);
	my $ipos=pop @sf;
	my $contname=join("_",@sf);
	my($poinit,$poend)=split(/\-/,$ipos);

	push(@listorfs,{'orf',=>$orf,'contig'=>$contname,'posinit'=>$poinit});
	}
@sortedorfs=sort {
	$a->{'contig'} cmp $b->{'contig'} ||
	$a->{'posinit'} <=> $b->{'posinit'}
	} @listorfs;

foreach my $orfm(@sortedorfs) { 
	my $orf=$orfm->{'orf'};
	my($cogprint,$keggprint);
	my $ctg=$orf;
	$ctg=~s/\_\d+\-\d+$//;
	# next if((!$mapping{'s22'}{$orf}{'fpkm'}) && (!$mapping{'st8'}{$orf}{'fpkm'})); 
	my $funcogm=$orfdata{$orf}{cog};
	my $funkeggm=$orfdata{$orf}{kegg};
	if($orfdata{$orf}{cogaver}) { $cogprint="$funcogm*"; } else { $cogprint="$funcogm"; }
	if($orfdata{$orf}{keggaver}) { $keggprint="$funkeggm*"; } else { $keggprint="$funkeggm"; }
	printf outfile1 "$orf\t$ctg\t$orfdata{$orf}{molecule}\t$orfdata{$orf}{method}\t$orfdata{$orf}{length}\t%.2f\t$orfdata{$orf}{name}\t$orfdata{$orf}{tax}\t$keggprint\t$kegg{$funkeggm}{fun}\t$kegg{$funkeggm}{path}\t$cogprint\t$cog{$funcogm}{fun}\t$cog{$funcogm}{path}\t$orfdata{$orf}{pfam}",$orfdata{$orf}{gc};

	#-- Abundance values

	foreach my $cnt(sort keys %mapping) { my $sdat=$mapping{$cnt}{$orf}{'rpkm'} || "0"; print outfile1 "\t$sdat"; }
	foreach my $cnt(sort keys %mapping) { my $sdat=$mapping{$cnt}{$orf}{'coverage'} || "0"; print outfile1 "\t$sdat"; }
	foreach my $cnt(sort keys %mapping) { my $sdat=$mapping{$cnt}{$orf}{'raw'} || "0"; print outfile1 "\t$sdat"; }
	foreach my $cnt(sort keys %mapping) { my $sdat=$mapping{$cnt}{$orf}{'rawbases'} || "0"; print outfile1 "\t$sdat"; }

	#-- aa sequences (if requested)

	if($seqsinfile) { print outfile1 "\t$orfdata{$orf}{aaseq}"; }
	print outfile1 "\n";
}
close outfile1;

print "============\nGENE TABLE CREATED: $mergedfile\n============\n\n";


#------------------- GC calculation

sub gc_count {
 my $seq=shift;
 my $corf=shift;
 my @m=($seq=~/G|C/gi);
 my $lseq=length $seq;
 if(!$lseq) { print "Zero length sequence found for $corf\n"; next; }
 my $gc=(($#m+1)/length $seq)*100;
 return $gc;
              }


