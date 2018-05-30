#!/usr/bin/perl

$|=1;

use strict;

my $file="$databasedir/uniprot_trembl.dat.gz";

my $command="wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz -nc -P $databasedir";
# system $command;


my($id,$seq,$acc,$inseq,$M);
my @kid;
open(infile1,"zcat $file|") || die "Cannot open $file\n";
while(<infile1>) {
	chomp;
	next if !$_;
	if($_=~/^AC\s+(.*)/) { $id=$1; chop $id; }
	if($_=~/^DR\s+KO/) { 
		my @f=split(/\; /,$_);
		$acc=$f[1];
		next if($acc!~/\w/);
		push(@kid,$acc);
		$m=1;
		}
	if($_=~/^SQ/) { $inseq=1; next; }
	if($_=~/\/\//) {
		if($m) {
			my $kegg=join(";",sort @kid);
			$seq=~s/\s+//g;
			print ">$id|$kegg\n$seq\n";
			}
		($inseq,$m)=0;
		@kid=();
		$seq="";
		}
	if($inseq) { $seq.=$_; }
	}   		     
close infile1;
