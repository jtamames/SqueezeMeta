#!/usr/bin/perl

#-- Creates STAMP profiles from sqm_(long)reads data
#-- (c) J Tamames (CSIC), Feb 2023, in BAE JCI, Antarctica

use strict;

my($header);
my $project=$ARGV[0];
my %rankequiv=('k','kingdom','p','phylum','c','class','o','order','f','family','g','genus','s','species');
my %allsamples;

keggcog($project);
tax($project);
metadata($project);

sub metadata {
	my $outfile="$project/$project.met";
	open(out,">$outfile") || die "Cannot open $outfile\n";
	print out "Sample\tGroup\n";
	foreach my $sam(sort keys %allsamples) {
		print out "$sam\t$sam\n";
		}
	close out;
	print "Metadata file created in $outfile\n";
	}
	
	

sub keggcog {
	my(@fils,@allfuns);
	opendir(indir,$project) || die "Cannot open directory $project\n";
	@fils=grep(/out\.allreads\.fun/,readdir(indir));
	closedir(indir);
	foreach my $allfun(@fils) {
		my @e=split(/\./,$allfun);
		my $fclass=$e[$#e];
		$fclass=~s/^fun//;
		push(@allfuns,$fclass);
		}
	foreach my $func(@allfuns) {
		my $infile="$project/$project.out.allreads.fun$func";
		my $outfile="$project/$project.$func.spf";
		my $lastdat;
		if(($func eq "kegg") || ($func eq "cog")) { $lastdat=2; } else { $lastdat=1; }
		open(in,$infile) || die "Cannot open $infile\n";
		open(out,">$outfile") || die "Cannot open $outfile\n";
		$header="";
		while(<in>) {	
			chomp;
			next if(!$_ || ($_=~/^\#/));
			if(!$header) {
				$header=$_;
				my @header=split(/\t/,$_);
				print out "$func ID";
				for(my $pos=2; $pos<=($#header)-$lastdat; $pos++) {
					print out "\t$header[$pos]";
					$allsamples{$header[$pos]}=1;
				}
				print out "\n";	
				next;
				}
			my @dat=split(/\t/,$_);	
			my $fun=$dat[$#dat-$lastdat+1];
			my $funclass;
			if($lastdat>1) { $funclass=$dat[$#dat]; }
			my @aclasses=split(/\|/,$funclass);
			next if($fun!~/[a-z]/);
			print out "$dat[0]:$fun";
			for(my $pos=2; $pos<=$#dat-$lastdat; $pos++) {
				print out "\t$dat[$pos]";
				}
			print out "\n";
			}
		close(in);
		close out;
		print "$func STAMP profile created in $outfile\n";
		}
	}	

sub tax {
	my $infile="$project/$project.out.allreads.mcount";
	foreach my $tcount('reads','orfs') {
		my $header;
		my @header;
		my $outfile="$project/$project.tax.$tcount.spf";
		open(in,$infile) || die "Cannot open $infile\n";
		open(out,">$outfile") || die "Cannot open $outfile\n";
		while(<in>) {	
			chomp;
			next if(!$_ || ($_=~/^\#/));
			if(!$header) {
				$header=$_;
				@header=split(/\t/,$header);
				print out "Rank\tTaxon";
				for(my $pos=4; $pos<=$#header; $pos++) {
					my $tlabel=$header[$pos];
					if($tlabel=~/$tcount/i) {
						$tlabel=~s/ reads| ORFs//;
						print out "\t$tlabel";
						}
					}
				print out "\n";	
				next;
				}
			my @dat=split(/\t/,$_);	
			my $rank=$rankequiv{$dat[0]};
			my @tx=split(/\;/,$dat[1]);
			my $taxon=$tx[$#tx];
			# $taxon=~s/..//;
			print out "$rank\t$taxon";
			for(my $pos=4; $pos<=$#dat; $pos++) {
				if($header[$pos]=~/$tcount/i) { print out "\t$dat[$pos]"; }
				}
			print out "\n";
			}
		close(in);
		close out;
		print "Tax $tcount STAMP profile created in $outfile\n";
		}	
	}
