#!/usr/bin/env perl

#-- Part of SqueezeMeta distribution. 01/05/2018 Original version, (c) Javier Tamames, CNB-CSIC
#-- Counts the coverage of all functions
#-- Modified 18/01/19 JT for working with new mapcount files

use strict;
use Tie::IxHash;
use Cwd;
use lib ".";

$|=1;

my $pwd=cwd();
my $taxreq=$ARGV[1];	#-- Invoke it with a name of taxon to get just functions for that taxon

my $projectpath=$ARGV[0];
if(!$projectpath) { die "Please provide a valid project name or project path\n"; }
if(-s "$projectpath/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $projectpath. Is the project path ok?"; }
do "$projectpath/SqueezeMeta_conf.pl";
our($projectname);
my $project=$projectname;

do "$projectpath/parameters.pl";

#-- Checking for version compatibility

our($installpath);
system("perl $installpath/utils/versionchange.pl $projectpath");
do "$projectpath/SqueezeMeta_conf.pl";

#-- Configuration variables from conf file

our($installpath,$datapath,$resultpath,$extpath,$kegglist,$coglist,$ntfile,$fun3tax,$fun3kegg,$fun3cog,$fun3tax_blastx,$fun3kegg_blastx,$fun3cog_blastx,$opt_db,$nokegg,$nocog,$mapcountfile,$doublepass,$minraw12,$syslogfile);

print "  Calculating coverage for functions\n";
my(%funs,%taxf,%validid,%tfun,%totalbases,%totalreads,%allsamples,%funstat,%longorfs,%taxcount,%optdb);

	#-- Reading KEGG functions and pathways

open(outsyslog,">>$syslogfile") || warn "Cannot open syslog file $syslogfile for writing the program log\n";

open(infile1,$kegglist) || warn "Missing KEGG equivalence file $kegglist\n";
while(<infile1>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$funs{kegg}{$t[0]}{name}=$t[1];
	$funs{kegg}{$t[0]}{fun}=$t[2];
	$funs{kegg}{$t[0]}{path}=$t[3];
	}
close infile1;

	#-- Reading COG functions and pathways

open(infile2,$coglist) || warn "Missing COG equivalence file $coglist\n";
while(<infile2>) {
	chomp;
	next if(!$_ || ($_=~/\#/));
	my @t=split(/\t/,$_);
	$funs{cog}{$t[0]}{fun}=$t[1];
	$funs{cog}{$t[0]}{path}=$t[2];
	}
close infile2;


	#-- Reading OPT_DB functions and names

if($opt_db) {
	open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
	while(<infile0>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($dbname,$extdb,$listf)=split(/\t/,$_);
		if(-e $listf) {
			open(infile3,$listf) || warn "Can't open names file for $opt_db\n";
			while(<infile3>) {
				chomp;
				next if(!$_ || ($_=~/\#/));
				my @t=split(/\t/,$_);
				$funs{$dbname}{$t[0]}{fun}=$t[1];
				$funs{$dbname}{$t[0]}{path}=$t[1];
				# print "*$dbname*$t[0]*path*$funs{$dbname}{$t[0]}{path}*\n";
				}
			close infile3;
			}
		}
	close infile0;	
	}


my %equival;
tie %equival,"Tie::IxHash";
%equival=('k','k','p','p','c','c','o','o','f','f','g','g','s','s');

	#-- Read the taxonomic assignment for each gene

my $taxfile;
if($doublepass) { $taxfile="$fun3tax_blastx.wranks"; } else { $taxfile="$fun3tax.wranks"; }
open(infile3,"$taxfile") || warn "Can't open wranks file $taxfile\n";
while(<infile3>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my @t=split(/\;/,$k[1]);
	
	#-- Loop for all ranks for the gene
	
	foreach my $cm(@t) {
		my ($rank,$ttax)=split(/\_/,$cm);
		my $cortax=$equival{$rank};
		next if(!$cortax);
		$taxf{$k[0]}{$cortax}=$ttax;
		if($taxreq && ($ttax eq $taxreq)) { $validid{$k[0]}=1; }
		# print "$k[0] $cortax $ttax\n";
		}
	}
close infile3;


	#-- Reading KEGG assignments

if(!$nokegg) {
	if($doublepass) { $fun3kegg=$fun3kegg_blastx; }
	open(infile4,$fun3kegg) || die "Can't open KEGG assignments in $fun3kegg\n";;
	while(<infile4>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		$tfun{$k[0]}{kegg}=$k[1];
		}
	close infile4;
	}
	
	#-- Reading COG assignments

if(!$nocog) {
	if($doublepass) { $fun3cog=$fun3cog_blastx; }
	open(infile5,$fun3cog) || die "Can't open COG assignments in $fun3cog\n";
	while(<infile5>) {
		chomp;
		next if(!$_ || ($_=~/^\#/));
		my @k=split(/\t/,$_);
		$tfun{$k[0]}{cog}=$k[1];
		}
	close infile5;
	}

	
	#-- Reading OPTDB assignments

if($opt_db) {
	open(infile0,$opt_db) || warn "Can't open EXTDB file $opt_db\n"; 
	while(<infile0>) {
		chomp;
		next if(!$_ || ($_=~/\#/));
		my($dbname,$extdb,$dblist)=split(/\t/,$_);
		$optdb{$dbname}=1;
		my $fun3opt="$resultpath/07.$project.fun3.$dbname";
		if($doublepass) { $fun3opt="$resultpath/08.$project.fun3.$dbname"; }
		open(infile10,$fun3opt) || warn "Can't open fun3 $dbname annotation file $fun3opt\n";
		while(<infile10>) { 
			chomp;
			next if(!$_ || ($_=~/\#/));
			my @k=split(/\t/,$_);
			$tfun{$k[0]}{$dbname}=$k[1];
			}
		close infile10;  
		}
	close infile0;          
}

	
	#-- Reading coverages for all genes

print "  Reading coverage in $mapcountfile\n";
open(infile6,$mapcountfile) || warn "Can't open coverages in $mapcountfile\n";
while(<infile6>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $cfun_kegg=$tfun{$k[0]}{kegg};	#-- Corresponding KEGG for this gene
	my $cfun_cog=$tfun{$k[0]}{cog}; 	#-- Corresponding COG for this gene
	# foreach my $odb(sort keys %optdb) { my $cfun_opt{$k[0]}{$odb}=$tfun{$k[0]}{$o};
	my $sample=$k[$#k];
	$longorfs{$k[0]}=$k[1];		#-- Length of the gene
	my $mapbases=$k[3];
	# print "$k[0] $cfun_kegg $cfun_cog $sample $longorfs{$k[0]}\n";
	$totalbases{$sample}+=$mapbases;
	next if((!$cfun_kegg) && (!$cfun_cog));
	next if($taxreq && (!$validid{$k[0]}));
	$allsamples{$sample}++;
	if($mapbases) {
	
		#-- Counting KEGGs
	
		if($cfun_kegg) {
			my @kegglist=split(/\;/,$cfun_kegg);	#-- Support for multiple COGS (in annotations such as COG0001;COG0002, all COGs get the counts)
			foreach my $tlist_kegg(@kegglist) {
				$funstat{kegg}{$tlist_kegg}{$sample}{copies}++;
				$funstat{kegg}{$tlist_kegg}{$sample}{length}+=$longorfs{$k[0]}; 
				$funstat{kegg}{$tlist_kegg}{$sample}{bases}+=$mapbases;
				foreach my $tk(keys %equival) {
					my $krank=$equival{$tk};
					my $itax=$taxf{$k[0]}{$krank};
					if($itax) { $taxcount{kegg}{$tlist_kegg}{$sample}{$krank}{$itax}++; }
					}
				}
			}
	
		#-- Counting OPT_DB
	
		foreach my $odb(sort keys %optdb) {
			my $cfun_opt=$tfun{$k[0]}{$odb};
			if($cfun_opt) {
				my @optlist=split(/\;/,$cfun_opt);	#-- Support for multiple COGS (in annotations such as COG0001;COG0002, all COGs get the counts)
				foreach my $tlist_opt(@optlist) {
					$funstat{$odb}{$tlist_opt}{$sample}{copies}++;
					$funstat{$odb}{$tlist_opt}{$sample}{length}+=$longorfs{$k[0]}; 
					$funstat{$odb}{$tlist_opt}{$sample}{bases}+=$mapbases;
					foreach my $tk(keys %equival) {
						my $krank=$equival{$tk};
						my $itax=$taxf{$k[0]}{$krank};
						if($itax) { $taxcount{$odb}{$tlist_opt}{$sample}{$krank}{$itax}++; }
						}
					}
				}
			}

		#-- Counting COGs
	
	if($cfun_cog) { 
		my @coglist=split(/\;/,$cfun_cog);	#-- Support for multiple COGS (in annotations such as COG0001;COG0002, all COGs get the counts)
		foreach my $tlist_cog(@coglist) {
			$funstat{cog}{$tlist_cog}{$sample}{copies}++;
			$funstat{cog}{$tlist_cog}{$sample}{length}+=$longorfs{$k[0]}; #-- Para leer las longitudes directamente de las secuencias (para ficheros coverage antiguos que no la tienen)
			$funstat{cog}{$tlist_cog}{$sample}{bases}+=$mapbases;
			# print "$k[0]*$cfun_cog*$sample*$longorfs{$k[0]}*$funstat{cog}{$cfun_cog}{$sample}{length}\n";
			foreach my $tk(keys %equival) {
				my $krank=$equival{$tk};
				my $itax=$taxf{$k[0]}{$krank};
				if($itax) { $taxcount{cog}{$tlist_cog}{$sample}{$krank}{$itax}++; }
 				}
			}
		}
	}
}
close infile6;


	#-- Reading RPKMs for all genes

print "  Reading rpkm in $mapcountfile\n";
open(infile7,$mapcountfile) || die "Can't open $mapcountfile\n";
while(<infile7>) {
	chomp;
	next if(!$_ || ($_=~/^\#/));
	my @k=split(/\t/,$_);
	my $cfun_kegg=$tfun{$k[0]}{kegg}; 
	my $cfun_cog=$tfun{$k[0]}{cog}; 
	my $sample=$k[$#k];
	$totalreads{$sample}+=$k[2];
	next if((!$cfun_kegg) && (!$cfun_cog));
	next if($taxreq && (!$validid{$k[0]}));
	if($k[2]) {
		my @kegglist=split(/\;/,$cfun_kegg);
		my @coglist=split(/\;/,$cfun_cog);
		foreach my $tlist_kegg(@kegglist) { 
			if($tlist_kegg) { $funstat{kegg}{$tlist_kegg}{$sample}{reads}+=$k[2]; }
			}
		foreach my $tlist_cog(@coglist) { 
			if($tlist_cog) { $funstat{cog}{$tlist_cog}{$sample}{reads}+=$k[2]; }  
			}
		foreach my $odb(sort keys %optdb) {
			my $cfun_opt=$tfun{$k[0]}{$odb};
			if($cfun_opt) { $funstat{$odb}{$cfun_opt}{$sample}{reads}+=$k[2]; }
			}
		}
	
	}
close infile7;

	#-- Creating output files

my $rawf;
my %rpk;
foreach my $classfun(sort keys %funstat) {
	$rawf="$resultpath/12.$project.$classfun.funcover";
	print "  Now creating $classfun coverage output in $rawf\n";
	print outsyslog "Creating $classfun coverage output in $rawf\n";
	open(outfile1,">$rawf") || die "Can't open $rawf for writing\n";
	print outfile1 "#-- Created by $0 from $mapcountfile, ",scalar localtime;
	if($taxreq) { print outfile1 ", for taxon $taxreq"; }
	print outfile1 "\n";
	# print outfile1 "# $classfun ID\tSample\tCopy number\tTotal length\tTotal bases\tCoverage\tNorm Coverage\tNorm Coverage per copy\tTPM\tDistribution\tName\tFunction\n";
	print outfile1 "# $classfun ID\tSample\tCopy number\tTotal length\tTotal bases\tCoverage\tTPM\tDistribution";
	if($classfun ne "cog") { print outfile1 "\tName"; }
	print outfile1 "\tFunction\n";
	
	#-- Calculation of coverage, norm coverage, and RPKM

	foreach my $samp(sort keys %allsamples) {
		my $accumrpk;
		foreach my $kid(sort keys %{ $funstat{$classfun} }) {		#-- For TPM calculation
			next if(!$funstat{$classfun}{$kid}{$samp}{length}); 
			my $longt=$funstat{$classfun}{$kid}{$samp}{length};
                        next if(!$funstat{$classfun}{$kid}{$samp}{copies});
                        my $copies=$funstat{$classfun}{$kid}{$samp}{copies};
                        my $avglongt=$longt/$copies;
			$rpk{$kid}=$funstat{$classfun}{$kid}{$samp}{reads}/$avglongt;
			$accumrpk+=$rpk{$kid}; 
			}
		$accumrpk/=1000000;
			
		foreach my $kid(sort keys %{ $funstat{$classfun} }) {
			next if(!$funstat{$classfun}{$kid}{$samp}{length}); 
			my $cover=$funstat{$classfun}{$kid}{$samp}{bases}/$funstat{$classfun}{$kid}{$samp}{length};
			my $normcover=(($funstat{$classfun}{$kid}{$samp}{bases}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalbases{$samp}));
			my $normcoverpercopy=$normcover*$funstat{$classfun}{$kid}{$samp}{copies};
			# print "$classfun*$kid*$samp*$funstat{$classfun}{$kid}{$samp}{length}*$totalreads{$samp}\n";
			my $rpkm=(($funstat{$classfun}{$kid}{$samp}{reads}*1000000000)/($funstat{$classfun}{$kid}{$samp}{length}*$totalreads{$samp}));
 			my $tpm=$rpk{$kid}/$accumrpk;
			my $stringtax=""; 
			foreach my $tk(keys %equival) {
				my $krank=$equival{$tk};
				my $countt=0;
				foreach my $tt(keys %{ $taxcount{$classfun}{$kid}{$samp}{$krank} }) { $countt++; }
				$stringtax.="$krank:$countt;";
				}
			chop $stringtax;
			printf outfile1 "$kid\t$samp\t$funstat{$classfun}{$kid}{$samp}{copies}\t$funstat{$classfun}{$kid}{$samp}{length}\t$funstat{$classfun}{$kid}{$samp}{bases}\t%.3f\t%.3f\t$stringtax",$cover,$tpm; 
 			if($classfun ne "cog") { print outfile1 "\t$funs{$classfun}{$kid}{name}"; }
			print outfile1 "\t$funs{$classfun}{$kid}{fun}\n";
			}
		}
close outfile1;
	}

foreach my $classfun(sort keys %funstat) {
	$rawf="$extpath/12.$project.$classfun.stamp";        #-- Creating STAMP files
	print outsyslog "Creating $classfun raw reads output in $rawf\n";
	print "  Now creating $classfun raw reads output in $rawf\n";
	open(outfile2,">$rawf") || die "Can't open $rawf for writing\n";
	if($classfun eq "cog") { print outfile2 "$classfun class\t$classfun ID"; }
        else { print outfile2 "$classfun ID"; }
foreach my $samp(sort keys %allsamples) { print outfile2 "\t$samp"; }
	print outfile2 "\n";
	foreach my $kid(sort keys %{ $funstat{$classfun} }) {
                my($cstring,$accum,$funid,$funclass);
                $funclass=$funs{$classfun}{$kid}{path};
                if(!$funclass || ($funclass=~/\|/)) { $funclass="Unclassified"; } 
                if($classfun eq "cog") { $funid="$funclass\t"; }
                $funid.="$kid";
		if($funs{$classfun}{$kid}{fun}) { $funid.=":$funs{$classfun}{$kid}{fun}"; }
		# print outfile2 "$funid";
		$cstring.="$funid";		
		foreach my $samp(sort keys %allsamples) { 
			my $rnum=$funstat{$classfun}{$kid}{$samp}{reads} || "0";
			$accum+=$rnum;
			$cstring.="\t$rnum"; 
			#print outfile2 "\t$funstat{$classfun}{$kid}{$samp}{reads}"; 
			}
		if($accum>=$minraw12) { print outfile2 "$cstring\n"; }
		# print outfile2 "\n";
		}
	close outfile2;
	}


