#!/usr/bin/env perl

use strict;

my $SQLITE_CACHE_SIZE           = 2500000;

my $DOWNLOAD_TEST               = 1;
my $DOWNLOAD_GENERAL            = 1;
my $DOWNLOAD_SILVA              = 1;
my $DOWNLOAD_KEGG               = 1;
my $DOWNLOAD_NR                 = 1;
my $DOWNLOAD_EGGNOG             = 1;
my $DOWNLOAD_PFAM               = 1;
my $DOWNLOAD_CHECKM2            = 1;
my $LCA_RECTAXA                 = 1;
my $LCA_NRINDEX                 = 1;
my $LCA_TAXID_TREE              = 1;
my $LCA_TAXID_DB                = 1;
my $LCA_PARENTS_DB              = 1;
my $REMOVE_NR                   = 1;
my $REMOVE_TAXDUMP              = 1;
my $REMOVE_LCA_TAX_INTERMEDIATE = 1;

if(!$ARGV[0]) { die "Please provide a download directory!   "; }

use Cwd 'abs_path';
my $download_dir = abs_path($ARGV[0]);
my $database_dir = "$download_dir/db";
my $lca_dir = "$database_dir/LCA_tax";

if(!$download_dir) { die "Usage: perl make_databases.pl <download dir>\n"; }

###scriptdir patch v2, Fernando Puente-Sánchez, 18-XI-2019
use File::Basename;

my $dbscriptdir;
if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $dbscriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $dbscriptdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$dbscriptdir/../..");
my $libpath = "$installpath/lib";
require "$libpath/install_utils/download_confirm.pl";
require "$libpath/install_utils/get_host.pl";
###

my $host = get_host();

system("rm $download_dir/test.tar.gz $libpath/classifier.tar.gz $download_dir/db.tar.gz $download_dir/kegg.dmnd.gz > /dev/null 2>&1");

### Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
if($DOWNLOAD_TEST){
	print "\nDownloading and unpacking test data...\n\n";
	download_confirm("test.tar.gz", "test.md5", $host, $download_dir);
}


### Download general db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
if($DOWNLOAD_GENERAL){
	print "Downloading and unpacking general database tarball...\n";
	download_confirm("db.tar.gz", "db.md5", $host, $download_dir);
}


### Download and unpack silva databases
if($DOWNLOAD_SILVA){
	print "Downloading and unpacking SILVA databases (https://www.arb-silva.de/silva-license-information)...\n";
	download_confirm("silva.nr_v132.align.gz", "silva.nr_v132.align.md5",  $host, $database_dir);
	download_confirm("silva.nr_v132.tax.gz", "silva.nr_v132.tax.md5",  $host, $database_dir);
}


### Download and create kegg db.
if($DOWNLOAD_KEGG){
	print "\nDownloading and creating kegg database...\n\n";
	download_confirm("kegg.db.gz", "kegg.db.md5", $host, $database_dir);
	my $command = "$installpath/bin/diamond makedb --in $database_dir/kegg.db -d $database_dir/keggdb -p 8";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, you ran out of disk space, or something is wrong with your DIAMOND binary\n\n"; }
	system("rm $database_dir/kegg.db");
}


### Download and create nr db.
if($DOWNLOAD_NR){
	print "Downloading and creating nr database. This may take a several hours...\n";
	my $command = "perl $libpath/install_utils/make_nr_db_2020.pl $database_dir";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, you ran out of disk space, or something is wrong with your DIAMOND binary\n\n"; }
}


### Download and create eggnog db.
if($DOWNLOAD_EGGNOG){
	print "\nDownloading and creating eggnog database...\n\n";
	my $command = "perl $libpath/install_utils/make_eggnog_db.pl $database_dir";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, you ran out of disk space, or something is wrong with your DIAMOND binary\n\n"; }
}


### Download and create Pfam db.
if($DOWNLOAD_PFAM){
	print "\nDownloading Pfam database...\n\n";
	my $command = "wget -P $database_dir ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; gunzip $database_dir/Pfam-A.hmm.gz";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted or you ran out of disk space\n\n"; }
}


### Download CheckM2 reference db.
if($DOWNLOAD_CHECKM2){
	print "\nDownloading Checkm2 reference...\n\n";
	download_confirm("uniref100.KO.1.dmnd.gz", "uniref100.KO.1.dmnd.md5",  $host, $database_dir);
}


### Create LCA database from nr data.
if($LCA_RECTAXA and $LCA_NRINDEX and $LCA_TAXID_TREE and $LCA_TAXID_DB and $LCA_PARENTS_DB){
	print "\nCreating lca.db...\n";
	system "rm $lca_dir";
	system "mkdir $lca_dir";
}


if($LCA_RECTAXA){
	print "\n  Running rectaxa.pl\n";
	my $command = "perl $libpath/install_utils/rectaxa.pl $lca_dir";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
}
if($REMOVE_TAXDUMP) { system("rm $lca_dir/*dmp $lca_dir/new_taxdump.tar.gz"); }


if($LCA_NRINDEX){
	print "\n  Running nrindex.pl\n";
	my $command = "perl $libpath/install_utils/nrindex.pl $database_dir";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
}
if($REMOVE_NR) { system ("rm -r $database_dir/nr.faa"); }


if($LCA_TAXID_TREE){
	print "\n  Running taxid_tree.pl\n";
	my $command = "perl $libpath/install_utils/taxid_tree.pl $lca_dir";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
	my $command = "sed -i \"s/['\\\"]//g\" $lca_dir/taxid_tree.txt"; # Remove quotes for sqlite3.
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
}


if($LCA_TAXID_DB){
	print "\n  Creating sqlite database taxid.db\n\n";
	if(-e "$lca_dir/taxid.db") { system("rm $lca_dir/taxid.db"); }
	my $command = "sqlite3 $lca_dir/taxid.db < $libpath/install_utils/taxid.sql";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
	my $command = "echo '.import $lca_dir/taxid_tree.txt taxid' | sqlite3 $lca_dir/taxid.db -cmd 'PRAGMA journal_mode=OFF;PRAGMA cache_size=$SQLITE_CACHE_SIZE;PRAGMA synchronous=0;PRAGMA temp_store=2;' -cmd '.separator \"\\t\"'";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
	my $textrows = `wc -l $lca_dir/taxid_tree.txt`;
	my $dbrows = `echo 'SELECT count(*) FROM taxid;' | sqlite3 $lca_dir/taxid.db`;
	if($textrows != $dbrows) { die "Error creating taxid.db, please contact us!\n\n" }
	system("md5sum $lca_dir/taxid.db > $lca_dir/taxid.md5");
}


if($LCA_PARENTS_DB){
	print "\n  Creating sqlite database parents.db\n\n";
	if(-e "$lca_dir/parents.db") { system("rm $lca_dir/parents.db"); }
	my $command = "sqlite3 $lca_dir/parents.db < $libpath/install_utils/parents.sql";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
	my $command = "echo '.import $lca_dir/parents.txt parents' | sqlite3 $lca_dir/parents.db -cmd 'PRAGMA journal_mode=OFF;PRAGMA cache_size=$SQLITE_CACHE_SIZE;PRAGMA synchronous=0;PRAGMA temp_store=2;' -cmd '.separator \"\\t\"'";
	my $ecode = system $command;
	if($ecode!=0) { die "Error running command:     $command\n\n"; }
}

if($REMOVE_LCA_TAX_INTERMEDIATE) { system("rm $lca_dir/nr.taxlist.db $lca_dir/taxid_tree.txt $lca_dir/taxatree.txt"); }


### Report db creation date.
my $timestamp = scalar localtime;
my $time_command = "echo \"Finished database creatin at > $database_dir";
system("echo  \"Finished database creation on $timestamp.\" > $database_dir/DB_BUILD_DATE");


### Finish configuration.
system("perl $installpath/utils/install_utils/configure_nodb.pl $database_dir");
