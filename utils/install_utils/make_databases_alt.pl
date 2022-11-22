#!/usr/bin/env perl

use strict;

my $REMOVE_NR=1;
my $REMOVE_TAXDUMP=1;
my $REMOVE_LCA_TAX_INTERMEDIATE=1;

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
###


system("rm $download_dir/test.tar.gz $libpath/classifier.tar.gz $download_dir/db.tar.gz $download_dir/kegg.dmnd.gz > /dev/null 2>&1");


### Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
print "\nDownloading and unpacking test data...\n\n";
system("wget -P $download_dir -O $download_dir/test.tar.gz https://saco.csic.es/index.php/s/qo264F9BAbjmp3e/download");
system("tar -xvzf $download_dir/test.tar.gz -C $download_dir; rm $download_dir/test.tar.gz");


### Download general db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading and unpacking general database tarball...\n";
system("wget -P $download_dir -O $download_dir/db.tar.gz https://saco.csic.es/index.php/s/isXq28Ms5Zk9NjR/download");
system("tar -xvzf $download_dir/db.tar.gz -C $download_dir; rm $download_dir/db.tar.gz");


### Download and unpack silva databases
print "Downloading and unpacking SILVA databases (https://www.arb-silva.de/silva-license-information)...\n";
system("wget -P $database_dir -O $database_dir/silva.nr_v132.align.gz https://saco.csic.es/index.php/s/9M6QpbkqfscATf2/download");
system("gunzip $database_dir/silva.nr_v132.align.gz");
system("wget -P $database_dir -O $database_dir/silva.nr_v132.tax.gz https://saco.csic.es/index.php/s/HNnF5knj4YgM577/download");
system("gunzip $database_dir/silva.nr_v132.tax.gz");


### Download and create kegg db.
print "\nDownloading and creating kegg database...\n\n";
system("wget -P $database_dir -O $database_dir/kegg.db.gz https://saco.csic.es/index.php/s/AkiPW2k2wtwRNcY/download");
my $command = "gunzip $database_dir/kegg.db.gz";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, or you ran out of disk space"; }
my $command = "$installpath/bin/diamond makedb --in $database_dir/kegg.db -d $database_dir/keggdb -p 8";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, you ran out of disk space, or something is wrong with your DIAMOND binary"; }
system("rm $database_dir/kegg.db");


### Download and create nr db.
print "Downloading and creating nr database. This may take a several hours...\n";
system "perl $libpath/install_utils/make_nr_db_2020.pl $database_dir";


### Download and create eggnog db.
print "\nDownloading and creating eggnog database...\n\n";
system "perl $libpath/install_utils/make_eggnog_db.pl $database_dir";


### Download and create Pfam db.
print "\nDownloading Pfam database...\n\n";
system "wget -P $database_dir ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; gunzip $database_dir/Pfam-A.hmm.gz";


### Create LCA database from nr data.
print "\nCreating lca.db...\n";
system "rm $lca_dir";
system "mkdir $lca_dir";

print "\n  Running rectaxa.pl\n";
system "perl $libpath/install_utils/rectaxa.pl $lca_dir";

print "\n  Running nrindex.pl\n";
system "perl $libpath/install_utils/nrindex.pl $database_dir";

print "\n  Running taxid_tree.pl\n";
system "perl $libpath/install_utils/taxid_tree.pl $lca_dir";
system "sed -i \"s/['\\\"]//g\" $lca_dir/taxid_tree.txt"; # Remove quotes for sqlite3.

print "\n  Creating sqlite databases\n\n";

my $command = "sqlite3 $lca_dir/taxid.db < $libpath/install_utils/taxid.sql";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n"; }

my $command = "echo '.import $lca_dir/taxid_tree.txt taxid' | sqlite3 $lca_dir/taxid.db -cmd '.separator \"\\t\"'";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n"; }

my $textrows = `wc -l $lca_dir/taxid_tree.txt`;
my $dbrows = `echo 'SELECT count(*) FROM taxid;' | sqlite3 $lca_dir/taxid.db`;
if($textrows != $dbrows) { die "Error creating taxid.db, please contact us!" }
system("md5sum $lca_dir/taxid.db > $lca_dir/taxid.md5");

my $command = "sqlite3 $lca_dir/parents.db < $libpath/install_utils/parents.sql";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n"; }

my $command = "echo '.import $lca_dir/parents.txt parents' | sqlite3 $lca_dir/parents.db -cmd '.separator \"\\t\"'";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n"; }

if($REMOVE_NR) { system ("rm -r $database_dir/nr.faa"); }
if($REMOVE_TAXDUMP) { system("rm $lca_dir/*dmp $lca_dir/new_taxdump.tar.gz"); }
if($REMOVE_LCA_TAX_INTERMEDIATE) { system("rm $lca_dir/nr.taxlist.db $lca_dir/taxid_tree.txt $lca_dir/taxatree.txt"); }


### Report db creation date.
my $timestamp = scalar localtime;
my $time_command = "echo \"Finished database creatin at > $database_dir";
system("echo  \"Finished database creation on $timestamp.\" > $database_dir/DB_BUILD_DATE");


### Finish configuration.
system("perl $installpath/utils/install_utils/configure_nodb_alt.pl $database_dir");
