#!/usr/bin/perl

use strict;

my $REMOVE_NR=1;
my $REMOVE_TAXDUMP=1;
my $REMOVE_LCA_TAX_INTERMEDIATE=1;

###scriptdir patch, Fernando Puente-Sánchez, 07-V-2018
use File::Basename;
our $dbscriptdir = dirname(__FILE__);
our $installpath = "$dbscriptdir/../..";
our $libpath = "$installpath/lib";
###

use Cwd 'abs_path';
my $download_dir = abs_path($ARGV[0]);
my $database_dir = "$download_dir/db";
my $lca_dir = "$database_dir/LCA_tax";

if(!$download_dir) { die "Usage: perl make_databases.pl <download dir>\n"; }

system("rm $download_dir/test.tar.gz $libpath/classifier.tar.gz $download_dir/db.tar.gz $download_dir/kegg.dmnd.gz");


###Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
print "\nDownloading and unpacking test data...\n\n";
system("wget -U '' -P $download_dir http://wwwuser.cnb.csic.es/~squeezem/test.tar.gz; tar -xvzf $download_dir/test.tar.gz -C $download_dir; rm $download_dir/test.tar.gz");


###Download rdp classifier.
print("Downloading and unpacking RDP classifier...\n");
system("wget -U '' -P $libpath http://wwwuser.cnb.csic.es/~squeezem/classifier.tar.gz; tar -xvzf $libpath/classifier.tar.gz -C $libpath; rm $libpath/classifier.tar.gz");


###Download general db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading and unpacking general database tarball...\n";
system("wget -U '' -P $download_dir http://wwwuser.cnb.csic.es/~squeezem/db.tar.gz; tar -xvzf $download_dir/db.tar.gz -C $download_dir; rm $download_dir/db.tar.gz");


###Download and create kegg db.
print "\nDownloading and creating kegg database...\n\n";
system("wget -U '' -P $database_dir http://wwwuser.cnb.csic.es/~squeezem/kegg.db.gz; gunzip $database_dir/kegg.db.gz");
system("$installpath/bin/diamond makedb --in $database_dir/kegg.db -d $database_dir/keggdb -p 8");
system("rm $database_dir/kegg.db");


###Download and create nr db.
print "Downloading and creating nr database. This will take a while (several hours)...\n";
system "perl $dbscriptdir/make_nr_db.pl $database_dir";


###Download and create eggnog db.
print "\nDownloading and creating eggnog database...\n\n";
system "perl $dbscriptdir/make_eggnog_db.pl $database_dir";


###Download and create Pfam db.
print "\nDownloading Pfam database...\n\n";
system "wget -P $database_dir ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; gunzip $database_dir/Pfam-A.hmm.gz";


###Create LCA database from nr data.
print "\nCreating lca.db...\n";
system "rm $lca_dir";
system "mkdir $lca_dir";

print "\n  Running rectaxa.pl\n";
system "perl $dbscriptdir/rectaxa.pl $lca_dir";

print "\n  Running nrindex.pl\n";
system "perl $dbscriptdir/nrindex.pl $database_dir";

print "\n  Running taxid_tree.pl\n";
system "perl $dbscriptdir/taxid_tree.pl $lca_dir";
system "sed -i \"s/['\\\"]//g\" $lca_dir/taxid_tree.txt"; # Remove quotes for sqlite3.

print "\n  Creating sqlite databases\n\n";

system "sqlite3 $lca_dir/taxid.db < $dbscriptdir/taxid.sql";
system "echo '.import $lca_dir/taxid_tree.txt taxid' | sqlite3 $lca_dir/taxid.db -cmd '.separator \"\\t\"'";
my $textrows = `wc -l $lca_dir/taxid_tree.txt`;
my $dbrows = `echo 'SELECT count(*) FROM taxid;' | sqlite3 $lca_dir/taxid.db`;
if($textrows != $dbrows) { die "Error creating taxid.db, please contact us!" }

system "sqlite3 $lca_dir/parents.db < $dbscriptdir/parents.sql";
system "echo '.import $lca_dir/parents.txt parents' | sqlite3 $lca_dir/parents.db -cmd '.separator \"\\t\"'";

if($REMOVE_NR) { system ("rm -r $database_dir/nr.faa"); }
if($REMOVE_TAXDUMP) { system("rm $lca_dir/*dmp $lca_dir/new_taxdump.tar.gz"); }
if($REMOVE_LCA_TAX_INTERMEDIATE) { system("rm $lca_dir/nr.taxlist.tsv $lca_dir/taxid_tree.txt $lca_dir/taxatree.txt"); }


###Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");


my $checkm_manifest = "{\"dataRoot\": \"$database_dir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

open(outfile2,">$installpath/scripts/SqueezeMeta_conf.pl") || die;
open(infile1, "$installpath/scripts/SqueezeMeta_conf_original.pl") || die;
while(<infile1>) {
	if($_=~/^\$databasepath/) { print outfile2 "\$databasepath=\"$database_dir\";\n"; }
	else { print outfile2 $_; }
	}
close infile1;
close outfile2;

print("Done\n");

