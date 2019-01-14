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
my $databasedir=abs_path($ARGV[0]);


if(!$databasedir) { die "Usage: perl make_databases.pl <database dir>\n"; }

system("rm $databasedir/test.tar.gz $libpath/classifier.tar.gz $databasedir/db.tar.gz $databasedir/kegg.dmnd.gz");

###Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
print "\nDownloading and unpacking test data...\n\n";
system("wget -U '' -P $databasedir http://wwwuser.cnb.csic.es/~squeezem/test.tar.gz; tar -xvzf $databasedir/test.tar.gz -C $databasedir; rm $databasedir/test.tar.gz");

###Download rdp classifier.
print("Downloading and unpacking RDP classifier...\n");
system("wget -U '' -P $libpath http://wwwuser.cnb.csic.es/~squeezem/classifier.tar.gz; tar -xvzf $libpath/classifier.tar.gz -C $libpath; rm $libpath/classifier.tar.gz");

###Download general db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading and unpacking general database tarball...\n";
system("wget -U '' -P $databasedir http://wwwuser.cnb.csic.es/~squeezem/db.tar.gz; tar -xvzf $databasedir/db.tar.gz -C $databasedir; rm $databasedir/db.tar.gz");


###From now on refer to the database folder.
my $databasedir="$databasedir/db";


###Download and create kegg db.
print "\nDownloading and creating kegg database...\n\n";
system("wget -U '' -P $databasedir http://wwwuser.cnb.csic.es/~squeezem/kegg.db.gz; gunzip $databasedir/kegg.db.gz");
system("$installpath/bin/diamond makedb --in $databasedir/kegg.db -d $databasedir/keggdb -p 8");
system("rm $databasedir/kegg.db");


###Download and create nr db.
print "Downloading and creating nr database. This will take a while (several hours)...\n";
system "perl $dbscriptdir/make_nr_db.pl $databasedir";



###Download and create eggnog db.
print "\nDownloading and creating eggnog database...\n\n";
system "perl $dbscriptdir/make_eggnog_db.pl $databasedir";


###Download and create Pfam db.
print "\nDownloading Pfam database...\n\n";
system "wget -P $databasedir ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz; gunzip $databasedir/Pfam-A.hmm.gz";

my $lca_dir = "$databasedir/LCA_tax";
###Create lca database from nr data.
print "\nCreating lca.db...\n";
system "rm $lca_dir";
system "mkdir $lca_dir";

print "\n  Running rectaxa.pl\n";
system "perl $dbscriptdir/rectaxa.pl $lca_dir";

print "\n  Running nrindex.pl\n";
system "perl $dbscriptdir/nrindex.pl $databasedir";

print "\n  Running taxid_tree.pl\n";
system "perl $dbscriptdir/taxid_tree.pl $lca_dir";
system "sed -i \"s/['\\\"]//g\" $lca_dir/taxid_tree.txt"; # Remove quotes for sqlite3.

print "\n  Creating sqlite databases\n\n";

system "sqlite3 $lca_dir/taxid.db < $dbscriptdir/taxid.sql";
system "echo '.import $lca_dir/taxid_tree.txt taxid' | sqlite3 $lca_dir/taxid.db -cmd '.separator \"\\t\"'";
my $textrows = system "wc -l $lca_dir/taxid_tree.txt";
my $dbrows = system "echo 'SELECT count(*) FROM taxid; | sqlite3 $lca_dir/taxid.db";
if($textrows != $dbrows) { die "Error creating taxid.db, please contact us!" }

system "sqlite3 $lca_dir/parents.db < $dbscriptdir/parents.sql";
system "echo '.import $lca_dir/parents.txt parents' | sqlite3 $lca_dir/LCA_tax/parents.db -cmd '.separator \"\\t\"'";

if($REMOVE_NR) { system ("rm -r $databasedir/nr.faa"); }
if($REMOVE_TAXDUMP) { system("rm $lca_dir/*dmp $lca_dir/new_taxdump.tar.gz"); }
if($REMOVE_LCA_TAX_INTERMEDIATE) { system("rm $lca_dir/nr.taxlist.tsv $lca_dir/LCA_tax/taxid_tree.txt $lca_dir/taxatree.txt"); }


###Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");


my $checkm_manifest = "{\"dataRoot\": \"$databasedir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

open(outfile2,">$installpath/scripts/SqueezeMeta_conf.pl") || die;
open(infile1, "$installpath/scripts/SqueezeMeta_conf_original.pl") || die;
while(<infile1>) {
	if($_=~/^\$databasepath/) { print outfile2 "\$databasepath=\"$databasedir\";\n"; }
	else { print outfile2 $_; }
	}
close infile1;
close outfile2;

print("Done\n");

