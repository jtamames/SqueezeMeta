#!/usr/bin/perl

use strict;

###scriptdir patch, Fernando Puente-Sánchez, 07-V-2018
use File::Basename;
our $dbscriptdir = dirname(__FILE__);
our $installpath = "$dbscriptdir/../..";
###

use Cwd 'abs_path';
my $databasedir=abs_path($ARGV[0]);


if(!$databasedir) { die "Usage: perl make_databases.pl <database dir>\n"; }

system("rm $databasedir/test.tar.gz $databasedir/db.tar.gz $databasedir/kegg.dmnd.gz");
###Download test data.
print "\nDownloading and unpacking test data...\n\n";
system("wget -U '' -P $databasedir http://wwwuser.cnb.csic.es/~squeezem/test.tar.gz; tar -xvzf $databasedir/test.tar.gz -C $databasedir; rm $databasedir/test.tar.gz");


###Download general db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading an unpacking general database tarball...\n";
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


###Create lca database from nr data.
print "\nCreating lca.db...\n";
system "mkdir $databasedir/LCA_tax";
print "\n  Running rectaxa.pl\n";
system "perl $dbscriptdir/rectaxa.pl $databasedir";
print "\n  Running nrindex.pl\n";
system "perl $dbscriptdir/nrindex.pl $databasedir";
print "\n  Running taxid_tree.pl\n";
system "perl $dbscriptdir/taxid_tree.pl $databasedir";

system ("rm $databasedir/nr.faa");

print "\n  Creating sqlite databases\n\n";

system "sqlite3 $databasedir/LCA_tax/taxid.db < $dbscriptdir/taxid.sql";
system "echo '.import $databasedir/LCA_tax/taxid_tree.txt taxid' | sqlite3 $databasedir/LCA_tax/taxid.db -cmd '.separator \"\\t\"'";

system "sqlite3 $databasedir/LCA_tax/parents.db < $dbscriptdir/parents.sql";
system "echo '.import $databasedir/LCA_tax/parents.txt parents' | sqlite3 $databasedir/LCA_tax/parents.db -cmd '.separator \"\\t\"'";

system("rm $databasedir/LCA_tax/nr.taxlist.tsv $databasedir/LCA_tax/taxid_tree.txt $databasedir/LCA_tax/taxatree.txt");


###Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");

#
my $checkm_manifest = "{\"dataRoot\": \"$databasedir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

open(outfile2,">$installpath/scripts/squeezeM_conf.pl") || die;
open(infile1, "$installpath/scripts/squeezeM_conf_original.pl") || die;
while(<infile1>) {
	if($_=~/^\$databasepath/) { print outfile2 "\$databasepath=\"$databasedir\";\n"; }
	else { print outfile2 $_; }
	}
close infile1;
close outfile2;

print("Done\n");

