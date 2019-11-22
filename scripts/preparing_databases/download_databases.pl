#!/usr/bin/env perl

use strict;

my $REMOVE_NR=1;
my $REMOVE_TAXDUMP=1;
my $REMOVE_LCA_TAX_INTERMEDIATE=1;

use Cwd 'abs_path';
my $download_dir = abs_path($ARGV[0]);
my $database_dir = "$download_dir/db";

if(!$download_dir) { die "Usage: perl download_databases.pl <download dir>\n"; }

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
###

system("rm $download_dir/test.tar.gz $libpath/classifier.tar.gz $download_dir/SqueezeMetaDB.tar.gz");


### Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
print "\nDownloading and unpacking test data...\n\n";
system("wget -U '' -P $download_dir http://wwwuser.cnb.csic.es/~squeezem/test.tar.gz; tar -xvzf $download_dir/test.tar.gz -C $download_dir; rm $download_dir/test.tar.gz");


### Download rdp classifier.
print("Downloading and unpacking RDP classifier...\n");
system("wget -U '' -P $libpath http://wwwuser.cnb.csic.es/~squeezem/classifier.tar.gz; tar -xvzf $libpath/classifier.tar.gz -C $libpath; rm $libpath/classifier.tar.gz");


### Download db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading and unpacking database tarball...\n";
system("wget -U '' -P $download_dir http://silvani.cnb.csic.es/SqueezeMeta/SqueezeMetaDB.tar.gz; tar -xvzf $download_dir/SqueezeMetaDB.tar.gz -C $download_dir; rm $download_dir/SqueezeMetaDB.tar.gz");


### Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");

my $checkm_manifest = "{\"dataRoot\": \"$database_dir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

open(outfile2,">$installpath/scripts/SqueezeMeta_conf.pl") || die;
open(infile1, "$installpath/scripts/SqueezeMeta_conf_original.pl") || die;
while(<infile1>) {
	if($_=~/^\$databasepath/) { print outfile2 "\$databasepath = \"$database_dir\";\n"; }
	else { print outfile2 $_; }
	}
close infile1;
close outfile2;

print("Done\n");


