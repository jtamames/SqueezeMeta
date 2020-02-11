#!/usr/bin/env perl

use strict;

if(!$ARGV[0]) { die "Please provide a database directory!   "; }
use Cwd 'abs_path';
my $databasedir=abs_path($ARGV[0]);

if(!$databasedir) { die "Usage: perl install_nodb.pl <database dir>\n"; }
print("Make sure that $databasedir contains all the database files (nr.dmnd, etc...)\n\n");


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


system("rm $libpath/classifier.tar.gz");

###Download rdp classifier.
print("Downloading and unpacking RDP classifier...\n");
system("wget -U '' -P $libpath http://wwwuser.cnb.csic.es/~squeezem/classifier.tar.gz; tar -xvzf $libpath/classifier.tar.gz -C $libpath; rm $libpath/classifier.tar.gz");
system("cd $installpath/bin/; ln -s $libpath/classifier/classifier.jar ."); # Add symlink

###Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");


my $checkm_manifest = "{\"dataRoot\": \"$databasedir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

open(outfile2,">$installpath/scripts/SqueezeMeta_conf.pl") || die;
open(infile1, "$installpath/scripts/SqueezeMeta_conf_original.pl") || die;
while(<infile1>) {
	if($_=~/^\$databasepath/) { print outfile2 "\$databasepath = \"$databasedir\";\n"; }
	else { print outfile2 $_; }
	}
close infile1;
close outfile2;

print("Done\n");

