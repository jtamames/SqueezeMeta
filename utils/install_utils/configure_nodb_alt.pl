#!/usr/bin/env perl

use strict;

if(!$ARGV[0]) { die "Please provide a database directory!   "; }
use Cwd 'abs_path';
my $databasedir=abs_path($ARGV[0]);

if(!$databasedir) { die "Usage: perl configure_nodb_alt.pl <database dir>\n"; }
print("\nMake sure that $databasedir contains all the database files (nr.dmnd, etc...)\n\n");


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

###Check that the databases are present in the requested dir.
unless(-e "$databasedir/nr.dmnd"){
        die "The file $databasedir/nr.dmnd does not exist. Please check that you provided the right database directory. If you are running configure_nodb.pl directly, this probably means that you did not provide the right database directory. If this happened while running download_databases.pl or make_databases.pl, this probably means that your download got interrupted.\n\n";
}

###Download rdp classifier.
system("rm $libpath/classifier.tar.gz > /dev/null 2>&1");
print("Downloading and unpacking RDP classifier...\n");
system("wget -P $libpath -O $libpath/classifier.tar.gz https://saco.csic.es/index.php/s/D46ieFfdFZirXK5/download");
my $command = "tar -xvzf $libpath/classifier.tar.gz -C $libpath; rm $libpath/classifier.tar.gz";
my $ecode = system $command;
if($ecode!=0) { die "Error running command:     $command\n\nThis probably means that your download got interrupted, or that you ran out of disk space"; }
system("cd $installpath/bin/; ln -s $libpath/classifier/classifier.jar . > /dev/null 2>&1"); # Add symlink

###Update configuration files to reflect new db path.
print("\nUpdating configuration...\n");


my $checkm_manifest = "{\"dataRoot\": \"$databasedir\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"localManifestName\": \".dmanifest\", \"remoteManifestName\": \".dmanifest\"}\n";

open(outfile1, ">$installpath/lib/checkm/DATA_CONFIG") || die;
print outfile1 $checkm_manifest;
close outfile1;

my $checkm2_dmnd_path = "{\"Type\": \"DIAMONDDB\", \"DBPATH\": \"$databasedir/uniref100.KO.1.dmnd\"}";                                                                                                                                                                                                                                                                                                                            open(outfile2, ">$installpath/lib/checkm2/version/diamond_path.json") || die;                                                                                                                                    print outfile2 $checkm2_dmnd_path;                                                                                                                                                                               close outfile2;                                                                                                                                                                                                                                                                                                                                                                                                                   open(outfile3,">$installpath/scripts/SqueezeMeta_conf.pl") || die;                                                                                                                                               open(infile1, "$installpath/scripts/SqueezeMeta_conf_original.pl") || die;                                                                                                                                       while(<infile1>) {                                                                                                                                                                                                       if($_=~/^\$databasepath/) { print outfile3 "\$databasepath = \"$databasedir\";\n"; }                                                                                                                             else { print outfile3 $_; }                                                                                                                                                                                      }                                                                                                                                                                                                        close infile1;                                                                                                                                                                                                   close outfile3;

print("Done\n");

###Test environment.
print("\nTesting environment...\n");
system("perl $installpath/utils/install_utils/test_install.pl");

