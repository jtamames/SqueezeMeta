#!/usr/bin/env perl

use strict;

if(!$ARGV[0]) { die "Please provide a download directory!   "; }

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
require "$libpath/install_utils/download_confirm.pl";
require "$libpath/install_utils/get_host.pl";
###

my $host = get_host();

system("rm $download_dir/test.tar.gz $libpath/classifier.tar.gz $download_dir/SqueezeMetaDB.tar.gz > /dev/null 2>&1");

### Download test data (-U '' so that we give the server an user agent string, it complains otherwise).
print "\nDownloading and unpacking test data...\n\n";
download_confirm("test.tar.gz", "test.md5", $host, $download_dir);

### Download db tarball. (-U '' so that we give the server an user agent string, it complains otherwise)
print "Downloading and unpacking database tarball...\n";
download_confirm("SqueezeMetaDB.tar.gz", "SqueezeMetaDB.md5", $host, $download_dir);


### Finish configuration.
system("perl $installpath/utils/install_utils/configure_nodb.pl $database_dir");
