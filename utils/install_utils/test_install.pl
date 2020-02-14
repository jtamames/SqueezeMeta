#!/usr/bin/env perl

use strict;
use warnings;

my $scriptdir;
my $ecode;



# Where are we running from?
use File::Basename;
use Cwd 'abs_path';

if(-l __FILE__)
        {
        my $symlinkpath = dirname(__FILE__);
        my $symlinkdest = readlink(__FILE__);
        $scriptdir = dirname(abs_path("$symlinkpath/$symlinkdest"));
        }
else
        {
        $scriptdir = abs_path(dirname(__FILE__));
        }
my $installpath = abs_path("$scriptdir/../..");
###

if(!-e "$installpath/scripts/SqueezeMeta_conf_original.pl") { die ("\nCRITICAL ERROR: Can not find SqueezeMeta_conf_original.pl. We actually thought this was impossible. If the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n"); }

our $warnings;

print("\n");
print("Checking for gcc\n");
check_command("gcc --help", "ERROR: The GCC compiler can not be found in this environment!!");


print("\n");
print("Checking that ruby is installed\n");
check_command("ruby -h", "ERROR: The ruby interpreter can not be found in this environment!!");


print("\n");
print("Checking that java is installed\n");
check_command("java -h", "ERROR: The java VM can not be found in this environment!!");


print("\n");
print("Checking that all the required perl libraries are available in this environment\n");
check_perl_library("Term::ANSIColor");
check_perl_library("DBI");
check_perl_library("DBD::SQLite::Constants");
check_perl_library("Time::Seconds");
check_perl_library("Tie::IxHash");
check_perl_library("Linux::MemInfo");
check_perl_library("Getopt::Long");  # If this script is running then this is 100% present but meh...
check_perl_library("File::Basename");
check_perl_library("DBD::SQLite");
check_perl_library("Data::Dumper");
check_perl_library("Cwd"); # If this script is running then this is 100% present but meh...
check_perl_library("XML::LibXML");


print("\n");
print("Checking that all the required python libraries are available in this environment\n");
$ecode = check_command("python3 -h", "The python3 interpreter can not be found in this environment!!");
if(!$ecode) {
	check_python_library("numpy");
	check_python_library("scipy");
	check_python_library("matplotlib");
	check_python_library("dendropy");
	check_python_library("pysam");
	check_python_library("pandas");
	check_python_library("cython");
	#check_python_library("madeToFail");
}


print("\n");
print("Checking that all the required R libraries are available in this environment\n");
$ecode = check_command("R -h", "The R interpreter can not be found in this environment!!");
if(!$ecode) {
	check_R_library("doMC");
	check_R_library("ggplot2");
	check_R_library("data.table");
	check_R_library("reshape2");
	check_R_library("pathview");
	check_R_library("DASTool");
	check_R_library("SQMtools");
	#check_R_library("madeToFail");
}

print("\n");
print("Checking that SqueezeMeta is properly configured\n");
if(!-e "$installpath/scripts/SqueezeMeta_conf.pl") {
	print("\n");
	warn("\tSqueezeMeta doesn't know where the databases are located!\n");
	print("\n");
	print("\tIf you didn't download them yet please run:\n");
	print("\n");
	print("\t* perl $installpath/utils/install_utils/make_databases.pl /download/path/\n");
	print("\t\t(To download the latest source data and compile the databases in your server)\n");
	print("\n");
	print("\t* perl $installpath/utils/install_utils/download_databases.pl /download/path/\n");
	print("\t\t(To download a pre-compiled version of the database, which is much quicker)\n");
	print("\n");
	print("\tIf the databases are already present in your server, you can configure this installation of SqueezeMeta with:\n");
	print("\n");
	print("\t* perl $installpath/utils/install_utils/configure_nodb.pl /path/to/db\n");
	print("\t\t(Where /path/to/db is the path to the \"db\" folder generated when downloading the databases)\n");
	print("\n");
	$warnings .= "\t- Databases were not installed, or are not configured in this installation of SqueezeMeta\n";
} else {
	do "$installpath/scripts/SqueezeMeta_conf.pl";
	our $databasepath;
	if(!-e "$databasepath/nr.dmnd") { warn("\tSqueezeMeta_conf.pl says that databases are located in $databasepath but we can't find nr.db there\n"); }
	else {
		print("\tDatabases OK\n");
		open( infile_, "$installpath/lib/checkm/DATA_CONFIG" ) || die ("\nCRITICAL ERROR: Can not find the checkm DATA_CONFIG file in $installpath/lib/checkm. This indicates a broken installation. If the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n");
		my $manifest = <infile_>;
		my @parsed_manifest = split(/\: |\, /, $manifest);
		my $checkm_databasepath = $parsed_manifest[1];
		$checkm_databasepath =~ s/\"//g;
		if($checkm_databasepath ne $databasepath) { die("CRITICAL ERROR: the database path in the checkM manifest does not match with the database path in the SqueezeMeta_conf.pl file. This indicates a broken installation. If the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n");  }
		else  { print("\tCheckM manifest OK\n"); }
	}
}

if($warnings) {
	print("\n");
	print("------------------------------------------------------------------------------\n");
	print("\n");
	print("WARNING: Some SqueezeMeta dependencies could not be found in your environment!\n");
	print($warnings);
} else {
	print("\n");
	print("All checks successful\n");
}


print("\n");


sub check_command {
	my $command = $_[0];
	my $msg = $_[1];
	my $ecode = system("$command > /dev/null 2>&1");
        if(!$ecode) { print("\t$command OK\n"    ); }
	else {
		warn("\t$command NOT OK\n"); 
		if($msg) { print("\t\t$msg\n"); }
		$warnings .= "\t- $msg\n";
	}
	return $ecode;
}


sub check_perl_library {
        my $library = $_[0];
        my $msg = $_[1];
	if(!$msg) { $msg = "Missing perl library \"$library\""; }
	my $command = "perl -e 'use $library'";
	return check_command($command, $msg);
}


sub check_python_library {
	my $library = $_[0];
	my $msg = $_[1];
	if(!$msg) { $msg = "Missing python library \"$library\""; }
	my $command = "python3 -c 'import $library'";
        return check_command($command, $msg);
}


sub check_R_library {
	my $library = $_[0];
	my $msg = $_[1];
	if(!$msg) { $msg = "Missing R library \"$library\""; }
        my $command = "R -e 'library($library)'";
        return check_command($command, $msg);
}
