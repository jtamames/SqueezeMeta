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
our $installpath = abs_path("$scriptdir/../..");

###

if(!-e "$installpath/scripts/SqueezeMeta_conf_original.pl") { die ("\nCRITICAL ERROR: Can not find SqueezeMeta_conf_original.pl. We actually thought this was impossible. If the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n"); }

do "$installpath/scripts/SqueezeMeta_conf_original.pl";

our($spades_soft, $metabat_soft, $jgi_summ_soft, $samtools_soft, $bwa_soft, $minimap2_soft, $diamond_soft, $hmmer_soft, $cdhit_soft, $kmerdb_soft, $aragorn_soft, $mothur_soft);

our $warnings;

#print("\n");
#print("Checking for gcc\n");
#check_command("gcc --help", "ERROR: The GCC compiler can not be found in this environment!!");

print("\n");
print("Checking the OS\n");
my $os = "$^O";
if($os eq "linux") { print("\t$os OK\n"); } else { die("The \"$os\" OS was detected, but SqueezeMeta only runs on Linux!!\n\n"); }


print("\n");
print("Checking that tree is installed\n");
check_command("tree --help", "ERROR: The tree program can not be found in this environment!!");


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
my $dbd_sqlite_error = check_perl_library("DBD::SQLite::Constants");
check_perl_library("Time::Seconds");
check_perl_library("Tie::IxHash");
check_perl_library("Linux::MemInfo");
check_perl_library("Getopt::Long");  # If this script is running then this is 100% present but meh...
check_perl_library("File::Basename");
check_perl_library("DBD::SQLite");
check_perl_library("Data::Dumper");
check_perl_library("Cwd"); # If this script is running then this is 100% present but meh...
check_perl_library("XML::LibXML");
check_perl_library("XML::Parser");
check_perl_library("Term::ANSIColor");


print("\n");
print("Checking that all the required python libraries are available in this environment\n");
$ecode = check_command("python3 -h", "The python3 interpreter can not be found in this environment!!");
if(!$ecode) {
	check_python_library("numpy");
	check_python_library("scipy");
	check_python_library("matplotlib");
	check_python_library("dendropy");
	check_python_library("pysam");
	check_python_library("Bio.Seq");
	check_python_library("pandas");
	check_python_library("sklearn");
	#check_python_library("nose");
	check_python_library("cython");
        check_python_library("future");
	check_python_library("tensorflow");
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
print("Checking binaries\n"); ### Some can fail if the user doesn't have write permissions in pwd!
my $tmpdir = `mktemp -d`;
chomp($tmpdir);
check_command("cd $tmpdir; $spades_soft --test");
check_command("$metabat_soft -h");
check_command("$jgi_summ_soft -h");
check_command("$samtools_soft --version");
check_command("$bwa_soft mem $scriptdir/test_data/ctgs.fasta $scriptdir/test_data/seqs.fq"); 
check_command("$minimap2_soft --version");
check_command("$diamond_soft version");
check_command("$hmmer_soft -h");
check_command("$cdhit_soft -i $scriptdir/test_data/ctgs.fasta -o $tmpdir/foo.fasta");
check_command("$kmerdb_soft -h");
check_command("$aragorn_soft -h");
check_command("$mothur_soft -h");
system("rm -r $tmpdir > /dev/null 2>&1");


print("\n");
print("Checking that SqueezeMeta is properly configured...");
if(!-e "$installpath/scripts/SqueezeMeta_conf.pl") {
	print("\n\n");
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
	print(" checking database in $databasepath\n");
	my $ecode = system("$installpath/bin/diamond dbinfo --db $databasepath/nr.dmnd >/dev/null 2>&1");
	if($ecode) {
		my $msg = "\tCRITICAL ERROR: SqueezeMeta_conf.pl says that databases are located in $databasepath but we can't find nr.db there, or it is corrupted\n\n";
		$warnings .= $msg;
		warn($msg);
	}
	else {
		print("\tnr.db OK\n");
		
		open( infile_, "$installpath/lib/checkm/DATA_CONFIG" ) || die ("\nCRITICAL ERROR: Can not find the checkm DATA_CONFIG file in $installpath/lib/checkm. This indicates a broken installation.\n\t\tIf the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n");
		my $manifest = <infile_>;
		my @parsed_manifest = split(/\: |\, /, $manifest);
		my $checkm_databasepath = $parsed_manifest[1];
		$checkm_databasepath =~ s/\"//g;
		if($checkm_databasepath ne $databasepath) { die("CRITICAL ERROR: the database path in the checkM manifest does not match with the database path in the SqueezeMeta_conf.pl file. This indicates a broken installation. If the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n");  }
		else  { print("\tCheckM manifest OK\n"); }
		
		my $checkm2_diamond_path_file = "$installpath/lib/checkm2/version/diamond_path.json";
		open( infile_, $checkm2_diamond_path_file) || die ("\nCRITICAL ERROR: Can not find the checkm2 diamond_path.json file in $installpath/lib/checkm2. This indicates a broken installation.\n\t\tIf the error persists after reinstalling from scratch please open an issue at http://github.com/jtamames/SqueezeMeta\n\n");
		$manifest = <infile_>;
		@parsed_manifest = split(/\: |\, /, $manifest);
		my $checkm2_diamond_path = $parsed_manifest[3];
		$checkm2_diamond_path =~ s/\"|}//g;
		my $ecode = system("$installpath/bin/diamond dbinfo --db $checkm2_diamond_path >/dev/null 2>&1");
		if($ecode) { die("CRITICAL ERROR: the $checkm2_diamond_path_file points to the checkm2 database being located at $checkm2_diamond_path, but we can't find it there, or it is corrupted\n\n");  }
		else  { print("\tCheckM database OK\n"); }
		
		if(!$dbd_sqlite_error) {
			my $ecode = system("perl $installpath/lib/install_utils/test_sqlite_db.pl $databasepath >/dev/null 2>&1");
			if($ecode) {
				my $msg = "\tCRITICAL ERROR: The LCA_tax/taxid.db database is not present in $databasepath, it is malformed, or there is other problem with your SQLite configuration\n";
				$warnings .= $msg;
		       		warn($msg);
			}
			else {
				print("\tLCA_tax DB OK\n");
		                my $db_build_date_file = "$databasepath/DB_BUILD_DATE";
				if(!-e $db_build_date_file) {
					my $msg = "\tCRITICAL ERROR: SqueezeMeta_conf.pl says that databases are located in $databasepath but we can't find the DB_BUILD_DATE file there.\n\t\tThis can happen if make_databases.pl failed to complete successfully\n\n";
					$warnings .= $msg;
					warn($msg);
				}
				else {
	                		open my $fh, "<", $db_build_date_file;
			                chomp(my $db_build_date = do { local $/; <$fh> }); # read the file in one line
					$db_build_date =~ s/Finished database creation on //;
					$db_build_date =~ s/\.$//; 
		                	print("\tDatabase was built on $db_build_date\n");
				}
			}
		}
	}
}

if($warnings) {
	print("\n");
	print("------------------------------------------------------------------------------\n");
	print("\n");
	print("WARNING: Some SqueezeMeta dependencies could not be found in your environment or failed to execute!\n");
	print($warnings);
	die("\n");
} else {
	print("\n");
	print("All checks successful\n");
}


print("\n");


sub check_command {
	my $command = $_[0];
	my $msg = $_[1];
	my $out;
	my @args;
	my $c2;
	if(!$msg) {
		@args = split ';', $command;
		$c2 = $args[-1];
		@args = split ' ', $c2;
		@args = grep !/PATH/, @args; # remove leading env variables before the actual command
		$out = basename($args[0]);
		$msg = "ERROR: Error running $command";
	}
	else { $out = $command; }
	my $ecode = system("$command > /dev/null 2>&1");
        if(!$ecode) { print("\t$out OK\n"    ); }
	else {
		warn("\t$out NOT OK\n"); 
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
