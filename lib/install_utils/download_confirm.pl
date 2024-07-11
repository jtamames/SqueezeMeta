use strict;
use warnings;
use Digest::MD5 qw(md5_hex);

my %host2urls;
$host2urls{"https://saco.csic.es"}{"SQMhere"}                 = 's/ZFLcptDpWGxAT98/download?path=%2F&files=SQMhere';
$host2urls{"https://saco.csic.es"}{"classifier.md5"}          = 's/ZFLcptDpWGxAT98/download?path=%2F&files=classifier.md5';
$host2urls{"https://saco.csic.es"}{"classifier.tar.gz"}       = 's/ZFLcptDpWGxAT98/download?path=%2F&files=classifier.tar.gz';
$host2urls{"https://saco.csic.es"}{"db.md5"}                  = 's/ZFLcptDpWGxAT98/download?path=%2F&files=db.md5';
$host2urls{"https://saco.csic.es"}{"db.tar.gz"}               = 's/ZFLcptDpWGxAT98/download?path=%2F&files=db.tar.gz';
$host2urls{"https://saco.csic.es"}{"kegg.db.md5"}             = 's/ZFLcptDpWGxAT98/download?path=%2F&files=kegg.db.md5';
$host2urls{"https://saco.csic.es"}{"kegg.db.tar.gz"}          = 's/ZFLcptDpWGxAT98/download?path=%2F&files=kegg.db.tar.gz';
$host2urls{"https://saco.csic.es"}{"silva.nr_v132.align.md5"} = 's/ZFLcptDpWGxAT98/download?path=%2F&files=silva.nr_v132.align.md5';
$host2urls{"https://saco.csic.es"}{"silva.nr_v132.align.gz"}  = 's/ZFLcptDpWGxAT98/download?path=%2F&files=silva.nr_v132.align.gz';
$host2urls{"https://saco.csic.es"}{"SqueezeMetaDB.md5"}       = 's/ZFLcptDpWGxAT98/download?path=%2F&files=SqueezeMetaDB.md5';
$host2urls{"https://saco.csic.es"}{"SqueezeMetaDB.tar.gz"}    = 's/ZFLcptDpWGxAT98/download?path=%2F&files=SqueezeMetaDB.tar.gz';
$host2urls{"https://saco.csic.es"}{"test.md5"}                = 's/ZFLcptDpWGxAT98/download?path=%2F&files=test.md5';
$host2urls{"https://saco.csic.es"}{"test.tar.gz"}             = 's/ZFLcptDpWGxAT98/download?path=%2F&files=test.tar.gz';
$host2urls{"https://saco.csic.es"}{"uniref100.KO.1.dmnd.md5"} = 's/ZFLcptDpWGxAT98/download?path=%2F&files=uniref100.KO.1.dmnd.md5';
$host2urls{"https://saco.csic.es"}{"uniref100.KO.1.dmnd.gz"}  = 's/ZFLcptDpWGxAT98/download?path=%2F&files=uniref100.KO.1.dmnd.gz';


sub download_confirm {
	my ($file_name, $md5_name, $host, $download_dir) = @_;

	system("rm $download_dir/$file_name > /dev/null 2>&1");
        system("rm $download_dir/$md5_name  > /dev/null 2>&1");
	
	run_command(build_wget_command($file_name, $host, $download_dir));
	run_command(build_wget_command($md5_name,  $host, $download_dir));

	open( my $QUERY, "$download_dir/$file_name" ) || die "Can't open $download_dir/$file_name\n";
	binmode($QUERY);
	my $md5 = Digest::MD5->new->addfile($QUERY)->hexdigest;
	close $QUERY;

	open(my $CONFIRM, "$download_dir/$md5_name") || die "Can't open $download_dir/$md5_name\n";
	$_ = <$CONFIRM>;
	chomp;
	close $CONFIRM;

	if($md5 eq $_) {
		if($file_name =~ /\.tar.gz\z/) {
			run_command("tar -xvzf $download_dir/$file_name -C $download_dir; rm $download_dir/$file_name $download_dir/$md5_name");
		} elsif ($file_name =~ /\.gz\z/) {
			run_command("gunzip $download_dir/$file_name");
		}
	} else { die "Wrong MD5 for $download_dir/$file_name. Either your download is truncated or we messed up. If this persists over time, please contact us at github.com/jtamames/SqueezeMeta\n" }
        
	}


sub build_wget_command {
	my ($file_name, $host, $download_dir) = @_;

	my $wget_command = "wget -U '' ";
        print("$host\n\n");
	if(!exists($host2urls{$host})) {
		# We can wget using the file name, we assume it is in $host/SqueezeMeta
		$wget_command .= "\"$host/SqueezeMeta/$file_name\"";
		}
	else {  # The wget url is non-standard
		$wget_command .= "\"$host/$host2urls{$host}{$file_name}\"";
		}
	$wget_command .= " -O $download_dir/$file_name";
	return $wget_command;
	}


sub run_command {
	my $command = shift;
	my $ecode = system($command);
	if($ecode) { die "Error running command: $command"; }
	}
1; # return a true value b/c perl is folkloric
