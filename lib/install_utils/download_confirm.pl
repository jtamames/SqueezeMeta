use strict;
use warnings;
use Digest::MD5 qw(md5_hex);

sub download_confirm {
	my ($file_name, $md5_name, $remote_path, $download_dir) = @_;

	system("rm $download_dir/$file_name > /dev/null 2>&1");
        system("rm $download_dir/$md5_name  > /dev/null 2>&1");
	
	run_command("wget -U '' -P $download_dir $remote_path/$file_name");
	run_command("wget -U '' -P $download_dir $remote_path/$md5_name");

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

sub run_command {
	my $command = shift;
	my $ecode = system($command);
	if($ecode) { die "Error running command: $command"; }
	}
1; # return a true value b/c perl is folkloric
