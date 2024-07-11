use strict;
use warnings;

sub get_host {
	my %host2checkfiles = ("http://andes.cnb.csic.es" => "SqueezeMeta/SQMhere",
                               "https://saco.csic.es" => "s/ZFLcptDpWGxAT98/download?path=%2F&files=SQMhere");
	my @goodhosts;
	foreach my $host (keys %host2checkfiles) {
		my $ecode = system("wget -T10 -t1 -O/dev/null -q \"$host/$host2checkfiles{$host}\"");
	        if(!$ecode) { push(@goodhosts, $host); }
	}
	if(!@goodhosts) { die "No host could be reached!" }
	my $host = $goodhosts[ rand @goodhosts ];
	return $host;
}

1; # return a true value b/c perl is folkloric
