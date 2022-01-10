use strict;
use warnings;

sub get_host {
	my @hosts = ("http://silvani.cnb.csic.es", "http://andes.cnb.csic.es");
	my @goodhosts;
	for my $host (@hosts) {
		my $ecode = system("wget -T10 -t1 -O/dev/null -q $host/SqueezeMeta/SQMhere");
	        if(!$ecode) { push(@goodhosts, $host); }
	}
	if(!@goodhosts) { die "No host could be reached!" }
	my $host = $goodhosts[ rand @goodhosts ];
	return $host;
}

1; # return a true value b/c perl is folkloric
