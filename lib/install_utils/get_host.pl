use strict;
use warnings;

sub get_host {
	my %preferredhosts = ("http://andes.cnb.csic.es" => "SqueezeMeta/SQMhere");
	my %fallbackhosts  = ("https://saco.csic.es"      => "s/ZFLcptDpWGxAT98/download?path=%2F&files=SQMhere");
	my @goodhosts;
	foreach my $host (keys %preferredhosts) {
		my $ecode = system("wget -T10 -t1 -O/dev/null -q \"$host/$preferredhosts{$host}\"");
	        if(!$ecode) { push(@goodhosts, $host); }
	}
	if(!@goodhosts) {
		        foreach my $host (keys %fallbackhosts) {
				my $ecode = system("wget -T10 -t1 -O/dev/null -q \"$host/$fallbackhosts{$host}\"");
			 	if(!$ecode) { push(@goodhosts, $host); }
			}
	}

	if(!@goodhosts) { die "No host could be reached!" }
	my $host = $goodhosts[ rand @goodhosts ];
	return $host;
}

1; # return a true value b/c perl is folkloric
