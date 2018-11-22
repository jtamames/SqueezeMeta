#!/usr/bin/perl -w
use strict;

#getsam($ARGV[0], $ARGV[1]);

# getsam(input sam file, output)
sub getsam
{
	my $MAX_MISMATCH = 10;

	my $sam = $_[0];
	my $out = $_[1];
	my $line;
	my @arr;
	my $tmp1;
	my $tmp2;
	my $i;
	my $j;
	my $mismatch;

	my %lenhash = ();
	my %sumhash = ();

	print "Reading SAM file to estimate abundance values...\n";
	open(ABUNDFILE, "<$sam") || die "Cannot open file $sam\n";
	while(defined($line = <ABUNDFILE>))
	{
		chomp($line);
		@arr = split(/\t/, $line);
		if ($arr[0] =~ /^\@/)
		{
			if ($arr[0] eq "\@SQ")
			{
				$tmp1 = substr($arr[1], 3);
				$tmp2 = substr($arr[2], 3);
				$lenhash{$tmp1} = $tmp2;
				$sumhash{$tmp1} = 0;
			}
		}
		elsif ($arr[2] ne "*")
		{
			$j = @arr;
			$mismatch = $MAX_MISMATCH + 1;
			for ($i = 10; $i < $j; $i++)
			{
				if ($arr[$i] =~ /XM:i:([0-9]+)/)
				{
					$mismatch = $1;
					last;
				}
			}
			$mismatch = $mismatch + getMismatchFromCIAGR($arr[5]);
			if ($mismatch <= $MAX_MISMATCH)
			{
				$sumhash{$arr[2]} += getCIGARLen($arr[5]);
			}
		}
	}
	close(ABUNDFILE);

	open(ABUNDOUTPUT, ">$out");
	foreach $line (sort keys %sumhash)
	{
		$tmp1 = $sumhash{$line} / $lenhash{$line};
		print ABUNDOUTPUT "$line\t$tmp1\n";
	}
	close(ABUNDOUTPUT);
}

sub getCIGARLen
{
	my @subarr = split(/([A-Z])/, $_[0]);
	my $m = @subarr;
	my $n;
	my $sum = 0;
	for ($n = 0; $n < $m - 1; $n++)
	{
		if ($subarr[$n + 1] eq "M")
		{
			$sum = $sum + $subarr[$n];
		}
	}
	return($sum);
}

sub getMismatchFromCIAGR
{
	my @CIAGRarr = split(/([A-Z=])/, $_[0]);
	my $CIAGRlen = @CIAGRarr;
	my $start_gap = 0;
	my $end_gap = 0;
	my $m;
	my $ret = 0;
	for ($m = 1; $m < $CIAGRlen; $m = $m + 2)
	{
		if ($CIAGRarr[$m] eq "I" || $CIAGRarr[$m] eq "D")
		{
			$ret = $ret + $CIAGRarr[$m - 1];
		}
		if ($CIAGRarr[$m] eq "S")
		{
			if ($m == 1 && $CIAGRarr[$m - 1] > 1)
			{
				$start_gap = 1;
			}
			elsif ($m == $CIAGRlen - 1 && $CIAGRarr[$m - 1] > 1)
			{
				$end_gap = 1;
			}
		}
	}
	if ($start_gap == 1 && $end_gap == 1)
	{
		#$ret = 10;
	}
	return $ret;
}

1;

