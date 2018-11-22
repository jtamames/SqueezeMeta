use strict;

#sepReads($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3], $ARGV[4]);

# main(read file, reads number (1 for first reads file, 2 for second, etc.), contig header (output destination of MaxBin), sam file, is_fastq--1 means yes)
sub sepReads
{
	my $MAX_MISMATCH = 10;
	my $readsfile = $_[0];
	my $readscount = "reads" . $_[1];
	my $contigheader = $_[2];
	my $samfile = $_[3];
	my $isfastq = $_[4];
	my $dir;
	my $f;
	my $f2;
	my $c;
	my $header;
	my $seq;
	my $tmpwrite;
	my $dest;
	my $i;
	my $j;
	my $count;
	my $totalline;
	my $line;
	my @arr;
	my %filehash;
	my %curfilehash;
	my %contighash;
	my %readhash;
	my %tmphash;
	my $mis;
	my $indel;
	my @filearr;
	my %writehash;
	my $isnoclass;

	# parse contig files through path
	@arr = split(/(\/)/, $contigheader);
	$j = @arr;
	$dir = "";
	for ($i = 0; $i < $j - 2; $i++)
	{
		$dir = $dir . $arr[$i];
	}
	$header = $arr[$j - 1];
	if ($dir eq "")
	{
		$dir = ".";
	}
	opendir(DIR, $dir);
	while($f = readdir(DIR))
	{
		if ($f =~ /$header.([0-9]+).fasta/)
		{
			$i = $1;
			open(FILE, "<$dir\/$f");
			while(defined($line = <FILE>))
			{
				if ($line =~ /^>/)
				{
					chomp($line);
					$line = substr($line, 1);
					$contighash{$line} = $i;
				}
			}
			close(FILE);

			my $tmp;
			if (-e "$contigheader.reads.$i")
			{
				open($tmp, ">>$contigheader.reads.$i");
			}
			else
			{
				open($tmp, ">$contigheader.reads.$i");
			}
			$filehash{$i} = $tmp;
		}
	}
	closedir(DIR);
	my $noclass;
	open($noclass, ">$contigheader.reads.noclass");

	# parse sam file
	print "Reading sam file $samfile...\n";
	open(FILE, "<$samfile") || die "Cannot open SAM file $samfile\n";
	while(defined($line = <FILE>))
	{
		@arr = split(/\t/, $line);
		if ($arr[0] ne "\@SQ" && $arr[0] ne "\@PG" && $arr[0] ne "\@HD" && exists $contighash{$arr[2]})
		{
			$mis = $MAX_MISMATCH + 1;
			for ($i = 11; $i < 18; $i++)
			{
				if (!(defined($arr[$i])))
				{
					last;
				}
				elsif ($arr[$i] =~ /XM:i:([0-9]+)/)
				{
					$mis = $1;
					last;
				}
			}
			$indel = getMismatchFromCIAGR($arr[5]);
			$mis = $mis + $indel;
			if ($mis > $MAX_MISMATCH)
			{
				next;
			}

			if (!(exists $readhash{$arr[0]}))
			{
				$readhash{$arr[0]} = $arr[2];
			}
			else
			{
				if (index($readhash{$arr[0]}, $arr[2]) == -1)
				{
					$readhash{$arr[0]} = $readhash{$arr[0]} . ";" . $arr[2];
				}
			}
		}
	}
	close(FILE);

	# get partial reads from parsed sam file
	print "Reading reads file $readsfile...\n";
	open(FILE, "<$readsfile");
	@filearr = ();
	$isnoclass = 0;
	$count = 0;
	$totalline = 0;
	while(defined($line = <FILE>))
	{
		$totalline++;
		if (($isfastq == 0 && $line =~ /^>/) || ($isfastq == 1 && $totalline % 4 == 1))
		{
			$count++;
			$i = index($line, " ");
			if ($i == -1)
			{
				$header = substr($line, 1);
			}
			else
			{
				$header = substr($line, 1, $i - 1);
			}
			chomp($header);
			if ($count % 2 == 1) # first pair
			{
				if (exists $readhash{$header})
				{
					$isnoclass = 0;
					@arr = split(/;/, $readhash{$header});
					@filearr = ();
					$tmpwrite = ">" . $readscount . "_" . $header . "\n";
					foreach $c (@arr)
					{
						$f = $filehash{$contighash{$c}};
						push(@filearr, $f);
					}
				}
				else
				{
					$isnoclass = 1;
					@filearr = ();
					$tmpwrite = ">" . $readscount . "_" . $header . "\n";
				}
			}
			else # second pair
			{
				if (exists $readhash{$header})
				{
					$isnoclass = 0;
					@arr = split(/;/, $readhash{$header});
					foreach $c (@arr)
					{
						$f = $filehash{$contighash{$c}};
						push (@filearr, $f);
					}
					%writehash = ();
					foreach $f (@filearr)
					{
						if (!(exists $writehash{$f}))
						{
							$writehash{$f} = 0;
							print $f $tmpwrite;
							print $f ">" . $readscount . "_" . $header . "\n";
						}
					}
				}
				else
				{
					if ($isnoclass == 1)
					{
						print $noclass $tmpwrite;
						print $noclass ">" . $readscount . "_" . $header . "\n";
					}
					else
					{
						$isnoclass = 0;
						%writehash = ();
						foreach $f (@filearr)
						{
							if (!(exists $writehash{$f}))
							{
								$writehash{$f} = 0;
								print $f $tmpwrite;
								print $f ">" . $readscount . "_" . $header . "\n";
							}
						}
					}
				}
			}
		}
		elsif (($isfastq == 1 && $totalline % 4 == 2) || ($isfastq == 0))
		{
			if ($count %2 == 1) # First pair
			{
				$tmpwrite = $tmpwrite . $line;
			}
			else
			{
				if ($isnoclass == 1)
				{
					print $noclass $line;
				}
				else
				{
					%writehash = ();
					foreach $f (@filearr)
					{
						if (!(exists $writehash{$f}))
						{
							$writehash{$f} = 0;
							print $f $line;
						}
					}
				}
			}
		}
	}

	# Close all files
	foreach $f (keys %filehash)
	{
		close($f);
	}
	close($noclass);
}

=cut
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
			#$ret = $ret + $CIAGRarr[$m - 1];
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
=cut

1;
