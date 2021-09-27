#!/usr/bin/env perl


use strict;
use Getopt::Long; #qw(:config gnu_getopt gnu_compat);
use Data::Dumper;
use threads;

## Globals
our $VERSION = "1.00";
our $AUTHOR = "SQM";
our $PROGRAM = $0;

## Defaults:
my $blastfile;
my $min_identity = 20;
my $min_alilong = 0;
my $min_overlap = 0.8;
my ($maxhits,$numthreads);
my $flgs = {
	    'help'          => 0,
	    'version'       => 0,
	    'best-only'     => 0,
	    'first-hit'     => 0,
	    'no-frames'     => 0,
	    'e-values'      => 0,
	    'no-identical'      => 0,
	    'show-bitscore'      => 0
	   };

GetOptions (
#	    "i|input-file=s"    => \$blastfile,
	    "h|help"          => \$flgs->{help},
	    "v|version"       => \$flgs->{version},
	    "o|overlap:f"     => \$min_overlap,
	    "d|identity:i"    => \$min_identity,
	    "l|alilong:i"    => \$min_alilong,
	    "m|max-hits:i"      => \$maxhits,
            "p|numthreads:i"      => \$numthreads,
	    "b|best-only"     => \$flgs->{'best-only'},
	    "f|first-hit"     => \$flgs->{'first-hit'},
	    "n|no-frames"     => \$flgs->{'no-frames'},
	    "e|e-values"      => \$flgs->{'e-values'},
	    "i|no-identical"      => \$flgs->{'no-identical'},
	    "s|show-bitscore"      => \$flgs->{'show-bitscore'}
	    );
$blastfile = pop @ARGV;
my $tempdir=$blastfile;
my @g=split(/\//,$blastfile);
pop @g;
my $tempdir=join("/",@g);
pop @g;
my $projectdir=join("/",@g);
my $syslogfile="$projectdir/syslog";
open(syslogfile,">>$syslogfile");

#-- Split Diamond file

$numthreads=splitfiles($numthreads);	#-- For setting number of threads equal to number of files, in case $numfile<$numthreads

#-- Launch threads

# print "  Starting multithread blastcollapse in $numthreads threads: ";
print syslogfile "  Starting blastcollapse LCA in $numthreads threads\n";
my $threadnum;
for($threadnum=1; $threadnum<=$numthreads; $threadnum++) {
	# print "$threadnum ";
	my $thr=threads->create(\&current_thread,$threadnum);
	}
# print "\n";
$_->join() for threads->list();

for(my $h=1; $h<=$numthreads; $h++) { 
	open(infile2,"$tempdir/collapsed.$h.m8") || die;
	while(<infile2>) { print; }
	close infile2;
	}

print syslogfile "  Removing temporaty collapsed files in $tempdir\n";
system("rm $tempdir/collapsed.*.m8");
system("rm $tempdir/diamond_collapse.*.m8");
system("rm $tempdir/wc");
close syslogfile;

sub splitfiles {
       # print "  Splitting Diamond file\n";
	my $numthreads=shift;
        print syslogfile "  Splitting Diamond file\n";
        system("wc -l $blastfile > $tempdir/wc");
        open(intemp,"$tempdir/wc");
        my $wc=<intemp>;
        close intemp;
        chomp $wc;
        $wc=~s/\s+.*//;    #-- Number of lines in the diamond result
        my $splitlines=int($wc/$numthreads);
        print syslogfile "Total lines in Diamond: $wc; Allocating $splitlines in $numthreads threads\n";

        my $nextp=$splitlines;
        my ($filelines,$splitorf);
        my $numfile=1;
        # print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
        open(outfiletemp,">$tempdir/diamond_collapse.$numfile.m8");
        open(infile2,$blastfile) || die "Can't open Diamond file $blastfile\n";
        while(<infile2>) {
                $filelines++;
                my @f=split(/\t/,$_);
                if($filelines==$nextp) { $splitorf=$f[0]; }
                if($filelines<=$nextp) { print outfiletemp $_; next; }
                elsif($f[0] ne $splitorf) {
                        close outfiletemp;
                        $numfile++;
                        print syslogfile "Opening file $numfile in line $filelines (estimated in $nextp)\n";
                        open(outfiletemp,">$tempdir/diamond_collapse.$numfile.m8");
                        print outfiletemp $_;
                        $nextp+=$splitlines;
			if($nextp<=$filelines) { $nextp=$filelines+1; }
                        }
                else { print outfiletemp $_; }
                }
        close infile2;
	return $numfile;	#-- For setting number of threads equal to number of files, in case $numfile<$numthreads
        }


sub current_thread {
	my $threadnum=shift;
	$blastfile = "$tempdir/diamond_collapse.$threadnum.m8";
	open(outfile1,">$tempdir/collapsed.$threadnum.m8")  || die;
	# print "#Created by $0, ",scalar localtime,"; Input:$blastfile; "; 
	# print "o:$min_overlap d:$min_identity l:$min_alilong b:$flgs->{'best-only'} f:$flgs->{'first-hit'} n:$flgs->{'no-frames'} e:$flgs->{'e-values'} i:$flgs->{'no-identical'} s:$flgs->{'show-bitscore'}\n"; 
	my $score_field = $flgs->{'e-values'} ? 10 : 11;
	$flgs->{version} && &print_version;
	$flgs->{help} && &print_help;
	$min_overlap /= 100 if ($min_overlap > 1);  ## Normalizing $min_overlap

	my $nextm8 = get_blastm8($blastfile || "piped");  # I/O Closure

	while (my $rec = $nextm8->()){

  	my @st = map {pop @$_
		      } sort {
			$a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]
	     	 } map {
			[
			 (sort {
		  	 $a <=> $b
			 } @{$_}[6,7]),$_
			]
	     	 } map { [split/\t/] } split /\n/,$rec;

	#  my @st = sort {$a->[6] <=> $b->[6] || $a->[7] <=> $b->[7]} map { [split/\t/] } split /\n/,$rec;

 	 my @datalines;
 	 my $prev_query;
 	 my ($coverage_i,$coverage_e);
 	 for (@st){
   		 my @fields=@$_;
   		 my ($query,$identity,$alilong) = @fields[0,2,3];
		 my ($querylong,$hitlong)= @fields[12,13];
   		 my ($frame,$alihitperc);
  		 my ($init1,$end1) = sort {$a<=>$b} @fields[6,7];
		 if($hitlong)  {
		 	$alihitperc=$alilong/$hitlong;
		 	next if(($alihitperc<$min_overlap) && ($init1>=10));   # Exlude partial hits not at the beggining or end of the sequence
			next if(($alihitperc<$min_overlap) && ($end1<=($querylong-10)));
   		 	}
		 next if($identity<$min_identity);
   		 next if($alilong<$min_alilong);
   		 if(($flgs->{'no-identical'}) && ($identity==100)) { next; }
    		#if(!$fields[12]) {      #-- Si ya trae el frame puesto, no la calcules
    		 $frame = ($fields[6]%3);
    		 $frame=3     if ($frame == 0);
    		 $frame *= -1 if ($fields[7]<$fields[6]);
   		 push @fields,$frame;
                #		     }
 
   		 # If we are in the first line:
   		 ($coverage_i,$coverage_e) = ($init1,$end1) if (! defined ($coverage_i) or ! defined ($coverage_e));
   		 $prev_query = $fields[0] if (! defined $prev_query);

   		 my $l1 = $end1 - $init1;
   		 my ($min_e,$max_e) = sort {$a<=>$b} ($end1,$coverage_e);
   		 my $olap = ($min_e - $init1) +1;
   		 my $perc1 = $olap / $l1;
   		 my $perc2 = $olap / ($coverage_e - $coverage_i);
   		 if (($prev_query eq $fields[0]) && (($perc1 >= $min_overlap) || ($perc2 >= $min_overlap))){
		#      push @datalines,[(@fields,$frame)];
     		 push @datalines,[@fields];
     		 $coverage_e += ($end1 - $coverage_e) if ($end1 > $coverage_e);
   		 } else {
     		 my $line=process (\@datalines);
		 print outfile1 $line;
     		 @datalines = [@fields];
		#      @datalines = [(@fields,$frame)];
     		 ($coverage_i,$coverage_e) = ($init1,$end1);
     		 $prev_query = $fields[0];
   		 }
 		 }
  		my $line=process (\@datalines);
		print outfile1 $line;	
	}
	close outfile1;
}

#### Esta versión de "process" sólo saca 1 record por hit (tiene en cuenta empates)
sub process_only1
  {
    my ($refarray) = @_;
    my $by_frames = build_struct($refarray);
	my $line;
	my $score_field = $flgs->{'e-values'} ? 10 : 11;

    for my $fr (keys %$by_frames){
      my @best_recs = ();
      my $best_ev;
      for my $r (@{${$by_frames}{$fr}}){
	my $curr_ev = $r->[$score_field];
	$best_ev = $curr_ev if (! defined $best_ev);
	if ($score_field == 10){
	  if ($best_ev > $curr_ev){
	    @best_recs = ($r);
	    $best_ev = $curr_ev;
	  } else {
	    push @best_recs,$r if ($best_ev == $curr_ev);
	  }
	} else {
	  if ($best_ev < $curr_ev){
	    @best_recs = ($r);
	    $best_ev = $curr_ev;
	  } else {
	    push @best_recs,$r if ($best_ev == $curr_ev);
	  }
	}
      }
	$line.=join ("\t",@$_),"\n" for (@best_recs);
      # print join ("\t",@$_),"\n" for (@best_recs);
    }
 return $line;
  }

sub process
  {
    goto (&process_only1) if $flgs->{'best-only'};
    my ($refarray) = @_;
    my $by_frames = build_struct($refarray);
	my $line;
	my $score_field = $flgs->{'e-values'} ? 10 : 11;
    my ($qinit,$qend,$hinit,$hend);
    my ($tosw,$ttem);

    for my $fr (keys %$by_frames){
      my @sorted;
      if ($score_field == 10){
	@sorted = sort {$a->[$score_field] <=> $b->[$score_field]} @{${$by_frames}{$fr}};
      } else {
	@sorted = sort {$b->[$score_field] <=> $a->[$score_field]} @{${$by_frames}{$fr}};
      }
      my @recj = @{shift @sorted};
      if(!$qinit) { ($qinit,$qend,$hinit,$hend)=@recj[6..9]; }
      $recj[11]=~s/\s+//g;
       if($flgs->{'show-bitscore'}) { $recj[1] .= "|$recj[11]|$recj[2]"; }
      my ($aux_min,$aux_max) = sort {$a <=> $b} @recj[6,7];
#      my $aux_max = $recj[7];
	my $numhits=1;
      for my $r (@sorted){
#	print STDERR Dumper $r;
	$numhits++;
	if($maxhits && ($numhits>=$maxhits)) { last; }
	my ($rec_min,$rec_max) = sort {$a <=> $b} @{$r}[6,7];
	my ($rec_minH,$rec_maxH) = sort {$a <=> $b} @{$r}[8,9];
        $r->[11]=~s/\s+//g;
	$recj[1] .= ";$r->[1]";
	if($flgs->{'show-bitscore'}) { $recj[1] .= "|$r->[11]|$r->[2]"; }
	if($flgs->{'first-hit'}) {  
	 $recj[6]=$qinit;
	 $recj[7]=$qend;
	 $recj[8]=$hinit;
	 $recj[9]=$hend;
	                }
	else {
	if($recj[6]>$recj[7]) { $tosw=1; }
 	$recj[6] = $rec_min if ($rec_min < $recj[6]);
 	$recj[7] = $rec_max if ($rec_max > $recj[7]);
 	$recj[8] = $rec_minH if ($rec_minH < $recj[8]);
 	$recj[9] = $rec_maxH if ($rec_maxH > $recj[9]);
	     }
      }
      if($tosw) { $ttem=$recj[6]; $recj[6]=$recj[7]; $recj[7]=$ttem; }
      if($recj[6]<$recj[7]) { $recj[0].="\_$recj[6]-$recj[7]"; } else { $recj[0].="\_$recj[7]-$recj[6]"; }
      $line.=join ("\t",@recj). "\n";
    }
 return $line;
  }


sub build_struct
  {
    my ($refarray) = @_;
    my %by_frames;
    for my $rec (@$refarray){
      if ($flgs->{'no-frames'}){
	push @{$by_frames{'-'}},$rec;
      } else {
	push @{$by_frames{@{$rec}[$#{$rec}]}},$rec;
      }
    }
    return \%by_frames;
  }

sub get_blastm8 {
  my ($F) = @_;
  my $m8;
  if ($F eq "piped"){
    $m8 = *STDIN;
  } else {
    open $m8, "<", $F or die "$F: $!";
  }
  my $last_seen;
  my $rec = '';
  return sub {
    while (1){
      if (eof $m8) { close $m8 unless ($F eq "piped"); return undef }
      return undef if (eof $m8);
      $_ = <$m8>;
      return $rec.$_ if (eof $m8);
      my @f = split /\t/;
      $last_seen = $f[0] if (! defined $last_seen);
      if ($last_seen eq $f[0]){
	$rec .= $_;
      } else {
	my $retrec = $rec;
	$rec = $_;
	$last_seen = $f[0];
	return $retrec;
      }
    }
  }
}

sub print_help
{
    print <<__EOH__;

Usage: $PROGRAM [options] <blast-output>

Options:

           [ -o | --overlap ]
                        Minimum percentage of overlap to collapse hits (defaults to 0.8)

	   [ -b | --best-only ]
	                Give only the best hit or concatenate the headers of the hits

	   [ -f | --first-hit ]
	                Get all data from first hit, do not adjust positions

           [ -e | --e-values ]
                        Use the e-values for selecting the best hit (by default the bit score is used).
                        If -b is present, this criteria is used for selecting the best hit.
                        If not, this criteria is used for sorting the headers to concatenate.


           [ -d | --identity ]
                        Identities below this value are discarded (default: 20)

           [ -l | --alilong ]
                        Alignmenth lengths below this value are discarded (default: 0)

           [ -n | --no-frames ]
                        By default, hits in each frame are treated independently unless this option is present
   
           [ -i | --no-identical ]
                        Excludes identical hits

           [ -s | --show-bitscore ]
                        Show the corresponding bitscore for all the collapsed hits
			
           [ -m | --maxhits ]
                        Maximum number of hits to collapse

           [-h | --help]
	                Print this help and exit

           [-v | --version]
	                Print program version and exit


__EOH__

exit();
}
sub print_version
{
    print <<__EOV__;
$PROGRAM v$VERSION
__EOV__
exit();
}
