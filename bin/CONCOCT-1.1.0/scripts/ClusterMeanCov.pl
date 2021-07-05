#!/usr/bin/perl

use strict;
use Getopt::Long;

my $tFile = '';
my $fFile = '';
my $cFile = '';
my $help = '';
my $quiet = '';

my $USAGE = <<"USAGE";
Usage: ./ClusterMeanCov.pl --cfile=clustering.csv --covfile=coverage.tsv --ffile=fastafile

Regular options:

--quiet             -- suppress variable names
--help      

USAGE

GetOptions("cfile=s"   => \$tFile,"covfile=s" => \$cFile, "ffile=s" => \$fFile, 'quiet' => \$quiet, 'help' => \$help) or die("Error in command line arguments\n");

if ($help ne '') {print $USAGE;}

die $USAGE unless ($tFile ne '' && $cFile ne '');

#read in clusterings
my %hashCluster = {}; #map contig id to cluster
my @clusters = (); #contigs in each cluster
my @sizes = (); #number of contigs in each cluster
my $maxt = 0; #largest cluster index observed

open(FILE, $tFile) or die "Can't open $tFile";

while(my $line = <FILE>){
  
  chomp($line);
  
  my @tokens = split(/,/,$line);

  my $name = $tokens[0];
  my $cluster = $tokens[1];
  
  $clusters[$cluster][$sizes[$cluster]]=$name;
  $sizes[$cluster]++;
  
  $hashCluster{$name} = $cluster;
  #print "$name $cluster\n";
  if($cluster > $maxt){
    $maxt = $cluster;
  }
}

close(FILE);

#number of clusters in data set
my $nClusters = $maxt + 1;

my @covs = ();
my %hashCov = {};
my $i = 0;
my $nSites = 0;

open(FILE, $cFile) or die "Can't open $cFile\n";

my $line = <FILE>;
chomp($line);
my $header=$line;
my @header = split(/\t/,$line);
while($line = <FILE>){
  chomp($line);
  my @tokens = split(/\t/,$line);
  my $id = shift(@tokens);

  if($i == 0){
    $nSites = scalar(@tokens);
  }
  $hashCov{$id} = $i;
  push(@covs,\@tokens);
  $i++;
}
my $nT = $i;

close(FILE);

my @Seq = ();
my @id       = ();

my $count = 0;

open(FILE, $fFile) or die "Can't open $fFile\n";

my $seq = "";

while($line = <FILE>){
    chomp($line);

    if($line =~ />(.*)/){

        $id[$count] = $1;

        if($seq ne ""){
            $Seq[$count - 1] = $seq;

            $seq = "";
        }

        $count++;
    }
    else{
        $seq .= $line;
    }
}

$Seq[$count - 1] = $seq;
my $total = $count;

my %hashLengths = {};

for($i = 0; $i < $total; $i++){
    $hashLengths{$id[$i]} = length($Seq[$i]);
}


my @mmeans = ();

#Normalise by sample sizes
my @sampleSums = ();

for(my $i = 0; $i < $nClusters; $i++){
  for(my $j = 0; $j < $sizes[$i]; $j++){
    my $iid = $clusters[$i][$j];
    for(my $k = 0; $k < $nSites; $k++){
      $sampleSums[$k] += $covs[$hashCov{$iid}][$k]*$hashLengths{$iid};
    }
  }
}

#for($i = 0; $i < $nT; $i++){
 # for(my $k = 0; $k < $nSites; $k++){
  #  $covs[$i][$k]/=$sampleSums[$k];
  #}
#}


#Calculate coverage sum for each cluster
my @csums = ();
my @clengths = ();
for($i = 0; $i < $nClusters; $i++){
  
  my @sums = ();
  for(my $j = 0; $j < $sizes[$i]; $j++){
    my $iid = $clusters[$i][$j];
    $clengths[$i] += $hashLengths{$iid};	
    for(my $k = 0; $k < $nSites; $k++){
      $sums[$k] += $covs[$hashCov{$iid}][$k]*$hashLengths{$iid};   
    }
  }
  
  for(my $k = 0; $k < $nSites; $k++){
    $csums[$i][$k] = $sums[$k];
  }
}

#calculate coverage means
my @nmeans = ();

for(my $i = 0; $i < $nClusters; $i++){
  for(my $j = 0; $j < $nSites; $j++){
    $nmeans[$i][$j] = $csums[$i][$j]/$clengths[$i];
  }
}

shift(@header);
my $hstring = join(",",@header);
print "Sample,$hstring\n";
for(my $i = 0; $i < $nClusters; $i++){
  my @temp = @{$nmeans[$i]};
  my $dstring = join(",",@temp);

  print "D$i,$dstring\n";
}

