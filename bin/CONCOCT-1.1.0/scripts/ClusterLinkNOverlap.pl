#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $tFile = '';
my $lFile = '';
my $cFile = '';
my $help = '';
my $quiet = '';
my $outFile = "Clustering_l.csv";

my $USAGE = <<"USAGE";
Usage: ./ClusterLinkOverlapN.pl --cfile=clustering.csv --lfile=linkage.tsv --covfile=coverage.tsv

Regular options:

--ofile    filename -- outputfile for confusion matrix default Conf.csv
--quiet             -- suppress variable names
--help      

USAGE

GetOptions("cfile=s"   => \$tFile,"lfile=s"  => \$lFile, "ofile=s" => \$outFile, "covfile=s" => \$cFile, 'quiet' => \$quiet, 'help' => \$help) or die("Error in command line arguments\n");

if ($help ne '') {print $USAGE;}

die $USAGE unless (($tFile ne '' && $lFile ne '') && $cFile ne '');




#constants
my $minLink = 10; #number of links between two contigs necessary for it to be counted

my $minNLinks = 0; #number of links between two clusters necessary for it to be counted

my $maxPercent = 0.05; #minimum percent of inter cluster links necessary for clustering

my $minOverlap = 0.7; #minimum overlap for clustering


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


my @tranArray = (); #number of links between two clusters

#read link file
open(FILE, $lFile) or die "Can't open $lFile";

my $line = <FILE>;
chomp($line);
my @tokens = split(/\t/,$line);
shift(@tokens);shift(@tokens);
my $nSamples = scalar(@tokens);
my @link_ids = @tokens;
while($line = <FILE>){
  chomp($line);

  my @tokens = split(/\t/,$line);

  my $contig1 = shift(@tokens);
  my $contig2 = shift(@tokens);

  my $link = 0;

  my $lLength = scalar(@tokens);
  
  #sum links between contigs ignoring those in interior
  for(my $i = 0; $i < $lLength; $i=$i+6){
    
    for(my $j = 0; $j < 3; $j++){
      $link+=$tokens[$i + $j];
      if($tokens[$i + $j] > 0){
	#print "$ids[$i + $j] $tokens[$i + $j] $link\n";
      }
    }
    
  }

  if($link > $minLink){

    my $c1 = $hashCluster{$contig1};
    my $c2 = $hashCluster{$contig2};

    #print "$contig1 $contig2 $c1 $c2 $link\n";
    if($c1 ne undef && $c2 ne undef){
      $tranArray[$c1][$c2]++;
    }
 
  }
  #print "$contig1 $contig2 $link\n";
}



#calculate total links associated with each cluster
my @totals = ();
for(my $i = 0; $i < $nClusters; $i++){
  for(my $j = 0; $j < $nClusters; $j++){
    $totals[$i]+=$tranArray[$i][$j];
  }
}

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

my @mmeans = ();

#Normalise by sample sizes
my @sampleSums = ();

for(my $i = 0; $i < $nClusters; $i++){
  for(my $j = 0; $j < $sizes[$i]; $j++){
    my $iid = $clusters[$i][$j];
    for(my $k = 0; $k < $nSites; $k++){
      $sampleSums[$k] += $covs[$hashCov{$iid}][$k];
    }
  }
}

for($i = 0; $i < $nT; $i++){
  for(my $k = 0; $k < $nSites; $k++){
    $covs[$i][$k]/=$sampleSums[$k];
  }
}


#Calculate coverage sum for each cluster
my @csums = ();
for($i = 0; $i < $nClusters; $i++){
  
  my @sums = ();
  for(my $j = 0; $j < $sizes[$i]; $j++){
    for(my $k = 0; $k < $nSites; $k++){
      my $iid = $clusters[$i][$j];
      $sums[$k] += $covs[$hashCov{$iid}][$k];   
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
    $nmeans[$i][$j] = $csums[$i][$j]/$sizes[$i];
  }
}

my @tmeans = ();
for(my $i = 0; $i < $nClusters; $i++){
  my $iTotal = 0.0;

  for(my $k = 0; $k < $nSites; $k++){
    $iTotal += $nmeans[$i][$k];
  }
  $tmeans[$i] = $iTotal;
}

for(my $i = 0; $i < $nClusters; $i++){
  #print "$i ";
  for(my $j = 0; $j < $nSites; $j++){
    $nmeans[$i][$j] /= $tmeans[$i];
  #  print "$nmeans[$i][$j] "
  }
  #print"\n";
}

#Calcaulate coverage overlaps

my @overlaps = ();

calcOverlapMatrix(\@nmeans,$nClusters);


#for($i = 0; $i < $nClusters; $i++){
  #print "$i ";
  #for(my $j = 0; $j < $nClusters; $j++){
   # print "$overlaps[$i][$j] "
  #}
  #print"\n";
#}

my $nClustered = $nClusters; #stores reduced cluster number
my @clustered = (); #new cluster mapping
my @clusterSizes = ();
my @newSizes = ();
#initialise clusters
for(my $i = 0; $i < $nClusters; $i++){
  $clusterSizes[$i] = 1;
  $clustered[$i][0] = $i;
  $newSizes[$i] = $sizes[$i];
}
#writeTranArray(\@tranArray,$nClusters);
#normalise link array and find largest diagonal element
my ($maxI, $maxJ, $maxIJ)= normMaxOff(\@tranArray,\@totals,$nClusters);
print "$maxIJ\n";
while($maxIJ > $maxPercent){
  #cluster pair with largest off diagonal
  my $iString = join(",",@{$clustered[$maxI]});
  my $jString = join(",",@{$clustered[$maxJ]});

  print "$nClustered $maxI->$iString $maxJ->$jString $maxIJ $totals[$maxI] $overlaps[$maxI][$maxJ] $tmeans[$maxI] $tmeans[$maxJ]\n";
 
  push(@{$clustered[$maxJ]},@{$clustered[$maxI]});
  $clusterSizes[$maxJ]+=$clusterSizes[$maxI];
  $clusterSizes[$maxI]=0;

  $newSizes[$maxJ] += $newSizes[$maxI];
  $newSizes[$maxI] = 0;

  $totals[$maxJ]+=$totals[$maxI];
  $totals[$maxI] = 0;

  $tranArray[$maxJ][$maxJ] += $tranArray[$maxI][$maxI];
  for(my $j = 0; $j < $nClustered; $j++){
    $tranArray[$maxJ][$j] += $tranArray[$maxI][$j]; 
  }
  $tranArray[$maxJ][$maxI] = 0;

  for(my $j = 0; $j < $nClustered; $j++){
    $tranArray[$j][$maxJ] += $tranArray[$j][$maxI]; 
  }
  $tranArray[$maxI][$maxJ] = 0;


  #update coverage matrix and overlaps

  for(my $k = 0; $k < $nSites; $k++){
    $csums[$maxJ][$k] += $csums[$maxI][$k];
    $csums[$maxI][$k] = 0;
  }

  for(my $k = 0; $k < $nSites; $k++){
    $nmeans[$maxJ][$k] =$csums[$maxJ][$k]/$newSizes[$maxJ];
    $nmeans[$maxI][$k] = 0.0;
  }

  $tmeans[$maxI] = 0.0; $tmeans[$maxJ] = 0.0;
  
  for(my $k = 0; $k < $nSites; $k++){
    $tmeans[$maxJ] += $nmeans[$maxJ][$k];
  }

  for(my $k = 0; $k < $nSites; $k++){
    $nmeans[$maxJ][$k] /= $tmeans[$maxJ];
  }
  
  #Calcaulate coverage overlaps
  calcOverlapMatrix(\@nmeans,$nClusters);

  $nClustered --;
  
  #writeTranArray(\@tranArray,$nClusters);

  #find newest largest off-diagonal
  ($maxI, $maxJ, $maxIJ)= normMaxOff(\@tranArray,\@totals,$nClusters);
}

print "$nClustered\n";
my $nc = 0;
for(my $i = 0; $i < $nClusters; $i++){
  if($clusterSizes[$i] > 0){
    print "D$i,$clusterSizes[$i],";

    my $dstring = join(",D",@{$clustered[$i]});

    print "E$nc,D$dstring\n";
    $nc++;
  }

}

my $nc = 0;

open(FILE,">$outFile") or die "Can't open $outFile\n";

for(my $i = 0; $i < $nClusters; $i++){
  if($clusterSizes[$i] > 0){
    
    foreach my $oc(@{$clustered[$i]}){
      foreach my $contig(@{$clusters[$oc]}){
	print FILE "$contig,$nc,$oc,$i\n";
      }

    }
    $nc++;
  }

}
close(FILE);


sub writeTranArray(){
  my @array = @{$_[0]};
  my $nClusters = $_[1];
  
  for(my $i = 0; $i < $nClusters; $i++){
    my @printArray = ();

    push(@printArray,@{$array[$i]});
    my $pString = join(",",@printArray);
    print "$i,$pString\n";
  }
}

sub calcOverlapMatrix(){ #Calculates coverage overlaps
  my @nmeans = @{$_[0]};
  my $nC = $_[1];

  for(my $i = 0; $i < $nC; $i++){

    for(my $j = 0; $j < $nC; $j++){
      
      my $overlapIJ = 0.0;
      #print "$j\n";
      for(my $k = 0; $k < $nSites; $k++){
	my $ik = $nmeans[$i][$k];
	my $jk = $nmeans[$j][$k];
	#print "$k $ik $jk\n";
	if($ik > $jk){
	  $overlapIJ += $jk;
	}
	else{
	  $overlapIJ += $ik;
	}
	#print "$overlapIJ\n";
      }
      
      $overlaps[$i][$j] = $overlapIJ;

    }
  }

}

sub normMaxOff(){
  my @tranArray = @{$_[0]};
  my @totals = @{$_[1]};
  my $nClusters = $_[2];

  my $maxIJ = 0.0;
  my $maxI = -1;
  my $maxJ = -1;
  
  my @pArray = ();
  for(my $i = 0; $i < $nClusters; $i++){
  
    for(my $j = 0; $j < $nClusters; $j++){
      if($totals[$i] > 0 && $totals[$j] > 0){
	$pArray[$i][$j] = $tranArray[$i][$j]/$totals[$i];
      }
      else{
	$pArray[$i][$j] = 0;
      }

      if($i != $j){
	if($overlaps[$i][$j] > $minOverlap){
	  if($pArray[$i][$j] > $maxIJ && $totals[$i] > $minNLinks){

	    $maxIJ = $pArray[$i][$j];
	    $maxI = $i;
	    $maxJ = $j;
	  }
	}
      }
    } 
  }

  return ($maxI, $maxJ,$maxIJ);
}
