#!/usr/bin/perl

#-- Recursive curation of bins. Attempts to remove contigs with duplicated marker genes to improve statistics
#-- (c) Javier Tamames, 2018

local $; = '#';
use strict;
use Cwd;
use lib ".";

my $pwd=cwd();
my $project=$ARGV[0];
$project=~s/\/$//; 
if(!$project) { die "Usage: $0 <project> [bin name]\n"; }
my $binreq=$ARGV[1];
if(-s "$project/SqueezeMeta_conf.pl" <= 1) { die "Can't find SqueezeMeta_conf.pl in $project. Is the project name ok?\n"; }
do "$project/SqueezeMeta_conf.pl";
do "$project/parameters.pl";

my $verbose=0;

our($alllog,%bindir,%dasdir,$installpath,$checkm_soft,$numthreads,$datapath,$tempdir,$taxlist);
my %branks=('k','domain','p','phylum','c','class','o','order','f','family','g','genus','s','species');
my(%tax,%consensus,%alltaxa,%goodseeds,%allc,%genes,%count,%contigs,%newcount,%newgenes,%provseeds,%removed,%taxcontig);
my($highscore,$provhighscore,$round,$score,$markers,$changes,$marker,$removed,$binname,$fastaname,$refined,$skip,$finalresult);
my @binlist;

my $markerdir="$datapath/checkm_markers";
my $checktemp="$tempdir/checkm_batch";
my $tempc="$tempdir/checkm_prov.txt";

if(-e $finalresult) { system "rm $finalresult"; }
if(-d $markerdir) {} else { system "mkdir $markerdir"; }
if(-d $checktemp) {} else { system "mkdir $checktemp"; print "Creating $checktemp\n";  }

my $bindir=$dasdir{DASTool};
if($binreq) { push(@binlist,$binreq); }
else {
	opendir(indir1,$bindir) || die "Cannot open binning directory $bindir\n";
	@binlist=(grep/\.fa$|\.fasta$/,readdir indir1);
	closedir indir1;
	}
	

 	#-- Read NCBI's taxonomy 

 	open(infile1,$taxlist) || die "Can't find taxonomy list $taxlist\n";
	print "Reading $taxlist\n";
	while(<infile1>) {
 		 chomp;
 		 next if !$_;
 	 	my @t=split(/\t/,$_);
 	 	$tax{$t[1]}=$t[2];
 	 	}
 	close infile1;

foreach my $bin(sort @binlist) {
	$finalresult="$tempdir/$bin\_checkm_nodupl.txt";
	$fastaname="$bindir/$bin";
	$refined=$bin;
	$refined=~s/\.fa$/\.refined\.fa/;
	$refined=~s/\.fasta$/\.refined\.fasta/;
	if($refined!~/\.refined/) { $refined.=".refined"; }
	$binname="$bindir/$bin\.tax";
	
	$skip=0;

	run_checkm($bin);	#-- Anulando esta llamada se haria sobre el ultimo bin evaluado por checkm
	if(!$skip) { evaluate($bin); }
	if(!$skip) { nebin($bin); }	#-- Extracts new bin and rerun checkm statistics
	}

sub run_checkm {

 my $bin=shift;
 #-- Reading the consensus taxa for the bin
 open(infile1,$binname) || die "Can't open $binname\n";
 print "\nFound bin $bin: ";
 while(<infile1>) { 
 	 chomp;
 	 if($_=~/Consensus/) {
 		 my($cons,$size,$chim,$chimlev)=split(/\t/,$_);
 		 $cons=~s/Consensus\: //;
 		 $size=~s/Total size\: //g;
 		 $consensus{$bin}=$cons;
 		 my @k=split(/\;/,$cons);
		 if($cons!~/k\_/) { print "No consensus tax. Skipping\n"; $skip=1; return; }
 		 print "$cons\n";
		 
 	 
 		 #-- We store the full taxonomy for the bin because not all taxa have checkm markers
 	 
 		 foreach my $ftax(reverse @k) { 
 			 my($ntax,$rank);
 			 if($ftax!~/\_/) { $ntax=$ftax; } else { ($rank,$ntax)=split(/\_/,$ftax); }
 			 $ntax=~s/unclassified //gi;
 			 $ntax=~s/ \<.*\>//gi; 
 			 if($tax{$ntax} && ($rank ne "n") && ($rank ne "s")) { 
 			 push( @{ $alltaxa{$binname} },"$branks{$rank}\_$ntax");
 			 #   print "$m\t$ntax\t$tax{$ntax}\n";
 				 }
 			 }
 		 }
	 else {
	 	my @u=split(/\s+/,$_);
		$u[0]=~s/\>//;
		$taxcontig{$u[0]}=join(" ",@u[4..$#u]);
		}
		 
 	 }
 close infile1;
 my $inloop=1;

 #-- We will find the deepest taxa with checkm markers
 while($inloop) {  
 	 my $taxf=shift(@{ $alltaxa{$binname}  }); 
 	 if(!$taxf) { last; $inloop=0; }
 	 my($rank,$tax)=split(/\_/,$taxf);
 	 $tax=~s/ \<.*//g;
 	 $tax=~s/\s+/\_/g;
 	 # my $rank=$equival{$grank};
 	 if($rank eq "superkingdom") { $rank="domain"; }
 	 print " Using profile for $rank rank : $tax\n";   
 	 $marker="$markerdir/$tax.ms"; 

 	 #-- Use already existing tax profile or create it

 	 if(-e $marker) {} else { 
 		 my $command="export PATH=\"$installpath/bin/pplacer\":\$PATH; $checkm_soft taxon_set $rank $tax $marker > /dev/null 2>&1"; #Override $PATH for external dependencies of checkm. (FPS).
 		 my $ecode = system $command;
 		 if($ecode!=0) { die "Error running command:	$command"; }
 		 }

 	 #-- If it was not possible to create the profile, go for the upper rank
 	 
 	 if(-e $marker) {} else { next; }

 	 my $fastafile=$binname;
 	 $fastafile=~s/\.tax//;
 	 $fastafile=~s/.*\///;
 	 # print ">>> $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null\n";
 	 # system("$checkm_soft analyze -t $numthreads -x $bins{$binname} $marker $bindir $checktemp > /dev/null");
 	 my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft analyze -t $numthreads -x $fastafile $marker $bindir $checktemp > /dev/null 2>&1";
 	 # print "$command\n";
 	 my $ecode = system $command;
 	 if($ecode!=0) { die "Error running command:	$command"; }

 	 my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft qa -t $numthreads $marker $checktemp -f $tempc > /dev/null 2>&1"; #Override $PATH for external dependencies of checkm. (FPS).
 	 # print "$command\n";
	 my $ecode = system $command;
 	 if($ecode!=0) { die "Error running command:	$command"; }
 	 $inloop=0;
 	 }
	}
	
	
sub evaluate() {	
   my $bin=shift;
   (%count,%genes,%contigs,%allc)=();
   open(infile2,"$checktemp/storage/bin_stats_ext.tsv") || die;
    while(<infile2>) {
    	    chomp;
    	    next if !$_;
    	    if($_=~/markers': (\d+)/) { $markers=$1; }
 	  }
    close infile2;

    open(infile3,"$checktemp/storage/marker_gene_stats.tsv") || die;
    while(<infile3>) {
    	    chomp;
    	    next if !$_;
    	    my($bin,$rest)=split(/\t/,$_);
    	    my @marks=split(/\}\,/,$rest);
    	    foreach my $gc(@marks) {
    		    $gc=~s/{|}|[|]|'| //g;
    		    my @mp=split(/\:|\,/,$gc);
    		    my $contigname=shift @mp;
    		    $contigname=~s/\_\d+$//;
    		    $allc{$contigname}=1;
    		    foreach my $gname(@mp) {
    			    next if($gname=~/\[|\]/);
    			    #  print "*$contigname*$gname*\n";
    			    $genes{$gname}{$contigname}=1;
    			    $count{$gname}++;
    			    $contigs{$contigname}{$gname}++;
    			    }
    		    } 
    	    }
    close infile3;
    my($present,$duplicated)=0;
    foreach my $g(keys %genes) {
    	    foreach my $cc(keys %{ $genes{$g} }) {
    	    if($genes{$g}{$cc}) { $present++; last; }
    		    }
    	    }
    my $completion=($present)/$markers;
    foreach my $p(keys %count) {
    	    if($count{$p}>1) { $duplicated++; }
    	    }
    if(!$duplicated) { $duplicated="0"; }
    $score=($present/$markers)-($duplicated/$markers);  			
    printf "Initial: $markers markers (M); $present single (S); $duplicated duplicated (D); %.3f score\n",$score;
    if(!$duplicated) { print "No duplicated markers\n"; $skip=1; return; }
    ($highscore,$provhighscore)=$score;

    #printf "Start: %.3f\n",$highscore;
    %goodseeds=('START',$highscore);
    $round=1;
    recursive($changes);
    }


sub recursive {
	my(%seen,%bad);
	(%newcount,%newgenes)=();
	$changes=1;
	my @allcontigs=keys %allc;
	while($changes) {
		$changes=0;
		my $stringp="";
		my %removed=();
		%provseeds=();
		foreach my $seed(keys %goodseeds) {
			if($verbose) { print "SEEDS:\n"; foreach my $u(keys %goodseeds) { print "--$u"; }; print "\n"; }
			foreach my $nextcontig(@allcontigs) {
				next if($seed=~/$nextcontig/);
				my @oldseeds=split(/\;/,$seed);
				if($verbose) { print " Evaluating contig $nextcontig\n"; }
				push(@oldseeds,$nextcontig);
				sort @oldseeds;
				my $currseed=join(";",sort @oldseeds);
				next if($seen{$currseed});
				next if($bad{$currseed});
				$seen{$currseed}=1;
				%newcount=%count; 
				foreach my $genname(keys %genes) {
					foreach my $contigid (keys %{ $genes{$genname} }) { 
						$newgenes{$genname}{$contigid}=$genes{$genname}{$contigid}; 
						# if($verbose && ($newcount{$genname}>1)) { print "    Starting $genname $contigid:$genes{$genname}{$contigid} $newcount{$genname}\n";  }
 						}                     
					}
					foreach my $removectg(@oldseeds) {
						next if($removectg eq "START");			      
						foreach my $ig(keys %{ $contigs{$removectg} }) { 
							# if($contigs{$removectg}{$ig} && $verbose) { print "    Removed $removectg because having $ig ($newcount{$ig})\n"; }
							$newcount{$ig}-=$contigs{$removectg}{$ig}; 
							$newgenes{$ig}{$removectg}=0; 
							$removed{$removectg}.="$ig;";
							#if($verbose) { print "   $removectg: Removing $ig\n"; }
							}
						}
						
					my($markers,$present,$duplicated,$score)=stats();
					if($verbose) { print "********SCORE= $score $highscore****\n"; }
					 
					              ############ Aqui debe incluirse seleccion por taxonomia o longitud de los contigs que den el mismo score
					
					if($score<$highscore) { $bad{$currseed}=1; } 
					else { 
 						# $provseeds{$currseed}=$score; 
						if(!$duplicated) { $duplicated="0"; }
 						if($verbose) { printf "   $currseed: Score %.5f This contig ($nextcontig): %.5f High %.5f; Removed $removed{$nextcontig}\n",$score,$goodseeds{$seed},$highscore; }
						$stringp="$markers M $present S $duplicated D";
						# print "$string";
						# foreach my $kk(sort keys %newcount) { print "  $kk\n"; }
						$changes=1;  
						# if($score>$provhighscore) { $provhighscore=$score; }
						if($score>$provhighscore) { 
							# print " **HIGH -> "; 
							$provhighscore=$score; 
							%provseeds=(); 
							$provseeds{$currseed}="$nextcontig ($taxcontig{$nextcontig}) SCORE $score"; 
							#foreach my $u(keys %provseeds) { print "$u"; }; 
							#print "\n"; 
							}
     						 }
  						#$string="";

					}
				}
			last if(!$changes);
			print "ROUND $round:"; 
			foreach my $uu(keys %provseeds) { 
				printf " Removing contig $provseeds{$uu} $stringp\n";
				$highscore=$provhighscore;
				%goodseeds=%provseeds;
				$removed=$uu;
				%provseeds=();
				$round++;
				#print "\n\n\nROUND $round\n\n";	
				# die if($round==3);
				}
			}		       
	}


sub nebin {
	my $bin=shift;
	my $into;
	my @rem=split(/\;/,$removed);
	map{ $removed{$_}=1; } @rem;
	open(outfile1,">$bindir/$refined") || die;
	open(infile4,$fastaname) || die;
	while(<infile4>) { 
		chomp;
		if($_=~/^\>([^ ]+)/) {
			my $tcontig=$1;
			if($removed{$tcontig}) { $into=1; } else { $into=0; }	
			}
		next if($into);
		print outfile1 "$_\n";
		}
	close infile4;
	close outfile1;
	
	#-- Rerun checkm for new statistics
	
	print "Recalculating checkm statistics\n";
	
	my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft analyze -t $numthreads -x $refined $marker $bindir $checktemp > /dev/null 2>&1";
 	my $ecode = system $command;
	# print "$command\n";
 	if($ecode!=0) { die "Error running command:	$command"; }
 	my $command = "export PATH=\"$installpath/bin\":\"$installpath/bin/hmmer\":\$PATH; $checkm_soft qa -t $numthreads $marker $checktemp -f $tempc > /dev/null 2>&1"; #Override $PATH for external dependencies of checkm. (FPS).
	# print "$command\n";
 	my $ecode = system $command;
	if(-e $finalresult) { system("cat $finalresult $tempc > $tempdir/gg; mv $tempdir/gg $finalresult;"); }
	else { system("mv $tempc $finalresult"); }
	}	
			


sub stats {
	my($present,$duplicated)=0;
	foreach my $g(keys %newgenes) {
		foreach my $cc(keys %{ $newgenes{$g} }) {
			if($newgenes{$g}{$cc}) { $present++; last; }
 			}
		}
	my $completion=$present/$markers;
	foreach my $p(keys %newcount) {
		if($newcount{$p}>1) { $duplicated++; }
		}
	my $score=($present/$markers)-($duplicated/$markers);			    
	if($verbose) {  printf "M $markers\tP $present\tD $duplicated\tS %.3f\n",$score; }
	return($markers,$present,$duplicated,$score);
	}		    

