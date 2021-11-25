#-- Internal parameters for scripts, that are NOT specified by the user when running SQM

#$blocksize=8;	#-- STEP4: Blocksize for Diamond (also for STEP8)
$evaluetax4=1e-03;	#-- STEP4: evalue for Diamond, for tax only (also for STEP8)
$minidentax4=40;	#-- STEP4: Minimum identity for Diamond, for tax only (also for STEP8)
$evaluefun4=1e-03;	#-- STEP4: evalue for Diamond, for functions only (also for STEP8)
$minidenfun4=30;	#-- STEP4: Minimum identity for Diamond, for functions only (also for STEP8)


$evaluehmmer5=1e-10;	#-- STEP5: Minimum evalue for hmmer

$scoreratio6=0.8;   #-- STEP6/8: Ratio first score/currsent score for the hit to be considered
$diffiden6=10;       #-- STEP6/8: Maximim identity difference with the first
$flex6=0.2;           #-- STEP6/8: Allows this PERCENTAGE (if less than one) or NUMBER (if greater than one) of hits from different taxa than LCA
$minhits6=2;        #-- STEP6/8: Minimum number of hits for the taxa (if there is only one valid hit, this value sets to one automatically
$noidfilter6=1;	#-- STEP6/8: Set to 1, creates a new set of results with no identity filters

$mindif7=0.1;		     #-- STEP7: Percentage of difference in scores for assigning bestaver
$maxhits7=3;	     #-- STEP7: Maximum number of hits for averaging
$minolap7=30;		     #-- STEP7: Minimum alignment percentage between query and hit

$mingenes9=1;			     #-- STEP9: Minimum number of genes for the contig 
$minconsperc_asig9=0.7;      #-- STEP9: Ratio genes for the taxon/sum(genes all taxa). Therefore it only considers assigned genes
$minconsperc_total9=0.5;     #-- STEP9: Ratio genes for the taxon/number of genes. Therefore it considers all (assigned+unassigned) genes
			     #-- THIS ONE can be set by arguments, because of its influence when using long reads
	
$keepsam10=1;  #-- STEP10: Set to one, it keeps SAM files. Set to zero, it deletes them when no longer needed

$minraw12=200;  	     #-- STEP12: Minimum number of raw counts to be included in the STAMP files

$seqsinfile13=0;     #-- STEP13: Put sequences in the output table (0=no, 1=yes)

$maxchimerism14=0.1;   #--STEP14:  Threshold for excluding chimeric contigs
$mingenes14=1;	     #--STEP14:  Threshold for excluding small contigs (few genes than this)
$smallnoannot14=1;     #-- STEP14: For excluding contigs with just one gene an no annotation

$score_tres15=0;    #-- STEP15: Score threshold for keeping bins (proxy for level of completeness)

$mincontigs16=1;  	     #-- STEP16: Minimum number of contigs for the bin 
$minconsperc_asig16=0.6;       #-- STEP16: Ratio contigs for the taxon/sum(genes all taxa). Therefore it only considers assigned contigs
$minconsperc_total16=0.3;      #-- STEP16: Ratio contigs for the taxon/number of contigs. Therefore it considers all (assigned+unassigned) contigs

$minsize17=20000;  #--STEP17: Minimum size of a bin to be considered

$minfraction20=0.1;    #--STEP20: Minimum percentage of genes from a pathway to be present
