
#-- Generic paths

#$installpath="";
$databasepath="$installpath/db";
$extdatapath="$installpath/data";
$scriptdir="$installpath/scripts";   #-- Scripts directory

#-- Paths relatives to the project

$projectname="";
$basedir=".";
$datapath="$basedir/$projectname/data";			#-- Directory containing all datafiles
$resultpath="$basedir/$projectname/results";		#-- Directory for storing results
$extpath="$basedir/$projectname/ext_tables";		#-- Directory for storing tables for further analysis
$tempdir="$basedir/$projectname/temp";			#-- Temp directory
$interdir="$basedir/$projectname/intermediate";			#-- Temp directory
%bindirs=("maxbin","$resultpath/maxbin","metabat2","$resultpath/metabat2");	#-- Directories for bins
%dasdir=("DASTool","$resultpath/DAS/$projectname\_DASTool\_bins");	#-- Directory for DASTool results


#-- Result files

$mappingfile="$datapath/00.$projectname.samples";       #-- Mapping file (samples -> fastq)
$contigsfna="$resultpath/01.$projectname.fasta";        #-- Contig file from assembly
$contigslen="$interdir/01.$projectname.lon";
$rnafile="$resultpath/02.$projectname.rnas";            #-- RNAs from barrnap
$gff_file="$resultpath/03.$projectname.gff";            #-- gff file from prodigal
$aafile="$resultpath/03.$projectname.faa";              #-- Aminoacid sequences for genes
$ntfile="$resultpath/03.$projectname.fna";              #-- Nucleotide sequences for genes
$taxdiamond="$interdir/04.$projectname.nr.diamond";	#-- Diamond result
$cogdiamond="$interdir/04.$projectname.eggnog.diamond";               #-- Diamond result, COGs
$keggdiamond="$interdir/04.$projectname.kegg.diamond";                #-- Diamond result, KEGG
$pfamhmmer="$interdir/05.$projectname.pfam.hmm";      #-- Hmmer result for Pfam
$fun3tax="$resultpath/06.$projectname.fun3.tax";	#-- Fun3 annotations, KEGG
$fun3kegg="$resultpath/07.$projectname.fun3.kegg";	#-- Fun3 annotations, KEGG
$fun3cog="$resultpath/07.$projectname.fun3.cog";	#-- Fun3 annotation, COGs
$fun3pfam="$resultpath/07.$projectname.fun3.pfam";	#-- Fun3 annotation, Pfams
$gff_file_blastx="$resultpath/08.$projectname.gff";            #-- gff file from prodigal & blastx
$fun3tax_blastx="$resultpath/08.$projectname.fun3.tax";	#-- Fun3 annotations prodigal & blastx, KEGG
$fun3kegg_blastx="$resultpath/08.$projectname.fun3.kegg";	#-- Fun3 annotations prodigal & blastx, KEGG
$fun3cog_blastx="$resultpath/08.$projectname.fun3.cog";	#-- Fun3 annotation prodigal & blastx, COGs
$fna_blastx="$interdir/08.$projectname.blastx.fna";	#-- Secuencias nt obtenidas por blastx
$allorfs="$tempdir/09.$projectname.allorfs";         #-- From summary_contigs.pl, allorfs file
$alllog="$interdir/09.$projectname.contiglog";	#-- From summary_contigs.pl, contiglog file (formerly alllog file)
$mapcountfile="$interdir/10.$projectname.mapcount";	#-- From mapsamples.pl, rpkm and coverage counts for all samples
$contigcov="$interdir/10.$projectname.contigcov";     #-- From mapbamsamples.pl, coverages of  for all samples
$mcountfile="$resultpath/11.$projectname.mcount";	#-- From mcount.pl, abundances of all taxa
$mergedfile="$resultpath/13.$projectname.orftable";	#-- GEN TABLE FILE
$bintax="$interdir/17.$projectname.bintax";           #-- From addtax2.pl
$bincov="$interdir/19.$projectname.bincov";           #-- Coverage of bins, from getbins.pl
$bintable="$resultpath/19.$projectname.bintable";       #-- Mapping of contigs in bins, from getbins.pl
$contigsinbins="$interdir/19.$projectname.contigsinbins";
$contigtable="$resultpath/20.$projectname.contigtable"; #-- From getcontigs.pl, CONTIGS TABLE

#-- Datafiles

$coglist="$extdatapath/coglist.txt";    #-- COG equivalence file (COGid -> Function -> Functional class)
$kegglist="$extdatapath/keggfun2.txt";  #-- KEGG equivalence file (KEGGid -> Function -> Functional class)
$pfamlist="$extdatapath/pfam.dat";      #-- PFAM equivalence file
$taxlist="$extdatapath/alltaxlist.txt"; #-- Tax equivalence file 
$nr_db="$databasepath/nr.dmnd";
$cog_db="$databasepath/eggnog";
$kegg_db="$databasepath/keggdb";
$lca_db="$databasepath/LCA_tax/taxid.db";
$bowtieref="$datapath/$projectname.bowtie";                     #-- Contigs formateados para Bowtie
$pfam_db="$databasepath/Pfam-A.hmm";

#-- Variables

$blocksize=8;
$nocog=0;
$nokegg=0;
$nopfam=0;
$euknofilter=0;
$nobins=0;
$doublepass=0;
$cleaning=0;
$cleaningoptions="LEADING:8 TRAILING:8 SLIDINGWINDOW:10:15 MINLEN:30";
$mapper="bowtie";

#-- External software

$metabat_soft="$installpath/bin/metabat2";
$maxbin_soft="$installpath/bin/MaxBin/run_MaxBin.pl";
$spades_soft="$installpath/bin/SPAdes/spades.py";
$barrnap_soft="$installpath/bin/barrnap";
$rdpclassifier_soft="java -jar $installpath/bin/classifier.jar";
$bowtie2_build_soft="$installpath/bin/bowtie2/bowtie2-build";
$bowtie2_x_soft="$installpath/bin/bowtie2/bowtie2";
$bwa_soft = "$installpath/bin/bwa";
$minimap2_soft = "$installpath/bin/minimap2";
$diamond_soft="$installpath/bin/diamond";
$hmmer_soft="$installpath/bin/hmmer/hmmsearch";
$megahit_soft="$installpath/bin/megahit/megahit";
$prinseq_soft="$installpath/bin/prinseq-lite.pl";
$prodigal_soft="$installpath/bin/prodigal";
$cdhit_soft="$installpath/bin/cd-hit-est";
$toamos_soft="$installpath/bin/AMOS/toAmos";
$minimus2_soft="$installpath/bin/AMOS/minimus2";
$checkm_soft="$installpath/bin/checkm";
$minpath_soft="python $installpath/bin/MinPath1.4.py";
$canu_soft="$installpath/bin/canu/canu";
$trimmomatic_soft="java -jar $installpath/bin/trimmomatic-0.38.jar";
$dastool_soft="$installpath/bin/DAS_Tool/DAS_Tool";
