cd $CONCOCT_EXAMPLE
cat $CONCOCT_TEST/reads/Sample*_R1.fa > All_R1.fa
cat $CONCOCT_TEST/reads/Sample*_R2.fa > All_R2.fa
velveth velveth_k71 71 -fasta -shortPaired -separate All_R1.fa All_R2.fa
velvetg velveth_k71 -ins_length 400 -exp_cov auto -cov_cutoff auto  

mkdir contigs
cp velveth_k71/contigs.fa contigs/velvet_71.fa
rm All_R1.fa
rm All_R2.fa


cd $CONCOCT_EXAMPLE
bowtie2-build contigs/velvet_71.fa contigs/velvet_71.fa

for f in $CONCOCT_TEST/reads/*_R1.fa; do
    mkdir -p map/$(basename $f);
    cd map/$(basename $f);
    bash $CONCOCT/scripts/map-bowtie2-markduplicates.sh -ct 1 -p '-f' $f $(echo $f | sed s/R1/R2/) pair $CONCOCT_EXAMPLE/contigs/velvet_71.fa asm bowtie2;
    cd ../..;
done

cd $CONCOCT_EXAMPLE/map
python $CONCOCT/scripts/gen_input_table.py --isbedfiles \
    --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
    ../contigs/velvet_71.fa */bowtie2/asm_pair-smds.coverage \
> concoct_inputtable.tsv
mkdir $CONCOCT_EXAMPLE/concoct-input
mv concoct_inputtable.tsv $CONCOCT_EXAMPLE/concoct-input/

cd $CONCOCT_EXAMPLE/map
python $CONCOCT/scripts/bam_to_linkage.py -m 8 \
    --regionlength 500 --fullsearch \
    --samplenames <(for s in Sample*; do echo $s | cut -d'_' -f1; done) \
    ../contigs/velvet_71.fa Sample*/bowtie2/asm_pair-smds.bam \
> concoct_linkage.tsv
mv concoct_linkage.tsv $CONCOCT_EXAMPLE/concoct-input/

cd $CONCOCT_EXAMPLE
cut -f1,3-26 concoct-input/concoct_inputtable.tsv > concoct-input/concoct_inputtableR.tsv

cd $CONCOCT_EXAMPLE
concoct -c 40 --coverage_file concoct-input/concoct_inputtableR.tsv --composition_file contigs/velvet_71.fa -b concoct-output/

cd $CONCOCT_EXAMPLE
mkdir evaluation-output
Rscript $CONCOCT/scripts/ClusterPlot.R -c concoct-output/clustering_gt1000.csv -p concoct-output/PCA_transformed_data_gt1000.csv -m concoct-output/pca_means_gt1000.csv -r concoct-output/pca_variances_gt1000_dim -l -o evaluation-output/ClusterPlot.pdf

cd $CONCOCT_EXAMPLE
cp $CONCOCT_TEST/evaluation-output/clustering_gt1000_s.csv evaluation-output/
$CONCOCT/scripts/Validate.pl --cfile=concoct-output/clustering_gt1000.csv --sfile=evaluation-output/clustering_gt1000_s.csv --ofile=evaluation-output/clustering_gt1000_conf.csv

$CONCOCT/scripts/ConfPlot.R  -c evaluation-output/clustering_gt1000_conf.csv -o  evaluation-output/clustering_gt1000_conf.pdf

$CONCOCT/scripts/PROKKA_RPSBLAST.sh -f annotations/proteins/velvet_71.faa -p -c 8 -r 1

mkdir $CONCOCT_EXAMPLE/annotations
mkdir $CONCOCT_EXAMPLE/annotations/proteins
mkdir $CONCOCT_EXAMPLE/annotations/cog-annotations 
cp $CONCOCT_TEST/annotations/proteins/* $CONCOCT_EXAMPLE/annotations/proteins/

cd $CONCOCT_EXAMPLE
$CONCOCT/scripts/COG_table.py -g annotations/proteins/velvet_71.gff -b annotations/cog-annotations/velvet_71.out -m $CONCOCT/scgs/scg_cogs_min0.97_max1.03_unique_genera.txt -c concoct-output/clustering_gt1000.csv -e mail@example.com > evaluation-output/clustering_gt1000_scg.tab

cd $CONCOCT_EXAMPLE
$CONCOCT/scripts/COGPlot.R -s evaluation-output/clustering_gt1000_scg.tab -o evaluation-output/clustering_gt1000_scg.pdf

