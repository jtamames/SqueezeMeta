#!/usr/bin/env bash

ctg2bin=$1
contigs=$2
bin_folder=$3

while read contig bin
do
echo $contig >> $bin_folder\/$bin\.DasToolBin.tmp
done < $ctg2bin

for i in $bin_folder\/*.DasToolBin.tmp
do
bin_id=$(echo $i | sed 's/\.DasToolBin\.tmp//g')
pullseq -i $contigs -n $i > $bin_id\.fa
rm $i
done
