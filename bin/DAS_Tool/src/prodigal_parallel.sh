#!/usr/bin/env bash

contigs=$1
proteins=$2
threads=$3
outputbasename=$proteins

grep ">" $contigs | sed 's/>//g;' | shuf > $outputbasename\_DasToolParallel.names
total_lines=$(wc -l <$outputbasename\_DasToolParallel.names)
((lines_per_file = (total_lines + threads - 1) / threads))
split -l $lines_per_file $outputbasename\_DasToolParallel.names $outputbasename\_DasToolParallel.tmp.
rm $outputbasename\_DasToolParallel.names

for i in $outputbasename\_DasToolParallel.tmp.*
do
  pullseq --input $contigs --names $i > $i\.contigs.fa
  prodigal -a $i\.prot.faa -i $i\.contigs.fa -p meta -m -q > /dev/null 2>&1 &
done
wait

cat $outputbasename\_*.prot.faa > $proteins
rm $outputbasename\_DasToolParallel.tmp.*
