#!/usr/bin/env bash

function display_help() {
    echo " "
    echo "Scaffolds2Bin_to_Fasta: Creates fasta files of bins given a metagenomics assembly (fasta format) and binning prediction (scaffolds2bin file) in tabular format."
    echo " (DAS Tool helper script)"
    echo " "
    echo "Usage: Scaffolds2Bin_to_Fasta.sh -i <scaffolds2bin.tsv> -a <assembly.fasta> -o <output_folder>"
    echo " "
    echo "   -i, --scaffolds2bin        Scaffolds to bin file."
    echo "   -a, --assembly             Assembly in fasta format."
    echo "   -o, --output_folder        Output folder. [default: <scaffolds2bin.tsv>_fasta]"
    echo "   -d, --delimiter            Delimiter of scaffolds to bin file. [default: \t]"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

[ $# -eq 0 ] && { display_help ; exit 1; }

delimiter="\t"

while [ "$1" != "" ]; do
    case $1 in
        -i | --scaffolds2bin )  shift
                                scaffolds2bin=$1
                                ;;
        -a | --assembly )       shift
                                assembly=$1
                                ;;
        -o | --output_folder )  shift
                                output_folder=$1
                                ;;
        -d | --delimiter )      shift
                                delimiter=$1
                                ;;
        -h | --help )           display_help
                                exit
                                ;;
        * )                     display_help
                                exit 1
    esac
    shift
done

if [ -z "$output_folder" ]
then
  output_folder=$scaffolds2bin\_fasta
fi

if [ -d "$output_folder" ]; then
  rm $output_folder\/*.ctgtmp
else
  mkdir $output_folder
fi

if [ "$delimiter" != "\t" ]; then
	perl -pe "s/$delimiter/\t/g;" $scaffolds2bin > $output_folder\/scaffolds2bin.tsv
	scaffolds2bin=$output_folder\/scaffolds2bin.tsv
fi


echo extracting bins to $output_folder
while read scaffold bin
do
echo $scaffold >> $output_folder\/$bin\.ctgtmp
done < $scaffolds2bin

for i in $output_folder\/*.ctgtmp
do
# bin_id=$(echo $i | perl -pe 's/\.ctgtmp//g')
bin_id=$(echo $i | sed 's/\.ctgtmp//g')
pullseq -i $assembly -n $i > $bin_id\.fasta
rm $i
done
