#!/bin/env bash 

function display_help() {
    echo " "
    echo "Fasta_to_Scaffolds2Bin: Converts genome bins in fasta format to scaffolds-to-bin table."
    echo " (DAS Tool helper script)"
    echo " "
    echo "Usage: Fasta_to_Scaffolds2Bin.sh -e fasta > my_scaffolds2bin.tsv"
    echo " "
    echo "   -e, --extension            Extension of fasta files. (default: fasta)"
    echo "   -i, --input_folder         Folder with bins in fasta format. (default: ./)"
    echo "   -h, --help                 Show this message."
    echo " "
    echo " "
    exit 1
}

# [ $# -eq 0 ] && { display_help ; exit 1; }

extension="fasta"
folder="."

while [ "$1" != "" ]; do
    case $1 in
        -e | --extension )      shift
                                extension=$1
                                ;;
        -i | --input_folder )   shift
                                folder=$1
                                ;;
        -h | --help )           display_help
                                exit
                                ;;
        * )                     display_help
                                exit 1
    esac
    shift
done

for i in $folder\/*.$extension
do
binname=$(echo $(basename $i) | sed "s/\\.$extension//g")
grep ">" $i | perl -pe "s/\n/\t$binname\n/g" | perl -pe "s/>//g" 
done


