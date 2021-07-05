#!/bin/bash
# ***************************************************************
# Name:      RPSBLAST.sh
# Purpose:   This script integrates with PROKKA and runs PROKKA_XXXXXXXX.faa with rpsblast against the  Cog database
# Version:   0.1-inodb_edit
# Authors:   Umer Zeeshan Ijaz (Umer.Ijaz@glasgow.ac.uk)
#                 http://userweb.eng.gla.ac.uk/umer.ijaz
# Contr:     Ino de Bruijn (ino.debruijn@scilifelab.se)
#
# Created:   2014-01-16
# License:   Copyright (c) 2014 Computational Microbial Genomics Group, University of Glasgow, UK
#
#            This program is free software: you can redistribute it and/or modify
#            it under the terms of the GNU General Public License as published by
#            the Free Software Foundation, either version 3 of the License, or
#            (at your option) any later version.
#
#            This program is distributed in the hope that it will be useful,
#            but WITHOUT ANY WARRANTY; without even the implied warranty of
#            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#            GNU General Public License for more details.
#
#            You should have received a copy of the GNU General Public License
#            along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **************************************************************/

HELPDOC=$( cat <<EOF
Usage:
    bash `basename $0` -f <fasta_file.faa> [options]

Options:
    -p Parallelize using GNU Parallel flag
    -c Number of cores to use (Default: 10)
    -t Number of threads to use (Default: 1)
    -r Number of reference matches (Default: 500)
    -e Expectation value (E) threshold for saving hits (Default: 0.00001)

EOF
)


# = Enable FP support ============== #
# By default, there is limited capability in bash to handle floating point
# operations. In this script bc is used to calculate the floating point operations.
# $float_scale parameter specifies the precision of the floating point.
# Reference: http://www.linuxjournal.com/content/floating-point-math-bash

float_scale=5
# Evaluate a floating point number expression.

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


# Evaluate a floating point number conditional expression.

function float_cond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}

function ceil () {
  echo "define ceil (x) {if (x<0) {return x/1} \
        else {if (scale(x)==0) {return x} \
        else {return x/1 + 1 }}} ; ceil($1)" | bc;
 }

# =/Enable FP support ============== #


# = Parameters to set ============== #
if [ ! -e $COGSDB_DIR ]; then
    echo "COGSDB_DIR doesn't exist. Set to directory of the COG database" >&2
    exit 1
fi
FASTA_FILE=""; # This field should be empty
PARALLELIZE_FLAG=0
NUMBER_OF_CORES=10
NUMBER_OF_THREADS=1
NUMBER_OF_REFERENCE_MATCHES=500
EVALUE=0.00001
# =/Parameters to set ============== #

# Parse options
while getopts ":phc:f:e:r:t:" opt; do
    case $opt in
        p)
            PARALLELIZE_FLAG=1
            ;;
        f)
            FASTA_FILE=$OPTARG
            ;;
        c)
            NUMBER_OF_CORES=$OPTARG
            ;;
        e)
            EVALUE=$OPTARG
            ;;
	t)
	    NUMBER_OF_THREADS=$OPTARG
	    ;;
        r)
            NUMBER_OF_REFERENCE_MATCHES=$OPTARG
            ;;
        h)
            echo "$HELPDOC"
            exit 0
            ;;
        \?)
            echo "$HELPDOC"
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
    esac
done

if [ -z $FASTA_FILE ]
then
        echo "$HELPDOC"
        exit 1
fi

fileName=`echo "$(basename $FASTA_FILE)" | cut -d'.' -f1`
if [ $PARALLELIZE_FLAG -eq 1 ]; then

    # Get the file size in KB
    sizeFileBytes=$(du -b $FASTA_FILE | sed 's/\([0-9]*\)\(.*\)/\1/')
    sizeChunks=$(ceil $(float_eval "$sizeFileBytes / ($NUMBER_OF_CORES * 1024)"))
    sizeChunksString="${sizeChunks}k"
    startTime=`date +%s`
    rpsblastOutFmt="\"6 qseqid sseqid evalue pident score qstart qend sstart send length slen\""
	
    cat $FASTA_FILE | parallel --block $sizeChunksString --recstart '>' --pipe rpsblast -outfmt $rpsblastOutFmt -max_target_seqs $NUMBER_OF_REFERENCE_MATCHES -evalue $EVALUE -db $COGSDB_DIR'/Cog' -num_threads $NUMBER_OF_THREADS  -query - > $fileName'.out' 
    echo "rpsblast using GNU parallel took $(expr `date +%s` - $startTime) seconds to generate $fileName.out from $FASTA_FILE"	
else
    startTime=`date +%s`
    rpsblast -outfmt "6 qseqid sseqid evalue pident score qstart qend sstart send length slen" -max_target_seqs $NUMBER_OF_REFERENCE_MATCHES -evalue $EVALUE -query $FASTA_FILE -db $COGSDB_DIR'/Cog' -out $fileName'.out' -num_threads $NUMBER_OF_THREADS
    echo "rpsblast took $(expr `date +%s` - $startTime) seconds to generate $fileName.out from $FASTA_FILE" 
fi
