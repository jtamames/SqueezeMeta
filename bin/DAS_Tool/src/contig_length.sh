#!/usr/bin/env bash

contigs=$1
length_table=$2

awk '/^>/ {if (x){print x}; print ;x=0;next; } {x += length($0)}END{print x}' $contigs | tr -d \\r | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n\([^>]\)/?\1/g' | tr '?' '\t' | sed 's/>//g' > $length_table
