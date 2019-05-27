#!/bin/bash -login
# this bash is used to get any pfam hmm file as $2 and $3 indicate
if [ $# -ne 2 ];then
  echo "get_hmm.sh <pfam_accession> <big hmm file>"
  exit
fi
	awk -v acc="$1" '
								{
									if ($1 ~ "HMMER3") {    # Tocado, JT
										head_line = $0
									} else if ($1 == "NAME") {
										pfam_name_line = $0
									} else if ($1 == "ACC" && $2 ~ acc) {
										flag = 1
										print head_line 
										print pfam_name_line
										print $0
									} else if (flag == 1) {
										print
										if ($1 == "//") {
											flag = 0
											exit 0
										}
									} 
								}
							' $2
										
						
	
