SQMTOOLS_DOCS_IN=../lib/SQMtools/man
SQMTOOLS_DOCS_OUT=source/SQMtools
for file in "$SQMTOOLS_DOCS_IN/"*Rd; do
	# Get the function name
	name=$(basename "$file" .Rd)
	# Get the path for the output rst file
	outfile=$SQMTOOLS_DOCS_OUT/$name.rst
	# Go from *Rd produced by devtools::document to html and then rst
	Rscript -e "tools::Rd2HTML(\"$file\")" | pandoc -f html -o $outfile
	# Add the function name as a title to the output rst file. Format is:
	# *****
	# Title
	# *****
	#
	# We need as many asterisks as characters has the function name
	size=${#name}
	# multiply the "*" character $size times
	heading_string=$(printf '*%.0s' $(seq 1 $size))
	# Use sed to write before the first line, in place
	sed -i "1i $heading_string\n$name\n$heading_string\n" $outfile 
done
