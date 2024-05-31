#' Print a named vector of sequences as a fasta-formatted string
#'
#' @param seqvec vector. The vector to be written as a fasta string.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @examples
#' data(Hadza)
#' seqvec2fasta(Hadza$orfs$seqs[1:10])
#' @export
seqvec2fasta = function(seqvec, output_name = "")
    {
    cat(sprintf('%s\n', paste(sprintf('>%s\n%s', names(seqvec), seqvec), collapse='\n')), file = output_name)
    }

