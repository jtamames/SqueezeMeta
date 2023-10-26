#' Print a named vector of sequences as a fasta-formatted string
#'
#' @param seqvec vector. The vector to be written as a fasta string.
#' @return a character string with the sequence/s written in the fasta format.
#' @examples
#' data(Hadza)
#' seqvec2fasta(Hadza$orfs$seqs[1:10])
#' @export
seqvec2fasta = function(seqvec)
    {
    cat(sprintf('%s\n', paste(sprintf('>%s\n%s', names(seqvec), seqvec), collapse='\n')))
    }

