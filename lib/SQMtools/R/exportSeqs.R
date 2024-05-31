#' Export the contigs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportContigs = function(SQM, output_name = "")
    {
    seqvec2fasta(SQM$contigs$seqs, output_name)
    }


#' Export the ORFs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportORFs = function(SQM, output_name = "")
    {
    seqvec2fasta(SQM$orfs$seqs, output_name)
    }

