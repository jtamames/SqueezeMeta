#' Export the contigs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportContigs = function(SQM, output_name = "")
    {
    if(!is.null(SQM$contigs$seqs))
        {
        seqvec2fasta(SQM$contigs$seqs, output_name)
    } else
        {
        warning('There are no contig sequences in your SQM project. Did you use `load_sequences = FALSE` when loading it?')
        }
    }


#' Export the ORFs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportORFs = function(SQM, output_name = "")
    {
    if(!is.null(SQM$orfs$seqs))
        {
        seqvec2fasta(SQM$orfs$seqs, output_name)
    } else
        {
        warning('There are no orf sequences in your SQM project. Did you use `load_sequences = FALSE` when loading it?')
        }
    }

