#' Export the contigs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportContigs = function(SQM, output_name = "")
    {
    exportSeqs(SQM, 'contigs', output_name)
    }


#' Export the ORFs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name. A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportORFs = function(SQM, output_name = "")
    {
    exportSeqs(SQM, 'orfs', output_name)
    }

exportSeqs = function(SQM, seqtype, output_name = "")
    {
    if(!inherits(SQM, c('SQM', 'SQMbunch'))) { stop('The first argument must be a SQM object') }
    if(inherits(SQM, 'SQM'))
        {
        projs = list(SQM)
        names(projs) = SQM$misc$project_name
    } else
        {
        projs = SQM$projects
        }
    seqvec = c()
    for(projname in names(projs))
        {
	seqs = projs[[projname]][[seqtype]]$seqs
        seqvec = c(seqvec, seqs)        
        if(is.null(seqs))
            {
            warning(sprintf('There are no sequences in project %s. Did you use `load_sequences = FALSE` when loading it?', projname))
            }
        }
    if(length(seqvec)) { seqvec2fasta(seqvec, output_name) } else { stop('No sequences were found') }
    }
 
