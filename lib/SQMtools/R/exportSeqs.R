#' Export the contigs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportContigs = function(SQM, output_name = "")
    {
    exportSeqs(SQM, 'contigs', output_name)
    }


#' Export the ORFs of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_name A connection, or a character string naming the file to print to. If "" (the default), sequences will be printed to the standard output connection.
#' @export
exportORFs = function(SQM, output_name = "")
    {
    exportSeqs(SQM, 'orfs', output_name)
    }


#' Export the bins of a SQM object
#'
#' @param SQM A SQM object.
#' @param output_dir Existing output directory to which to write the bins.
#' @export
exportBins = function(SQM, output_dir = "")
    {
    if(!inherits(SQM, c('SQM', 'SQMbunch'))) { stop('The first argument must be a SQM or SQMbunch object') }
    if(inherits(SQM, 'SQM'))
        {
        projs = list(SQM)
        names(projs) = SQM$misc$project_name
    } else
        {
        projs = SQM$projects
        }
    seen_bins = c()
    for(projname in names(projs))
        {
        proj = projs[[projname]]
        if(is.null(proj$bins))
            {
            warning(sprintf('Project %s does not contain bins', projname))
            next
            }
        bins = unique(proj$contigs$bins[,1])
        bins = bins[bins != 'No_bin']
        for(b in bins)
            {
            contigs  = rownames(proj$contigs$bins)[proj$contigs$bins[,1]==b]
            if(b %in% seen_bins)
                {
                warning(sprintf('Bin %s is present in project %s but a bin with a similar name was found in another project in this bunch. Skipping', b, projname))
                next
                }
            seqvec2fasta(proj$contigs$seqs[contigs], sprintf('%s/%s.fasta', output_dir, b))
            seen_bins = c(seen_bins, b)
            }
        }
    }

#' @importFrom Biostrings writeXStringSet
exportSeqs = function(SQM, seqtype, output_name = "")
    {
    if(!inherits(SQM, c('SQM', 'SQMbunch'))) { stop('The first argument must be a SQM or SQMbunch object') }
    if(inherits(SQM, 'SQM'))
        {
        projs = list(SQM)
        names(projs) = SQM$misc$project_name
    } else
        {
        projs = SQM$projects
        }
    seqvec = NA
    for(projname in names(projs))
        {
	seqs = projs[[projname]][[seqtype]]$seqs
        duplic = intersect(names(seqs), names(seqvec))
        if(length(duplic))
            {
            warning(sprintf('Project %s contains %s sequence names that were also present in other projects in this object', projname, length(duplic)))
            }
        if(is.null(seqs))
            {
            warning(sprintf('There are no sequences in project %s. Did you use `load_sequences = FALSE` when loading it?', projname))
            }
        # concatenating an empty vector and a XStringSet would return a list instead of a larger XStringSet, so we do this
        # seqvec starts being NA so on the first project we read we overwrite with the sequences, THEN we concatenate
        if(is.na(seqvec[1])) { seqvec = seqs } else { seqvec = c(seqvec, seqs) }
        }
    if(length(seqvec)) { writeXStringSet(seqvec, output_name) } else { stop('No sequences were found') }
    }
 
