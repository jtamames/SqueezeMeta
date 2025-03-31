#' @importFrom utils globalVariables 
globalVariables(c('RecA', 'MGOGs', 'MGKOs', 'USiCGs'))

#' @importFrom stats median
get_median_single_copy_cov = function(SQM)
    {
    scg = SQM$misc$single_copy_genes

    # Set up source annotations and single copy genes
    if(scg == 'RecA')
        {
        source_annot = 'COG'
        genes = RecA
    } else if(scg == 'MGOGs')
        {
        source_annot = 'COG'
        genes = MGOGs
    } else if(scg == 'MGKOs')
        {
        source_annot = 'KEGG'
        genes = MGKOs
    } else if(scg == 'USiCGs')
        {
        source_annot = 'KEGG'
        genes = USiCGs
    } else
        {
        stop(sprintf('Unknown single copy gene list: %s', scg))
        }

    # Check that we have the source annotations we need
    if(! source_annot %in% names(SQM$functions))
        {
        warning(sprintf('Your project does not contain %s annotations, can not calculate %s copy numbers', source_annot, scg))
        return(NA)
        }

    # For RecA, check whether we have it or not in our dataset
    if(scg == 'RecA' & (! 'COG0468' %in% rownames(SQM$functions$COG$cov) ) )
        {
        warning('RecA is not present in this dataset, can not calculate RecA copy numbers')
        return(NA)
        }
   
    # Calculate median single copy gene coverage
    cov_source = SQM$functions[[source_annot]]$cov
    present_genes = intersect(genes, rownames(cov_source))
    if(!length(present_genes))
        {
        warning('The selected single copy genes were missing from your dataset, can not calculate copy numbers')
        return(NA)
        }
    cov = apply(cov_source[present_genes,,drop=FALSE], 2, median)

    if(any(cov==0))
        {
        warning('Median single copy gene coverage is zero abundance in at least one sample, can not calculate copy numbers')
        return(NA)
        }

    return(cov)
    }

