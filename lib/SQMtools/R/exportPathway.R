#' Export the functions of a SQM object into KEGG pathway maps
#'
#' This function is a wrapper for the pathview package (Luo \emph{et al.}, 2017. \emph{Nucleic acids research}, 45:W501-W508). It will generate annotated KEGG pathway maps showing which reactions are present in the different samples. 
#'
#' @param SQM A SQM or SQMlite object.
#' @param pathway_id character. The five-number KEGG pathway identifier. A list of all pathway identifiers can be found in \url{https://www.genome.jp/kegg/pathway.html}.
#' @param samples character. A vector with the names of the samples to export. If absent, all samples will be exported (default \code{NULL}).
#' @param split_samples logical. Generate a different output file for each sample (default \code{FALSE}).
#' @param output_suffix character. Suffix to be added to the output files (default \code{"pathview"}).
#' @seealso \code{\link[plotFunctions]{plotFunctions}} for plotting the most functions taxa of a SQM object.
#' @examples
#' data(Hadza)
#' exportPathway(Hadza, "00910", output_suffix = "nitrogen_metabolism")
#' @export
exportPathway = function(SQM, pathway_id, samples = NULL, split_samples = F, output_suffix = 'pathview')
    {
    # Check params.
    if(!class(SQM) %in% c('SQM', 'SQMlite')) { stop('The first argument must be a SQM or a SQMlite object') }
    if(!is.null(samples))
        {
	missing_samples = setdiff(samples, SQM$misc$samples)
        if(length(missing_samples) > 0)
            {
            str = paste(missing_samples, collapse = '", "')
	    stop(sprintf('Samples "%s" are not present in this SQM object', str))
	    }
        }
    library(pathview)
    mat = SQM$functions$KEGG$abund
    if(!is.null(samples)) { mat = mat[,samples] }
    mat[mat==0] = NA

    if(split_samples)
        {
        for(sample in colnames(mat))
            {
            sample_suffix = sprintf('%s.%s', output_suffix, sample)
       	    pathview(gene.data = mat[,sample], pathway.id = pathway_id, species = 'ko', out.suffix = sample_suffix, plot.col.key=F, multi.state=F)
	    }
    }else
        {
	pathview(gene.data = mat, pathway.id = pathway_id, species = 'ko', out.suffix = output_suffix, plot.col.key=F, multi.state=T)
        }
    }
