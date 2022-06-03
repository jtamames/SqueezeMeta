#' Combine several SQM or SQMlite objects
#'
#' Combine an arbitrary number of SQM or SQMlite objects into a single SQMlite object. This function accepts objects originating from different projects (i.e. different SqueezeMeta runs).
#' @param ... an arbitrary number of SQM or SQMlite objects. Alternatively, a single list containing an arbitrary number of SQMlite objects.
#' @return A SQMlite object
#' @seealso \code{\link[subsetFun]{subsetFun}}, \code{\link[subsetTax]{subsetTax}}, \code{\link[combineSQM]{combineSQM}}
#' @examples
#' \dontrun{
#' data(Hadza)
#' # Load data coming from a different run
#' other = loadSQMlite("/path/to/other/project/tables") # e.g. if the project was run using sqm_reads
#' # (We could also use loadSQM to load the data as long as the data comes from a SqueezeMeta run)
#' combined = combineSQMlite(Hadza, other)
#' plotTaxonomy(combined, 'family') # Now we can plot together the samples from Hadza and the second project.
#' }
#' @export
combineSQMlite = function(...)
    {
    inSQMlite = list(...)
    # if there is only one argument and this argument is a list, treat it as a list containing SQM objects
    if(length(inSQMlite) == 1 & class(inSQMlite[[1]]) == 'list') { inSQMlite = list(...)[[1]] }
    return(Reduce(combineSQMlite_, inSQMlite))
    }


combineSQMlite_ = function(SQM1, SQM2)
    {

    if(!class(SQM1) %in% c('SQM', 'SQMlite') | !class(SQM2) %in% c('SQM', 'SQMlite')) { stop('This function only accepts SQM or SQMlite objects') }

    combSQM = list()

    ### Combine taxa
    combSQM$taxa = list()
    for(rank in names(SQM1$taxa))
        {
        combSQM$taxa[[rank]] = list()
	for(count in names(SQM1$taxa[[rank]]))
            {
            combSQM$taxa[[rank]][[count]] = merge.numeric.matrices(SQM1$taxa[[rank]][[count]], SQM2$taxa[[rank]][[count]])
            }
        }

    ### Combine functions
    combSQM$functions = list()
    common_methods = intersect(names(SQM1$functions), names(SQM2$functions))
    common_counts = intersect(names(SQM1$functions[[1]]), names(SQM2$functions[[1]]))
    for(method in common_methods)
        {
        for(count in common_counts)
            {
            combSQM$functions[[method]][[count]] = merge.numeric.matrices(SQM1$functions[[method]][[count]], SQM2$functions[[method]][[count]])
            }
        }

    ### Combine misc info
    combSQM$misc = list()
    combSQM$misc$project_name = c(SQM1$misc$project_name, SQM2$misc$project_name)
    combSQM$misc$samples = c(SQM1$misc$samples, SQM2$misc$samples)
    combSQM$misc$tax_names_long = list()
    for(rank in names(SQM1$misc$tax_names_long))
        {
        combSQM$misc$tax_names_long[[rank]] = named.unique(c(SQM1$misc$tax_names_long[[rank]], SQM2$misc$tax_names_long[[rank]]))
	}
    combSQM$misc$tax_names_short = named.unique(c(SQM1$misc$tax_names_short, SQM2$misc$tax_names_short))
    combSQM$misc$KEGG_names = named.unique(c(SQM1$misc$KEGG_names, SQM2$misc$KEGG_names))
    combSQM$misc$KEGG_paths = named.unique(c(SQM1$misc$KEGG_paths, SQM2$misc$KEGG_paths))
    combSQM$misc$COG_names  = named.unique(c(SQM1$misc$COG_names, SQM2$misc$COG_names))
    combSQM$misc$COG_paths  = named.unique(c(SQM1$misc$COG_paths, SQM2$misc$COG_paths))

    combSQM$misc$ext_annot_sources = intersect(SQM1$misc$ext_annot_sources, SQM2$misc$ext_annot_sources)
    for(method in combSQM$ext_annot_sources)
        {
        fieldn = sprintf('%s_names', method)
        combSQM$misc[[fieldn]] = named.unique(c(SQM1$misc[[fieldn]], SQM2$misc[[fieldn]]))
        }

    ### Combine total reads
    combSQM$total_reads = c(SQM1$total_reads, SQM2$total_reads)

    class(combSQM) = 'SQMlite'
    return(combSQM)
    }

