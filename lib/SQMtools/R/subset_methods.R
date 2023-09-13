#' Filter results by sample
#'
#' Create a SQM object containing only samples specified by the user, and the ORFs, contigs, bins, taxa and functions present in those samples.
#' @param SQM SQM object to be subsetted.
#' @param samples character. Samples to be included in the subset.
#' @param remove_missing bool. If \code{TRUE}, ORFs, contigs, bins, taxa and functions absent from the selected samples will be removed from the subsetted object (default \code{TRUE}).
#' @seealso \code{\link{subsetTax}}, \code{\link{subsetFun}}, \code{\link{subsetORFs}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @return SQM object containing only the requested samples.
#' @export
subsetSamples = function(SQM, samples, remove_missing = TRUE)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    check.samples(SQM, samples)
    subSQM = SQM
    subSQM$misc$samples = samples
    subSQM$total_reads = subSQM$total_reads[samples]
    
    ### ORFs
    if(remove_missing) {
        presentORFs = rownames(subSQM$orfs$abund)[ rowSums(subSQM$orfs$abund[,samples,drop=FALSE]) > 0 ]
    } else {
        presentORFs = rownames(subSQM$orfs$abund)
    }
    subSQM$orfs$table       = subSQM$orfs$table[intersect(presentORFs, rownames(subSQM$orfs$table )),       ,drop=FALSE]
    subSQM$orfs$abund       = subSQM$orfs$abund[intersect(presentORFs, rownames(subSQM$orfs$abund )),samples,drop=FALSE]
    subSQM$orfs$bases       = subSQM$orfs$bases[intersect(presentORFs, rownames(subSQM$orfs$bases )),samples,drop=FALSE]
    subSQM$orfs$cov         = subSQM$orfs$cov  [intersect(presentORFs, rownames(subSQM$orfs$cov   )),samples,drop=FALSE]
    subSQM$orfs$cpm         = subSQM$orfs$cpm  [intersect(presentORFs, rownames(subSQM$orfs$cpm   )),samples,drop=FALSE]
    subSQM$orfs$tpm         = subSQM$orfs$tpm  [intersect(presentORFs, rownames(subSQM$orfs$tpm   )),samples,drop=FALSE]
    subSQM$orfs$seqs        = subSQM$orfs$seqs [intersect(presentORFs,    names(subSQM$orfs$seqs  )),        drop=FALSE]
    subSQM$orfs$tax         = subSQM$orfs$tax  [intersect(presentORFs, rownames(subSQM$orfs$tax   )),       ,drop=FALSE]
    if('tax16S' %in% names(subSQM$orfs))
        {
        subSQM$orfs$tax16   = SQM$orfs$tax16S  [intersect(presentORFs, rownames(subSQM$orfs$tax16S)),       ,drop=FALSE]
        }
    if('markers' %in% names(subSQM$orfs))
        {
        subSQM$orfs$markers = SQM$orfs$markers [intersect(presentORFs, names(subSQM$orfs$markers  )),        drop=FALSE]
        }
    ### Contigs
    if(remove_missing) {
        presentContigs = rownames(subSQM$contigs$abund)[ rowSums(subSQM$contigs$abund[,samples,drop=FALSE]) > 0 ]
    } else {
        presentContigs = rownames(subSQM$contigs$abund)
    }
    subSQM$contigs$table = subSQM$contigs$table[intersect(presentContigs, rownames(subSQM$contigs$table)),       ,drop=FALSE]
    subSQM$contigs$abund = subSQM$contigs$abund[intersect(presentContigs, rownames(subSQM$contigs$abund)),samples,drop=FALSE]
    subSQM$contigs$bases = subSQM$contigs$bases[intersect(presentContigs, rownames(subSQM$contigs$bases)),samples,drop=FALSE]
    subSQM$contigs$cov   = subSQM$contigs$cov  [intersect(presentContigs, rownames(subSQM$contigs$cov  )),samples,drop=FALSE]
    subSQM$contigs$cpm   = subSQM$contigs$cpm  [intersect(presentContigs, rownames(subSQM$contigs$cpm  )),samples,drop=FALSE]
    subSQM$contigs$tpm   = subSQM$contigs$tpm  [intersect(presentContigs, rownames(subSQM$contigs$tpm  )),samples,drop=FALSE]
    subSQM$contigs$seqs  = subSQM$contigs$seqs [intersect(presentContigs, rownames(subSQM$contigs$seqs )),        drop=FALSE]
    subSQM$contigs$tax   = subSQM$contigs$tax  [intersect(presentContigs, rownames(subSQM$contigs$tax  )),       ,drop=FALSE]
    subSQM$contigs$bins  = subSQM$contigs$bins [intersect(presentContigs, rownames(subSQM$contigs$bins )),       ,drop=FALSE]
    
    ### Bins
    if('bins' %in% names(subSQM)) {
        if(remove_missing) {
            presentBins = rownames(subSQM$bins$abund)[ rowSums(subSQM$bins$abund[,samples,drop=FALSE]) > 0 ]
        } else {
            presentBins = rownames(subSQM$bins$abund)
        }
        subSQM$bins$table   = subSQM$bins$table  [intersect(presentBins, rownames(subSQM$bins$table  )),       ,drop=FALSE]
        subSQM$bins$length  = subSQM$bins$length [intersect(presentBins, rownames(subSQM$bins$length )),        drop=FALSE]
        subSQM$bins$abund   = subSQM$bins$abund  [intersect(presentBins, rownames(subSQM$bins$abund  )),samples,drop=FALSE]
	subSQM$bins$percent = subSQM$bins$percent[intersect(presentBins, rownames(subSQM$bins$percent)),samples,drop=FALSE]
        subSQM$bins$bases   = subSQM$bins$bases  [intersect(presentBins, rownames(subSQM$bins$bases  )),samples,drop=FALSE]
        subSQM$bins$cov     = subSQM$bins$cov    [intersect(presentBins, rownames(subSQM$bins$cov    )),samples,drop=FALSE]
        subSQM$bins$cpm     = subSQM$bins$cpm    [intersect(presentBins, rownames(subSQM$bins$cpm    )),samples,drop=FALSE]
        subSQM$bins$tax     = subSQM$bins$tax    [intersect(presentBins, rownames(subSQM$bins$tax    )),       ,drop=FALSE]
    }

    ### Taxa
    for(rank in names(subSQM$taxa)) {
        if(remove_missing) {
            presentTaxa = rownames(subSQM$taxa[[rank]]$abund)[ rowSums(subSQM$taxa[[rank]]$abund[,samples,drop=FALSE]) > 0 ]
        } else {
            presentTaxa = rownames(subSQM$taxa[[rank]]$abund)
	}
        for(count in names(subSQM$taxa[[rank]])) {
	    ta = subSQM$taxa[[rank]][[count]]
	    subSQM$taxa[[rank]][[count]] = ta[intersect(presentTaxa, rownames(ta)), samples, drop=FALSE]
	}
    }

    ### Functions
    for(method in names(subSQM$functions)) {
        if(remove_missing) {
            presentFuns = rownames(subSQM$functions[[method]]$abund)[ rowSums(subSQM$functions[[method]]$abund[,samples,drop=FALSE]) > 0 ]
        } else {
            presentFuns = rownames(subSQM$functions[[method]]$abund)
        }
        for(count in names(subSQM$functions[[method]])) {
	    ta = subSQM$functions[[method]][[count]]
            subSQM$functions[[method]][[count]] = ta[intersect(presentFuns, rownames(ta)), samples, drop=FALSE]
        }
    }
    return(subSQM)
}


#' Filter results by function
#'
#' Create a SQM object containing only the ORFs with a given function, and the contigs and bins that contain them.
#' @param SQM SQM object to be subsetted.
#' @param fun character. Pattern to search for in the different functional classifications.
#' @param columns character. Restrict the search to the provided column names from \code{SQM$orfs$table}. If not provided the search will be performed in all the columns containing functional information (default \code{NULL}).
#' @param ignore_case logical Make pattern matching case-insensitive (default \code{TRUE}).
#' @param fixed logical. If \code{TRUE}, pattern is a string to be matched as is. If \code{FALSE} the pattern is treated as a regular expression (default \code{FALSE}).
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin stats and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @seealso \code{\link{subsetTax}}, \code{\link{subsetORFs}}, \code{\link{subsetSamples}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @return SQM object containing only the requested function.
#' @examples
#' data(Hadza)
#' Hadza.iron = subsetFun(Hadza, "iron")
#' Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
#' # Search for multiple patterns using regular expressions
#' Hadza.twoKOs = subsetFun(Hadza, "K00812|K00813", fixed=FALSE)
#' @export
subsetFun = function(SQM, fun, columns = NULL, ignore_case=TRUE, fixed=FALSE, trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE,
		     rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = TRUE)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }

    fun = c(fun) # This suddenly became necessary when testing it in Ubuntu's R 3.6, and now I want to cut myself

    if(is.null(columns))
        { columns = c('Gene name', 'KEGG ID', 'KEGGFUN', 'KEGGPATH', 'COG ID', 'COGFUN', 'COGPATH', 'PFAM')
        for(method in SQM$misc$ext_annot_sources) { columns = c(columns, method, sprintf('%s NAME', method)) }
        }

    goodRows = rep(FALSE, nrow(SQM$orfs$table))
    for(col in columns) { goodRows = goodRows | grepl(fun, SQM$orfs$table[,col], ignore.case = ignore_case, fixed=fixed) }

    goodORFs = rownames(SQM$orfs$table)[goodRows]

    return ( subsetORFs(SQM, goodORFs, tax_source = 'orfs',
                        trusted_functions_only = trusted_functions_only,
                        ignore_unclassified_functions=ignore_unclassified_functions,
                        rescale_tpm         = rescale_tpm,
                        rescale_copy_number = rescale_copy_number,
                        recalculate_bin_stats = recalculate_bin_stats)
           )

    }



#' Filter results by taxonomy
#'
#' Create a SQM object containing only the contigs with a given consensus taxonomy, the ORFs contained in them and the bins that contain them.
#' @param SQM SQM object to be subsetted.
#' @param rank character. The taxonomic rank from which to select the desired taxa (\code{superkingdom}, \code{phylum}, \code{class}, \code{order}, \code{family}, \code{genus}, \code{species})
#' @param tax character. The taxon to select.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object. By default it is set to \code{TRUE}, which means that the returned TPMs will be scaled \emph{by million of reads of the selected taxon}.
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers for each function will represent the average copy number of that function \emph{per genome of the selected taxon}.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin stats and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @return SQM object containing only the requested taxon.
#' @seealso \code{\link{subsetFun}}, \code{\link{subsetContigs}}, \code{\link{subsetSamples}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @examples
#' data(Hadza)
#' Hadza.Prevotella = subsetTax(Hadza, "genus", "Prevotella")
#' Hadza.Proteobacteria = subsetTax(Hadza, "phylum", "Proteobacteria")
#' @export
subsetTax = function(SQM, rank, tax, trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE, recalculate_bin_stats = TRUE)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(!rank %in% colnames(SQM$contigs$tax)) { stop(sprintf('Valid taxonomic ranks are %s', paste(colnames(SQM$contigs$tax), collapse = ', '))) }
    goodContigs = rownames(SQM$contigs$tax)[SQM$contigs$tax[,rank] == tax]
    return ( subsetContigs(SQM, goodContigs,
                           trusted_functions_only = trusted_functions_only,
                           ignore_unclassified_functions=ignore_unclassified_functions,
                           rescale_tpm = rescale_tpm,
                           rescale_copy_number = rescale_copy_number,
	                   recalculate_bin_stats = recalculate_bin_stats)
           )
   
    }


#n Select bins
#'
#' Create a SQM object containing only the requested bins, and the contigs and ORFs contained in them.
#' @param SQM SQM object to be subsetted.
#' @param bins character. Vector of bins to be selected.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object. By default it is set to \code{TRUE}, which means that the returned TPMs will be scaled \emph{by million of reads of the selected bins}.
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers for each function will represent the average copy number of that function \emph{per genome of the selected bins}.
#' @return SQM object containing only the requested bins.
#' @seealso \code{\link{subsetContigs}}, \code{\link{subsetORFs}}
#' @examples 
#' data(Hadza)
#' # Which are the most complete bins?
#' topBinNames = rownames(Hadza$bins$table)[order(Hadza$bins$table[,"Completeness"],
#'                                          decreasing=TRUE)][1:2]
#' # Subset with the most complete bin.
#' topBin = subsetBins(Hadza, topBinNames[1])
#' @export
subsetBins = function(SQM, bins, trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    goodContigs = rownames(SQM$contigs$bins)[SQM$contigs$bins %in% bins]
    return ( subsetContigs(SQM, goodContigs,  
                           trusted_functions_only = trusted_functions_only,
                           ignore_unclassified_functions=ignore_unclassified_functions,
                           rescale_tpm = rescale_tpm,
                           rescale_copy_number = rescale_copy_number,
	                   recalculate_bin_stats = FALSE)
           )
    }


#' Select contigs
#'
#' Create a SQM object containing only the requested contigs, the ORFs contained in them and the bins that contain them.
#' @param SQM SQM object to be subsetted.
#' @param contigs character. Vector of contigs to be selected.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin stats and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @return SQM object containing only the selected contigs.
#' @seealso \code{\link{subsetORFs}}
#' @examples
#' data(Hadza)
#' # Which contigs have a GC content below 40?
#' lowGCcontigNames = rownames(Hadza$contigs$table[Hadza$contigs$table[,"GC perc"]<40,])
#' lowGCcontigs = subsetContigs(Hadza, lowGCcontigNames)
#' hist(lowGCcontigs$contigs$table[,"GC perc"])
#' @export
subsetContigs = function(SQM, contigs, trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = TRUE)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    goodORFs = rownames(SQM$orfs$table)[SQM$orfs$table[,'Contig ID'] %in% contigs]
    return ( subsetORFs(SQM, goodORFs, tax_source = 'contigs',
                        trusted_functions_only = trusted_functions_only,
                        ignore_unclassified_functions=ignore_unclassified_functions,
                        rescale_tpm = rescale_tpm,
                        rescale_copy_number = rescale_copy_number,
			recalculate_bin_stats = recalculate_bin_stats,
	                contigs_override = contigs)
           )
    }


#' Select random ORFs
#'
#' Create a random subset of a SQM object.
#' @param SQM SQM object to be subsetted.
#' @param N numeric. number of random ORFs to select.
#' @return SQM object containing a random subset of ORFs.
#' @seealso \code{\link{subsetORFs}}
#' @export
subsetRand = function(SQM, N)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    goodORFs = sample(rownames(SQM$orfs$table), N)
    return ( subsetORFs(SQM, goodORFs) )
    }


#' Select ORFs
#'
#' Create a SQM object containing only the requested ORFs, and the contigs and bins that contain them. Internally, all the other subset functions in this package end up calling \code{subsetORFs} to do the work for them.
#' @param SQM SQM object to be subsetted.
#' @param orfs character. Vector of ORFs to be selected.
#' @param tax_source character. Features used for calculating aggregated abundances at the different taxonomic ranks. Either \code{"orfs"} or \code{"contigs"} (default \code{"orfs"}).
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin stats and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @param contigs_override character. Optional vector of contigs to be included in the subsetted object.
#' @return SQM object containing the requested ORFs.
#' @section A note on contig/bins subsetting:
#' While this function selects the contigs and bins that contain the desired orfs, it DOES NOT recalculate contig/bin abundance and statistics based on the selected ORFs only. This means that the abundances presented in tables such as \code{SQM$contig$abund} or \code{SQM$bins$tpm} will still refer to the complete contigs and bins, regardless of whether only a fraction of their ORFs are actually present in the returned SQM object. This is also true for the statistics presented in \code{SQM$contigs$table} and \code{SQM$bins$table}.
#' @examples
#' data(Hadza)
#' # Select the 100 most abundant ORFs in our dataset.
#' mostAbundantORFnames = names(sort(rowSums(Hadza$orfs$tpm), decreasing=TRUE))[1:100]
#' mostAbundantORFs = subsetORFs(Hadza, mostAbundantORFnames)
#' @importFrom stats aggregate
#' @export
subsetORFs = function(SQM, orfs, tax_source = 'orfs', trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, 
		      rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = TRUE, contigs_override = NULL)
    {

    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(length(orfs)==0) { stop('No ORFs were selected. Perhaps the subset query yielded no results?') }
    if(!tax_source %in% c('contigs', 'orfs')) { stop('tax_source must be "orfs" or "contigs"') }
   
    orfs    = rownames(SQM$orfs$table[orfs,,drop=FALSE]) # Make sure it will work if orfs is a bool vector too.
    if(is.null(contigs_override)) { contigs = unique(SQM$orfs$table[orfs,'Contig ID'])
    } else { contigs = contigs_override } # so we can include contigs without ORFs if required
    bins    = unique( unlist(SQM$contigs$bins[contigs,]) )
    bins    = bins[bins!='No_bin']
    
    subSQM = SQM
    ### ORFs
    #    Table
    subSQM$orfs$table                 = SQM$orfs$table  [orfs    ,,drop=FALSE]
    #    Abundances
    subSQM$orfs$abund                 = SQM$orfs$abund  [orfs    ,,drop=FALSE]

    subSQM$orfs$bases                 = SQM$orfs$bases  [orfs    ,,drop=FALSE]
    subSQM$orfs$cov                   = SQM$orfs$cov    [orfs    ,,drop=FALSE]
    subSQM$orfs$cpm                   = SQM$orfs$cpm    [orfs    ,,drop=FALSE]
    subSQM$orfs$tpm                   = SQM$orfs$tpm    [orfs    ,,drop=FALSE]
    #    Sequences
    subSQM$orfs$seqs                  = SQM$orfs$seqs   [orfs]
    subSQM$orfs$tax                   = SQM$orfs$tax    [orfs    ,,drop=FALSE]
    #    Taxonomy
    if('tax16S' %in% names(subSQM$orfs))
        {
        subSQM$orfs$tax16S             = SQM$orfs$tax16S[orfs                ]
        }
    #    Markers
    if('markers' %in% names(subSQM$orfs))
        {
        subSQM$orfs$markers           = SQM$orfs$markers [orfs]
        }
    ### Contigs
    #    Table
    subSQM$contigs$table              = SQM$contigs$table[contigs,,drop=FALSE]
    #    Abundances
    subSQM$contigs$abund              = SQM$contigs$abund[contigs,,drop=FALSE]
    subSQM$contigs$bases              = SQM$contigs$bases[contigs,,drop=FALSE]
    subSQM$contigs$cov                = SQM$contigs$cov  [contigs,,drop=FALSE]
    subSQM$contigs$cpm                = SQM$contigs$cpm  [contigs,,drop=FALSE]
    subSQM$contigs$tpm                = SQM$contigs$tpm  [contigs,,drop=FALSE]
    #    Sequences
    subSQM$contigs$seqs               = SQM$contigs$seqs [contigs]
    #    Taxonomy
    subSQM$contigs$tax                = SQM$contigs$tax  [contigs,,drop=FALSE]
    #    Binning info
    if('bins' %in% names(subSQM))
        {
	subSQM$contigs$bins           = SQM$contigs$bins[contigs ,,drop=FALSE]
	### Bins
	if(recalculate_bin_stats & (!('tax16S' %in% names(subSQM$orfs)) & 'markers' %in% names(subSQM$orfs)))
            {
            warning('You requested to recalculate bin stats but 16S or marker gene info are missing. Will not recalculate bin stats')
            }
	# Table and Taxonomy
        if(recalculate_bin_stats & ('tax16S' %in% names(subSQM$orfs) & 'markers' %in% names(subSQM$orfs)))
            {
            bin_stats                 = get.bin.stats(subSQM)
            subSQM$bins$table         = bin_stats[['table']]
            subSQM$bins$tax           = bin_stats[['tax']] 
        }else
            {
            subSQM$contigs$bins       = SQM$contigs$bins[contigs ,,drop=FALSE]
            subSQM$bins$tax           = subSQM$bins$tax[bins     ,,drop=FALSE]
            }
	#    Abundances
	bin_abunds                    = get.bin.abunds(subSQM)
	subSQM$bins$abund             = bin_abunds[['abund']]
	subSQM$bins$percent           = bin_abunds[['percent']]
	subSQM$bins$bases             = bin_abunds[['bases']]
	subSQM$bins$length            = bin_abunds[['length']]
	subSQM$bins$cov               = bin_abunds[['cov']]
	subSQM$bins$cpm               = bin_abunds[['cpm']]
        }

    ### Taxonomy
    subSQM$taxa$superkingdom$abund    = aggregate_taxa(subSQM, 'superkingdom', tax_source)
    subSQM$taxa$phylum$abund          = aggregate_taxa(subSQM, 'phylum'      , tax_source)
    subSQM$taxa$class$abund           = aggregate_taxa(subSQM, 'class'       , tax_source)
    subSQM$taxa$order$abund           = aggregate_taxa(subSQM, 'order'       , tax_source)
    subSQM$taxa$family$abund          = aggregate_taxa(subSQM, 'family'      , tax_source)
    subSQM$taxa$genus$abund           = aggregate_taxa(subSQM, 'genus'       , tax_source)
    subSQM$taxa$species$abund         = aggregate_taxa(subSQM, 'species'     , tax_source)

    subSQM$taxa$superkingdom$percent  = 100 * t(t(subSQM$taxa$superkingdom$abund) / subSQM$total_reads)
    subSQM$taxa$phylum$percent        = 100 * t(t(subSQM$taxa$phylum$abund)       / subSQM$total_reads)
    subSQM$taxa$class$percent         = 100 * t(t(subSQM$taxa$class$abund)        / subSQM$total_reads)
    subSQM$taxa$order$percent         = 100 * t(t(subSQM$taxa$order$abund)        / subSQM$total_reads)
    subSQM$taxa$family$percent        = 100 * t(t(subSQM$taxa$family$abund)       / subSQM$total_reads)
    subSQM$taxa$genus$percent         = 100 * t(t(subSQM$taxa$genus$abund)        / subSQM$total_reads)
    subSQM$taxa$species$percent       = 100 * t(t(subSQM$taxa$species$abund)      / subSQM$total_reads)

    ### Functions
    if('KEGG' %in% names(subSQM$functions))
        {
        KEGG                          = aggregate_fun(subSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$KEGG$abund   = KEGG$abund
        subSQM$functions$KEGG$bases   = KEGG$bases
        subSQM$functions$KEGG$cov     = KEGG$cov
	subSQM$functions$KEGG$cpm     = t(t(KEGG$cov) /  (subSQM$total_reads/1000000))
        }
    
    if('COG' %in% names(subSQM$functions))
        {
        COG                           = aggregate_fun(subSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$COG$abund    = COG$abund
        subSQM$functions$COG$bases    = COG$bases
        subSQM$functions$COG$cov      = COG$cov
	subSQM$functions$COG$cpm      = t(t(COG$cov) /  (subSQM$total_reads/1000000))
        }


    if('PFAM' %in% names(subSQM$functions))
        {
        PFAM                          = aggregate_fun(subSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$PFAM$abund   = PFAM$abund
        subSQM$functions$PFAM$bases   = PFAM$bases
        subSQM$functions$PFAM$cov     = PFAM$cov
	subSQM$functions$PFAM$cpm     = t(t(PFAM$cov) /  (subSQM$total_reads/1000000))
        }

    ext_annots = list()
    for(method in subSQM$misc$ext_annot_sources)
        {
        ext_annots[[method]]          = aggregate_fun(subSQM, method, trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions[[method]]$abund = ext_annots[[method]]$abund
        subSQM$functions[[method]]$bases = ext_annots[[method]]$bases
        subSQM$functions[[method]]$cov   = ext_annots[[method]]$cov
	subSQM$functions[[method]]$cpm   = t(t(subSQM$functions[[method]]$cov) /  (subSQM$total_reads/1000000))
        }

    if(rescale_tpm)
        {
        if('KEGG' %in% names(subSQM$functions)) { subSQM$functions$KEGG$tpm = KEGG$tpm_rescaled }
        if('COG'  %in% names(subSQM$functions)) { subSQM$functions$COG$tpm  = COG$tpm_rescaled  }
        if('PFAM' %in% names(subSQM$functions)) { subSQM$functions$PFAM$tpm = PFAM$tpm_rescaled }

        for(method in subSQM$misc$ext_annot_sources)
            {
            subSQM$functions[[method]]$tpm = ext_annots[[method]]$tpm_rescaled
            }
        subSQM$orfs$tpm               = 1000000 * t(t(subSQM$orfs$tpm)   /colSums(subSQM$orfs$tpm)   )
        subSQM$contigs$tpm            = 1000000 * t(t(subSQM$contigs$tpm)/colSums(subSQM$contigs$tpm))

	for(method in names(subSQM$functions))
            {
            subSQM$misc$coding_fraction[[method]]        = rep(1, ncol(subSQM$orfs$tpm))
	    names(subSQM$misc$coding_fraction[[method]]) = names(subSQM$orfs$tpm)
            }

    }else
        {
        if('KEGG' %in% names(subSQM$functions)) { subSQM$functions$KEGG$tpm = KEGG$tpm }
        if('COG'  %in% names(subSQM$functions)) { subSQM$functions$COG$tpm  = COG$tpm  }
        if('PFAM' %in% names(subSQM$functions)) { subSQM$functions$PFAM$tpm = PFAM$tpm }
        for(method in subSQM$misc$ext_annot_sources)
            {
            subSQM$functions[[method]]$tpm = ext_annots[[method]]$tpm
            }
    }

    if(!is.null(subSQM$misc$RecA_cov))
        {
        if(rescale_copy_number)
            {
            if('COG0468' %in% rownames(COG$cov))
                {
                if(all(COG$cov['COG0468',]>0))
	            {
                    RecA = COG$cov['COG0468',]
		}else
	            {
	            warning('RecA has zero abundance in at least one sample in this subset. Will not rescale copy numbers.')
		    RecA = SQM$misc$RecA_cov
		    }
            }else
                {
                warning('RecA is not present in this subset. Will not rescale copy numbers.')
                RecA = SQM$misc$RecA_cov
                } 
        }else
           {
           RecA = SQM$misc$RecA_cov
           }
        if('KEGG' %in% names(subSQM$functions)) { subSQM$functions$KEGG$copy_number = t(t(KEGG$cov) / RecA) }
        if('COG'  %in% names(subSQM$functions)) { subSQM$functions$COG$copy_number  = t(t(COG$cov ) / RecA) }
        if('PFAM' %in% names(subSQM$functions)) { subSQM$functions$PFAM$copy_number = t(t(PFAM$cov) / RecA) }
        for(method in subSQM$misc$ext_annot_sources)
            {
            subSQM$functions[[method]]$copy_number = t(t(ext_annots[[method]]$cov) / RecA)
            }
        subSQM$misc$RecA_cov              = RecA
        }


    #subSQM$total_reads		     = colSums(subSQM$contigs$abund)

    ### Total reads
    return(subSQM)
    }

