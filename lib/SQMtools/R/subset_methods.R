subsetDispatch = function(f, SQM, ...)
    {
    if(!inherits(SQM, c('SQM', 'SQMbunch'))) { stop('The first argument must be a SQM or SQMbunch object') }
    if(!inherits(SQM, 'SQMbunch'))
        {
        args = c(list(SQM=SQM), list(...))
        subSQM = do.call(f, args)
    } else {
        projs = list()
        for(p in SQM$projects)
            {
            args = c(list(SQM=p), list(...))
	    args$allow_empty = TRUE
            projs = c(projs, list(do.call(f, args)))
            }
        subSQM = combineSQM(projs)
        }
    return(subSQM) 
    }


fix_samples_tax_fun = function(la, samples, new_sample_names, do_subset, remove_missing, do_rename)
    {
    res = la
    for(method in names(res))
        {
        for(count in names(res[[method]]))
            {
            ta = res[[method]][[count]]
            features = NULL
            if(remove_missing)
                {
                features = rownames(ta)[ rowSums(ta[,samples,drop=FALSE]) > 0 ]
                }
            res[[method]][[count]] = fix_samples(res[[method]][[count]], features, samples, new_sample_names,
                                                 do_subset, remove_missing, do_rename)
            }
        }
    return(res)
    }

fix_samples = function(ta, features, samples, new_sample_names, do_subset, remove_missing, do_rename)
    {
    if(remove_missing) { stopifnot(!is.null(features)) }
    if(do_rename) { stopifnot(length(samples) == length(new_sample_names)) }
    res = ta
    if('matrix' %in% class(res) | 'data.frame' %in% class(res)) # this is a matrix or data.frame-like
        {
        has_samples = all(samples %in% colnames(res))
        if(do_subset)
            {
            if(remove_missing)
                {
                features = intersect(features, rownames(res))
                res = res[features,,drop=FALSE] # We need the intersect for bins, as the abundance table
                }                               #  (used to get the features) contains No bin and Unmapped
            if(has_samples) { res = res[,samples,drop=FALSE] }
            }
        if(do_rename & all(samples %in% colnames(res)))
            {
            colnames(res) = new_sample_names
            }
    } else # assume this is a vector or vector-like
        {
        if(do_subset)
            {
            if(remove_missing)
                {
                features = intersect(features, names(res))
                res = res[features] # We need the intersect for orfs$tax16S, as it doesn't contain all the ORFs
                }
            }
        }
    return(res)
    }


#' Change sample names
#'
#' Change the sample names of a SQM or SQMlite object
#'
#' @param SQM SQM or SQMlite object
#' @param new_sample_names character. New sample names
#' @return SQM or SQMlite object with the new sample names
#' @export
renameSamples = function(SQM, new_sample_names)
    {
    return(subsetSamples(SQM, SQM$misc$samples, remove_missing = FALSE, new_sample_names = new_sample_names))
    }


#' Filter results by sample
#'
#' Create a SQM or SQMlite object containing only samples specified by the user, and the ORFs, contigs, bins, taxa and functions present in those samples.
#' @param SQM SQM or SQMlite object to be subsetted.
#' @param samples character. Samples to be included in the subset.
#' @param remove_missing bool. If \code{TRUE}, ORFs, contigs, bins, taxa and functions absent from the selected samples will be removed from the subsetted object (default \code{TRUE}).
#' @param new_sample_names character. New sample names to be included in the subset (default \code{NULL}).
#' @seealso \code{\link{subsetTax}}, \code{\link{subsetFun}}, \code{\link{subsetORFs}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @return SQM or SQMlite object containing only the requested samples.
#' @export
subsetSamples = function(SQM, samples, remove_missing = TRUE, new_sample_names = NULL)
    {
    if(!inherits(SQM, 'SQM') & !inherits(SQM, 'SQMlite'))
        {
        stop('The first argument must be a SQM or SQMlite object')
        }
    do_subset = TRUE
    do_rename = FALSE
    if(!is.null(new_sample_names)) { do_rename = TRUE }
    if(!is.null(new_sample_names))
        {
        if(length(new_sample_names) != length(samples)) # rename only
            {
            do_subset = FALSE
            remove_missing = FALSE
            if(do_rename)
                {
                stop('`samples` and `new_sample_names` must have the same length')
                }
            }
        }
    check.samples(SQM, samples)
    if(identical(SQM$misc$sample_names, samples)) { do_subset = FALSE } # we are just changing sample names
    if(!do_rename & !do_subset) { return(SQM) } # do nothing

    subSQM = SQM
    subSQM$misc$samples = samples
    subSQM$total_reads = subSQM$total_reads[samples]
    if(do_rename)
        {
        subSQM$misc$samples = new_sample_names
        if(!is.null(subSQM$total_reads)) { names(subSQM$total_reads) = samples }
        }
        
    ### orfs, contigs, bins
    for(feat in c('orfs', 'contigs', 'bins'))
        {
        if(!is.null(SQM[[feat]]))
            {
            features = NULL
            if(remove_missing)
                {
                features = rownames(subSQM[[feat]]$abund)[ rowSums(subSQM[[feat]]$abund[,samples,drop=FALSE]) > 0 ]
                }
            for(tn in c('table', 'abund', 'bases', 'cov', 'cpm', 'tpm', 'seqs', 'tax', 'tax16S', 'markers',
                        'bins', 'length', 'percent', 'tax_gtdb'))
                {
                if(!is.null(subSQM[[feat]][[tn]]))
                    {
                    subSQM[[feat]][[tn]] = fix_samples(subSQM[[feat]][[tn]], features, samples, new_sample_names,
                                                       do_subset, remove_missing, do_rename)
                    }
                }
            for(tabund in c('tax_abund', 'tax_abund_gtdb'))
                {
                if(!is.null(subSQM[[feat]][[tabund]]))
                    {
                    subSQM[[feat]][[tabund]] = fix_samples_tax_fun(subSQM[[feat]][[tabund]], samples, new_sample_names,
                                                                   do_subset, remove_missing, do_rename)
                    }
                }
                
            }
        }

    ### Taxa and functions
    for(feat in c('taxa', 'functions'))
        {
        if(!is.null(subSQM[[feat]]))
            {
            subSQM[[feat]] = fix_samples_tax_fun(subSQM[[feat]], samples, new_sample_names,
                                                 do_subset, remove_missing, do_rename)
            }
        }
    
    return(subSQM)
    }


#' Filter results by function
#'
#' Create a SQM or SQMbunch object containing only the ORFs with a given function, and the contigs and bins that contain them.
#' @param SQM SQM or SQMbunch object to be subsetted.
#' @param fun character. Pattern to search for in the different functional classifications.
#' @param columns character. Restrict the search to the provided column names from \code{SQM$orfs$table}. If not provided the search will be performed in all the columns containing functional information (default \code{NULL}).
#' @param ignore_case logical Make pattern matching case-insensitive (default \code{TRUE}).
#' @param fixed logical. If \code{TRUE}, pattern is a string to be matched as is. If \code{FALSE} the pattern is treated as a regular expression (default \code{FALSE}).
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{FALSE}).
#' @param allow_empty (internal use only).
#' @seealso \code{\link{subsetTax}}, \code{\link{subsetORFs}}, \code{\link{subsetSamples}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @return SQM or SQMbunch object containing only the requested function.
#' @examples
#' data(Hadza)
#' Hadza.iron = subsetFun(Hadza, "iron")
#' Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
#' # Search for multiple patterns using regular expressions
#' Hadza.twoKOs = subsetFun(Hadza, "K00812|K00813", fixed=FALSE)
#' @export
subsetFun = function(SQM, fun, columns = NULL, ignore_case=TRUE, fixed=FALSE,
                     trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE,
                     rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = FALSE, allow_empty = FALSE)
    {
    return(subsetDispatch(subsetFun_, SQM, fun, columns=columns, ignore_case = ignore_case, fixed = fixed,
                     trusted_functions_only = trusted_functions_only, ignore_unclassified_functions = ignore_unclassified_functions,
                     rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number, recalculate_bin_stats = recalculate_bin_stats,
		     allow_empty = allow_empty)
          )
    }
subsetFun_ = function(SQM, fun, columns, ignore_case, fixed,
                      trusted_functions_only, ignore_unclassified_functions,
		      rescale_tpm, rescale_copy_number, recalculate_bin_stats, allow_empty)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(SQM$misc$onlybins) { stop('This function can not be run on projects generated with the `--onlybins` flag') }

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
                        recalculate_bin_stats = recalculate_bin_stats,
		        allow_empty = allow_empty)
           )

    }



#' Filter results by taxonomy
#'
#' Create a SQM or SQMbunch object containing only the contigs/bins with a given consensus taxonomy, as well as the ORFs contained in them.
#' @param SQM SQM object to be subsetted.
#' @param rank character. The taxonomic rank from which to select the desired taxa (\code{superkingdom}, \code{phylum}, \code{class}, \code{order}, \code{family}, \code{genus}, \code{species})
#' @param tax character. A taxon or vector of taxa to be selected.
#' @param tax_source character, source data used for feature selection, and to generate the taxonomy tables present in \code{SQM$taxa}, either \code{"orfs"}, \code{"contigs"}, \code{"bins"} (GTDB bin taxonomy if available, SQM bin taxonomy otherwise), \code{"bins_gtdb"} (GTDB bin taxonomy) or \code{"bins_sqm"} (SQM bin taxonomy). When \code{"bins"}, \code{"bins_gtdb"} or \code{"bins_sqm"}, this function will select the bins from the desired taxa, otherwise it will select the contigs from the desired taxa. If using \code{"bins_gtdb"}, note that GTDB taxonomy may differ from the NCBI taxonomy used throughout the rest of SqueezeMeta. Default \code{"contigs"}, unless the project was created with the `--onlybins` flag, where it will be \code{"bins_gtdb"} if GTDB taxonomy is available for the bins.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object. By default it is set to \code{TRUE}, which means that the returned TPMs will be scaled \emph{by million of reads of the selected taxon}.
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers for each function will represent the average copy number of that function \emph{per genome of the selected taxon}.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @param allow_empty (internal use only).
#' @return SQM or SQMbunch object containing only the requested taxon.
#' @seealso \code{\link{subsetFun}}, \code{\link{subsetContigs}}, \code{\link{subsetSamples}}, \code{\link{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be selected with \code{\link{mostAbundant}}.
#' @examples
#' data(Hadza)
#' Hadza.Prevotella = subsetTax(Hadza, "genus", "Prevotella")
#' Hadza.Bacteroidota = subsetTax(Hadza, "phylum", "Bacteroidota")
#' @export
subsetTax = function(SQM, rank, tax, tax_source = NULL,
                     trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE,
                     recalculate_bin_stats = FALSE, allow_empty = FALSE)
    {
    return(subsetDispatch(subsetTax_, SQM, rank, tax, tax_source = tax_source,
			  trusted_functions_only = trusted_functions_only, ignore_unclassified_functions = ignore_unclassified_functions,
			  rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number, recalculate_bin_stats = recalculate_bin_stats,
			  allow_empty = allow_empty)
          )
    }
subsetTax_ = function(SQM, rank, tax, tax_source, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number, recalculate_bin_stats, allow_empty)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(is.null(tax_source)) { if(SQM$misc$onlybins) { tax_source = 'bins_gtdb' } else { tax_source = 'contigs' } }
    if(tax_source %in% c('bins', 'bins_gtdb', 'bins_sqm'))
        {
        subSQM = subsetBins(SQM, bins = NULL, rank = rank, tax = tax, min_completeness = NULL, max_contamination = NULL, tax_source = tax_source,
                            trusted_functions_only = trusted_functions_only,
                            ignore_unclassified_functions=ignore_unclassified_functions,
                            rescale_tpm = rescale_tpm,
                            rescale_copy_number = rescale_copy_number,
                            allow_empty = allow_empty)
    } else
        {
        if(!rank %in% colnames(SQM$contigs$tax)) { stop(sprintf('Valid taxonomic ranks are %s', paste(colnames(SQM$contigs$tax), collapse = ', '))) }
        goodContigs = rownames(SQM$contigs$tax)[SQM$contigs$tax[,rank] %in% tax]
        subSQM = subsetContigs(SQM, goodContigs, tax_source = tax_source,
                               trusted_functions_only = trusted_functions_only,
                               ignore_unclassified_functions=ignore_unclassified_functions,
                               rescale_tpm = rescale_tpm,
                               rescale_copy_number = rescale_copy_number,
	                       recalculate_bin_stats = recalculate_bin_stats,
	                       allow_empty = allow_empty)
        }
    return(subSQM)
    }


#n Select bins
#'
#' Create a SQM object containing only the requested bins, and the contigs and ORFs contained in them.
#' @param SQM SQM object to be subsetted.
#' @param bins character. Vector of bins to be selected. If provided, will override \code{rank}, \code{tax}, \code{min_completeness} and \code{max_contamination}.
#' @param rank character. The taxonomic rank from which to select the desired taxa (\code{superkingdom}, \code{phylum}, \code{class}, \code{order}, \code{family}, \code{genus}, \code{species})
#' @param tax character. A taxon or vector of taxa to be selected.
#' @param min_completeness numeric. Discard bins with completeness lower than this value (default \code{NULL}).
#' @param max_contamination numeric. Discard bins with contamination higher than this value (default \code{NULL}).
#' @param tax_source character, source data used for taxonomic subsetting (if \code{rank} and \code{tax} are provided) and for the aggregate taxonomy tables present in \code{SQM$taxa}, either \code{"orfs"}, \code{"contigs"}, \code{"bins"} (GTDB bin taxonomy if available, SQM bin taxonomy otherwise), \code{"bins_gtdb"} (GTDB bin taxonomy) or \code{"bins_sqm"} (SQM bin taxonomy). If using \code{bins_gtdb}, note that GTDB taxonomy may differ from the NCBI taxonomy used throughout the rest of SqueezeMeta. Default \code{"bins"}.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object. By default it is set to \code{TRUE}, which means that the returned TPMs will be scaled \emph{by million of reads of the selected bins}.
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers for each function will represent the average copy number of that function \emph{per genome of the selected taxon}.
#' @param allow_empty (internal use only).
#' @return SQM object containing only the requested bins.
#' @seealso \code{\link{subsetContigs}}, \code{\link{subsetORFs}}
#' @examples 
#' data(Hadza)
#' # Which are the most complete bins?
#' topBinNames = rownames(Hadza$bins$table)[order(Hadza$bins$table[,"Completeness"],
#'                                          decreasing=TRUE)][1:2]
#' # Subset with the most complete bin.
#' topBin = subsetBins(Hadza, topBinNames[1])
#' 
#' # Subset with all the bins over 90% completeness
#' over90 = subsetBins(Hadza, min_completeness = 90)
#'
#' # Subset with bins from the Phascolarctobacterium genus using SqueezeMeta's taxonomy
#' phasco = subsetBins(Hadza, tax_source = "bins", rank = "genus", tax = "Phascolarctobacterium")
#'
#' # Subset with binsfrom the Bacteroidota phylum using GTDB taxonomy
#' bact = subsetBins(Hadza, tax_source = "bins_gtdb", rank = "phylum", tax = "p__Bacteroidota")
#'
#' @export
subsetBins = function(SQM, bins = NULL, rank = NULL, tax = NULL, min_completeness = NULL, max_contamination = NULL, tax_source = 'bins',
                      trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE, allow_empty = FALSE)
    {
    return(subsetDispatch(subsetBins_, SQM, bins, rank, tax, min_completeness, max_contamination, tax_source = tax_source,
			  trusted_functions_only = trusted_functions_only,
                          ignore_unclassified_functions = ignore_unclassified_functions,
			  rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number, allow_empty = allow_empty)
          )
    }
subsetBins_ = function(SQM, bins, rank, tax, min_completeness, max_contamination, tax_source,
                       trusted_functions_only, ignore_unclassified_functions,
                       rescale_tpm, rescale_copy_number, allow_empty)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    
    if(!is.null(bins))
        {
        goodBins = bins
    } else
        {
        goodBins = rownames(SQM$bins$table)
        if(!is.null(rank) & !is.null(tax)) # do taxonomic filtering
            {
            if(!tax_source %in% c('bins', 'bins_gtdb', 'bins_sqm')) { stop('`tax_source` has to be "bins", "bins_gtdb" or "bins_sqm", in order to subset bins based on taxonomy') }
            tax_table = get_tax_table(SQM, tax_source)
            if(!rank %in% colnames(tax_table)) { stop(sprintf('Valid taxonomic ranks are %s', paste(colnames(SQM$contigs$tax), collapse = ', '))) }
            goodBins_tax = rownames(tax_table)[tax_table[,rank] == tax]
            goodBins = intersect(goodBins, goodBins_tax)
            }
        if(!is.null(min_completeness))
            {
            goodBins_comp = rownames(SQM$bins$table)[SQM$bins$table$Completeness >= min_completeness]
            goodBins = intersect(goodBins, goodBins_comp)
            }
        if(!is.null(max_contamination))
            {
            goodBins_cont = rownames(SQM$bins$table)[SQM$bins$table$Contamination <= max_contamination]
            goodBins = intersect(goodBins, goodBins_cont)
            }
        }
    goodContigs = rownames(SQM$contigs$bins)[SQM$contigs$bins %in% goodBins]
    return ( subsetContigs(SQM, goodContigs, tax_source = tax_source,
                           trusted_functions_only = trusted_functions_only,
                           ignore_unclassified_functions=ignore_unclassified_functions,
                           rescale_tpm = rescale_tpm,
                           rescale_copy_number = rescale_copy_number,
	                   recalculate_bin_stats = FALSE,
	                   allow_empty = allow_empty)
           )
    }


#' Select contigs
#'
#' Create a SQM object containing only the requested contigs, the ORFs contained in them and the bins that contain them.
#' @param SQM SQM object to be subsetted.
#' @param contigs character. Vector of contigs to be selected.
#' @param tax_source character, source data used for the taxonomy tables present in \code{SQM$taxa}, either \code{"bins"} (GTDB bin taxonomy if available, SQM bin taxonomy otherwise), \code{"bins_gtdb"} (GTDB bin taxonomy) or \code{"bins_sqm"} (SQM bin taxonomy). Default \code{"contigs"}.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @param allow_empty (internal use only).
#' @return SQM object containing only the selected contigs.
#' @seealso \code{\link{subsetORFs}}
#' @examples
#' data(Hadza)
#' # Which contigs have a GC content below 40?
#' lowGCcontigNames = rownames(Hadza$contigs$table[Hadza$contigs$table[,"GC perc"]<40,])
#' lowGCcontigs = subsetContigs(Hadza, lowGCcontigNames)
#' hist(lowGCcontigs$contigs$table[,"GC perc"])
#' @export
subsetContigs = function(SQM, contigs, tax_source = 'contigs',
                         trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE,
                         rescale_tpm = FALSE, rescale_copy_number = FALSE,
                         recalculate_bin_stats = TRUE, allow_empty = FALSE)
    {
    return(subsetDispatch(subsetContigs_, SQM, contigs, tax_source = tax_source,
			  trusted_functions_only = trusted_functions_only,
                          ignore_unclassified_functions = ignore_unclassified_functions,
			  rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number,
                          recalculate_bin_stats = recalculate_bin_stats,
			  allow_empty = allow_empty)
          )
    }
subsetContigs_ = function(SQM, contigs, tax_source, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number, recalculate_bin_stats, allow_empty)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(SQM$misc$onlybins)
        {
        if(!tax_source %in% c('bins', 'bins_gtdb', 'bins_sqm')) { tax_source = 'bins' }
        subSQM = subsetSQM_(SQM, 'contigs', contigs, tax_source = tax_source,
                        trusted_functions_only = trusted_functions_only,
                        ignore_unclassified_functions=ignore_unclassified_functions,
                        rescale_tpm = rescale_tpm,
                        rescale_copy_number = rescale_copy_number,
                        recalculate_bin_stats = recalculate_bin_stats,
                        contigs_override = NULL,
                        allow_empty = allow_empty)
    } else
        {
        goodORFs = rownames(SQM$orfs$table)[SQM$orfs$table[,'Contig ID'] %in% contigs]
        subSQM = subsetORFs(SQM, goodORFs, tax_source = tax_source,
                            trusted_functions_only = trusted_functions_only,
                            ignore_unclassified_functions=ignore_unclassified_functions,
                            rescale_tpm = rescale_tpm,
                            rescale_copy_number = rescale_copy_number,
			    recalculate_bin_stats = recalculate_bin_stats,
	                    contigs_override = contigs,
	                    allow_empty = allow_empty)
        }
    return(subSQM)
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
#' @param tax_source character, source data used for the taxonomy tables present in \code{SQM$taxa}, either \code{"orfs"}, \code{"contigs"}, \code{"bins"} (GTDB bin taxonomy if available, SQM bin taxonomy otherwise), \code{"bins_gtdb"} (GTDB bin taxonomy) or \code{"bins_sqm"} (SQM bin taxonomy). Default \code{"orfs"}.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{FALSE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{FALSE}, which means that the returned copy numbers for each function will represent the average copy number of that function per genome in the parent object.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @param contigs_override character. Optional vector of contigs to be included in the subsetted object.
#' @param allow_empty (internal use only).
#' @return SQM object containing the requested ORFs.
#' @section A note on contig/bins subsetting:
#' While this function selects the contigs and bins that contain the desired orfs, it DOES NOT recalculate contig abundance and statistics based on the selected ORFs only. This means that the abundances presented in tables such as \code{SQM$contig$abund} will still refer to the complete contigs, regardless of whether only a fraction of their ORFs are actually present in the returned SQM object. This is also true for the statistics presented in \code{SQM$contigs$table}. Bin statistics may be recalculated if \code{rescale_copy_number} is set to \code{TRUE}, but recalculation will be based on contigs, not ORFs.
#' @examples
#' data(Hadza)
#' # Select the 100 most abundant ORFs in our dataset.
#' mostAbundantORFnames = names(sort(rowSums(Hadza$orfs$tpm), decreasing=TRUE))[1:100]
#' mostAbundantORFs = subsetORFs(Hadza, mostAbundantORFnames)
#' @importFrom stats aggregate
#' @export
subsetORFs = function(SQM, orfs, tax_source = 'orfs', trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE,
                      rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = TRUE,
                      contigs_override = NULL, allow_empty = FALSE)
    {
    return(subsetDispatch(subsetORFs_, SQM, orfs, tax_source = tax_source,
			  trusted_functions_only = trusted_functions_only,
                          ignore_unclassified_functions = ignore_unclassified_functions,
                          rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number,
                          recalculate_bin_stats = recalculate_bin_stats,
			  contigs_override = contigs_override, allow_empty = allow_empty)
          ) 
    }
subsetORFs_ = function(SQM, orfs, tax_source = 'orfs', trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, 
		      rescale_tpm = FALSE, rescale_copy_number = FALSE, recalculate_bin_stats = TRUE,
		      contigs_override = NULL, allow_empty = FALSE)
    {
    if(SQM$misc$onlybins) { stop('This function can not be run on projects generated with the `--onlybins` flag') }
    subSQM = subsetSQM_(SQM, 'orfs', orfs, tax_source = tax_source,
                        trusted_functions_only = trusted_functions_only,
                        ignore_unclassified_functions = ignore_unclassified_functions,
                        rescale_tpm = rescale_tpm, rescale_copy_number = rescale_copy_number,
                        recalculate_bin_stats = recalculate_bin_stats,
                        contigs_override = contigs_override, allow_empty = allow_empty)
    return(subSQM)
    }


subsetSQM_ = function(SQM, feature_type, features, tax_source, trusted_functions_only, ignore_unclassified_functions,
                      rescale_tpm, rescale_copy_number, recalculate_bin_stats, contigs_override, allow_empty)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(!feature_type %in% c('orfs', 'contigs')) { stop('`feature_type` has to be either "orfs" or "contigs"') }
    if(!tax_source %in% c('contigs', 'orfs', 'bins', 'bins_gtdb', 'bins_sqm')) { stop('tax_source must be "orfs", "contigs", "bins", "bins_gtdb" or "bins_sqm"') }

    ### If we got no ORFs, either fail or return an empty object
    if(length(features)==0)
        {
        if(!allow_empty) {
            if(feature_type == 'orfs') { tx = 'ORFs' } else { tx = 'contigs' }
            stop(sprintf('No %s were selected. Perhaps the subset query yielded no results?', tx))
        } else
            {
            subSQM = SQM
            subSQM$orfs = NULL
            subSQM$contigs = NULL
            subSQM$bins = NULL
            for(rank in names(subSQM$taxa))
                {
                for(count in names(subSQM$taxa[[rank]]))
                    {
                    subSQM$taxa[[rank]][[count]] = subSQM$taxa[[rank]][[count]][0,,drop=F]
                    }
                }
            for(method in names(subSQM$functions))
                {
                for(count in names(subSQM$functions[[method]]))
                    {
                    subSQM$functions[[method]][[count]] = subSQM$functions[[method]][[count]][0,,drop=F]
                    }
                }
            subSQM$misc = list(project_name = subSQM$misc$project_name,
                               samples = subSQM$misc$samples,
                               ext_annot_sources = subSQM$misc$ext_annot_sources,
                               tax_source = tax_source,
                               onlybins = subSQM$misc$onlybins,
                               single_copy_genes = subSQM$misc$single_copy_genes)
            return(subSQM)
            }
        }

    ### Now go for it
    subSQM = SQM
    # Ignore ORFs if we want to subset contigs directly
    if(feature_type == 'contigs')
        {
        contigs = features
    } else # Otherwise work with ORFs and derive the contigs from there (plus contigs_override)
        { 
        orfs    = rownames(SQM$orfs$table[features,,drop=FALSE]) # Make sure it will work if orfs is a bool vector too.
        if(is.null(contigs_override)) { contigs = unique(SQM$orfs$table[orfs,'Contig ID'])
        } else { contigs = contigs_override } # so we can include contigs without ORFs if required
    
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
        if(!is.null(subSQM$orfs$seqs))
            { 
            subSQM$orfs$seqs              = SQM$orfs$seqs   [orfs]
            }
        #    Taxonomy
        subSQM$orfs$tax                   = SQM$orfs$tax    [orfs    ,,drop=FALSE]
        subSQM$orfs$tax_abund             = aggregate_taxa(subSQM, 'orfs')
        if('tax16S' %in% names(subSQM$orfs))
            {
            subSQM$orfs$tax16S            = SQM$orfs$tax16S[intersect(orfs, names(SQM$orfs$tax16S))]
            }
        #    Markers
        if('markers' %in% names(subSQM$orfs))
            {
            subSQM$orfs$markers           = SQM$orfs$markers [orfs]
            }
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
    if(!is.null(subSQM$contigs$seqs))
        {
        subSQM$contigs$seqs           = SQM$contigs$seqs [contigs]
        }
    #    Taxonomy
    if(!is.null(subSQM$contigs$tax))
        {
        subSQM$contigs$tax            = SQM$contigs$tax  [contigs,,drop=FALSE]
        subSQM$contigs$tax_abund      = aggregate_taxa(subSQM, 'contigs')
        }
    #    Binning info
    if(!is.null(subSQM$contigs$bins))
        {
	subSQM$contigs$bins           = SQM$contigs$bins [contigs ,,drop=FALSE]
        bins                          = unique(subSQM$contigs$bins[,1])
        bins                          = bins[bins!='No bin']
	### Bins
        if(!length(bins))
            {
            msg = 'The requested subset contains no bins'
            if(subSQM$misc$onlybins)
                {
                stop(msg)
            } else
                {
                warning(msg)
                subSQM$bins = NULL
                }
        } else
            {
	    if(recalculate_bin_stats & (!('tax16S' %in% names(subSQM$orfs)) & 'markers' %in% names(subSQM$orfs)))
                {
                warning('You requested to recalculate bin stats but 16S or marker gene info are missing. Will not recalculate completeness/contamination')
                }
	    # Table and Taxonomy
            if(recalculate_bin_stats & ('tax16S' %in% names(subSQM$orfs) & 'markers' %in% names(subSQM$orfs)))
                {
                bin_stats             = get.bin.stats(subSQM)
                subSQM$bins$table     = bin_stats[['table']]
                if(!is.null(subSQM$bins$tax))
                    {
                    subSQM$bins$tax   = bin_stats[['tax']]
                    }
            }else
                {
                subSQM$bins$table     = SQM$bins$table[bins,,drop=FALSE]
                if(!is.null(subSQM$bins$tax))
                    {
                    subSQM$bins$tax   = SQM$bins$tax  [bins,,drop=FALSE]
                    }
                }
	    # Abundances
            if(recalculate_bin_stats)
                {
	        bin_abunds            = get.bin.abunds(subSQM)
	        subSQM$bins$abund     = bin_abunds[['abund']]
	        subSQM$bins$percent   = bin_abunds[['percent']]
	        subSQM$bins$bases     = bin_abunds[['bases']]
	        subSQM$bins$length    = bin_abunds[['length']]
	        subSQM$bins$cov       = bin_abunds[['cov']]
	        subSQM$bins$cpm       = bin_abunds[['cpm']]
            } else
                {
                subSQM$bins$abund     = SQM$bins$abund  [bins,,drop=FALSE]
                subSQM$bins$percent   = SQM$bins$percent[bins,,drop=FALSE]
                subSQM$bins$bases     = SQM$bins$bases  [bins,,drop=FALSE]
                subSQM$bins$length    = SQM$bins$length [bins            ]
                subSQM$bins$cov       = SQM$bins$cov    [bins,,drop=FALSE]
                subSQM$bins$cpm       = SQM$bins$cpm    [bins,,drop=FALSE]
                }
            if(!is.null(subSQM$bin$tax))
                {
                subSQM$bins$tax_abund = aggregate_taxa(subSQM, 'bins_sqm', allow_missing_annots = TRUE)
                }
            # GTDB-tax
            if(!is.null(SQM$bins$tax_gtdb))
                {
                subSQM$bins$tax_gtdb  = SQM$bins$tax_gtdb[bins,,drop=FALSE]
                subSQM$bins$tax_abund_gtdb = aggregate_taxa(subSQM, 'bins_gtdb', allow_missing_annots = TRUE)
                }
            }
        }
    ### Taxonomy
    subSQM$misc$tax_source            = tax_source
    if(!is.null(subSQM$taxa)) { subSQM$taxa = get_preferred_tax(subSQM) }
  
    ### Functions 
    if(!is.null(subSQM$functions))
        {
        if(rescale_copy_number) # get a new single_copy_gene coverage
            {
            scg = get_median_single_copy_cov(subSQM) # this will emit a warning and return NA if we can't reliably calculate single copy gene coverage
            if(!any(is.na(scg))) # we want to rescale AND can get get single copy coverage in all samples
                {
                subSQM$misc$single_copy_cov = scg
                }
            }
        for(method in names(subSQM$functions))
            {
            abunds                           = aggregate_fun(subSQM, method, trusted_functions_only, ignore_unclassified_functions)
            subSQM$functions[[method]]$abund = abunds$abund
            subSQM$functions[[method]]$bases = abunds$bases
            subSQM$functions[[method]]$cov   = abunds$cov
	    subSQM$functions[[method]]$cpm   = t(t(abunds$cov) /  (subSQM$total_reads/1000000))
            subSQM$functions[[method]]$tpm   = abunds$tpm
            if(rescale_tpm)
                {
                subSQM$functions[[method]]$tpm = abunds$tpm_rescaled
                }
            if(has_copy_numbers(subSQM))
                {
                subSQM$functions[[method]]$copy_number = t(t(subSQM$functions[[method]]$cov) / subSQM$misc$single_copy_cov)
                }
            }
        }

    return(subSQM)
    }

