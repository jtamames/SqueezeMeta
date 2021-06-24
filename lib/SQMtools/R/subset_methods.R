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
#' @seealso \code{\link[subsetTax]{subsetTax}}, \code{\link[subsetORFs]{subsetORFs}}, \code{\link[combineSQM]{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be eselected with \code{\link[mostAbundant]{mostAbundant}}.
#' @return SQM object containing only the requested function.
#' @examples
#' data(Hadza)
#' Hadza.iron = subsetFun(Hadza, "iron")
#' Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
#' # Search for multiple patterns using regular expressions
#' Hadza.twoKOs = subsetFun(Hadza, "K00812|K00813", fixed=F)
#' @export
subsetFun = function(SQM, fun, columns = NULL, ignore_case=T, fixed=F, trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = F, rescale_copy_number = F)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
	
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
                        rescale_copy_number = rescale_copy_number)
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
#' @return SQM object containing only the requested taxon.
#' @seealso \code{\link[subsetFun]{subsetFun}}, \code{\link[subsetContigs]{subsetContigs}}, \code{\link[combineSQM]{combineSQM}}. The most abundant items of a particular table contained in a SQM object can be eselected with \code{\link[mostAbundant]{mostAbundant}}.
#' @examples
#' data(Hadza)
#' Hadza.Escherichia = subsetTax(Hadza, "genus", "Escherichia")
#' Hadza.Bacteroidetes = subsetTax(Hadza, "phylum", "Bacteroidetes")
#' @export
subsetTax = function(SQM, rank, tax, trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm =T, rescale_copy_number = T)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    goodContigs = rownames(SQM$contigs$tax)[SQM$contigs$tax[,rank] == tax]
    return ( subsetContigs(SQM, goodContigs,
                           trusted_functions_only = trusted_functions_only,
                           ignore_unclassified_functions=ignore_unclassified_functions,
                           rescale_tpm = rescale_tpm,
                           rescale_copy_number = rescale_copy_number)
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
#' @seealso \code{\link[subsetContigs]{subsetContigs}}, \code{\link[subsetORFs]{subsetORFs}}
#' @examples 
#' data(Hadza)
#' # Which are the two most complete bins?
#' topBinNames = rownames(Hadza$bins$table)[order(Hadza$bins$table[,"Completeness"], decreasing=T)][1:2]
#' topBins = subsetBins(Hadza, topBinNames)
#' @export
subsetBins = function(SQM, bins, trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = T, rescale_copy_number = T)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    goodContigs = rownames(SQM$contigs$bins)[SQM$contigs$bins %in% bins]
    return ( subsetContigs(SQM, goodContigs,  
                           trusted_functions_only = trusted_functions_only,
                           ignore_unclassified_functions=ignore_unclassified_functions,
                           rescale_tpm = rescale_tpm,
                           rescale_copy_number = rescale_copy_number)
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
#' @return SQM object containing only the selected contigs.
#' @seealso \code{\link[subsetORFs]{subsetORFs}}
#' @examples
#' data(Hadza)
#' # Which contigs have a GC content below 40?
#' lowGCcontigNames = rownames(Hadza$contigs$table[Hadza$contigs$table[,"GC perc"]<40,])
#' lowGCcontigs = subsetContigs(Hadza, lowGCcontigNames)
#' hist(lowGCcontigs$contigs$table[,"GC perc"])
#' @export
subsetContigs = function(SQM, contigs, trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = F, rescale_copy_number = F)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    goodORFs = rownames(SQM$orfs$table)[SQM$orfs$table[,'Contig ID'] %in% contigs]
    return ( subsetORFs(SQM, goodORFs, tax_source = 'contigs',
                        trusted_functions_only = trusted_functions_only,
                        ignore_unclassified_functions=ignore_unclassified_functions,
                        rescale_tpm = rescale_tpm,
                        rescale_copy_number = rescale_copy_number,
	                contigs_override = contigs)
           )
    }


#' Select random ORFs
#'
#' Create a random subset of a SQM object.
#' @param SQM SQM object to be subsetted.
#' @param N numeric. number of random ORFs to select.
#' @return SQM object containing a random subset of ORFs.
#' @seealso \code{\link[subsetORFs]{subsetORFs}}
#' @export
subsetRand = function(SQM, N)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
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
#' @return SQM object containing the requested ORFs.
#' @section A note on contig/bins subsetting:
#' While this function selects the contigs and bins that contain the desired orfs, it DOES NOT recalculate contig/bin abundance and statistics based on the selected ORFs only. This means that the abundances presented in tables such as \code{SQM$contig$abund} or \code{SQM$bins$tpm} will still refer to the complete contigs and bins, regardless of whether only a fraction of their ORFs are actually present in the returned SQM object. This is also true for the statistics presented in \code{SQM$contigs$table} and \code{SQM$bins$table}.
#' @examples
#' data(Hadza)
#' # Select the 100 most abundant ORFs in our dataset.
#' mostAbundantORFnames = names(sort(rowSums(Hadza$orfs$tpm), decreasing=T))[1:100]
#' mostAbundantORFs = subsetORFs(Hadza, mostAbundantORFnames)
#' @export
subsetORFs = function(SQM, orfs, tax_source = 'orfs', trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = F, rescale_copy_number = F, contigs_override = NULL)
    {

    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    if(length(orfs)==0) { stop('No ORFs were selected. Perhaps the subset query yielded no results?') }
    if(!tax_source %in% c('contigs', 'orfs')) { stop('tax_source must be "orfs" or "contigs"') }
   
    orfs    = rownames(SQM$orfs$table[orfs,,drop=F]) # Make sure it will work if orfs is a bool vector too.
    if(is.null(contigs_override)) { contigs = unique(SQM$orfs$table[orfs,'Contig ID'])
    } else { contigs = contigs_override } # so we can include contigs without ORFs if required
    bins    = unique( unlist(SQM$contigs$bins[contigs,]) )
    bins    = bins[bins!='No_bin']
    
    subSQM = SQM

    subSQM$orfs$table                 = SQM$orfs$table[orfs      ,,drop=F]
    subSQM$orfs$abund                 = SQM$orfs$abund[orfs      ,,drop=F]
    subSQM$orfs$bases                 = SQM$orfs$bases[orfs      ,,drop=F]
    subSQM$orfs$cov                   = SQM$orfs$cov[orfs        ,,drop=F]
    subSQM$orfs$tpm                   = SQM$orfs$tpm[orfs        ,,drop=F]
    subSQM$orfs$seqs                  = SQM$orfs$seqs[orfs]
    subSQM$orfs$tax                   = SQM$orfs$tax[orfs        ,,drop=F]

    subSQM$contigs$table              = SQM$contigs$table[contigs,,drop=F]
    subSQM$contigs$abund              = SQM$contigs$abund[contigs,,drop=F]
    subSQM$contigs$cov                = SQM$contigs$cov[contigs  ,,drop=F]
    subSQM$contigs$tpm                = SQM$contigs$tpm[contigs  ,,drop=F]
    subSQM$contigs$seqs               = SQM$contigs$seqs[contigs]
    subSQM$contigs$tax                = SQM$contigs$tax[contigs  ,,drop=F]
    if('bins' %in% names(subSQM))
        {
        subSQM$contigs$bins           = SQM$contigs$bins[contigs ,,drop=F]
        subSQM$bins$table             = subSQM$bins$table[bins   ,,drop=F]
        subSQM$bins$tpm               = subSQM$bins$tpm[bins     ,,drop=F]
        subSQM$bins$tax               = subSQM$bins$tax[bins     ,,drop=F]
        }

    subSQM$taxa$superkingdom$abund    = aggregate.taxa(subSQM, 'superkingdom', tax_source)
    subSQM$taxa$phylum$abund          = aggregate.taxa(subSQM, 'phylum'      , tax_source)
    subSQM$taxa$class$abund           = aggregate.taxa(subSQM, 'class'       , tax_source)
    subSQM$taxa$order$abund           = aggregate.taxa(subSQM, 'order'       , tax_source)
    subSQM$taxa$family$abund          = aggregate.taxa(subSQM, 'family'      , tax_source)
    subSQM$taxa$genus$abund           = aggregate.taxa(subSQM, 'genus'       , tax_source)
    subSQM$taxa$species$abund         = aggregate.taxa(subSQM, 'species'     , tax_source)

    subSQM$taxa$superkingdom$percent  = 100 * t(t(subSQM$taxa$superkingdom$abund) / subSQM$total_reads) #colSums(subSQM$taxa$superkingdom$abund))
    subSQM$taxa$phylum$percent        = 100 * t(t(subSQM$taxa$phylum$abund)       / subSQM$total_reads) #colSums(subSQM$taxa$phylum$abund))
    subSQM$taxa$class$percent         = 100 * t(t(subSQM$taxa$class$abund)        / subSQM$total_reads) #colSums(subSQM$taxa$class$abund))
    subSQM$taxa$order$percent         = 100 * t(t(subSQM$taxa$order$abund)        / subSQM$total_reads) #colSums(subSQM$taxa$order$abund))
    subSQM$taxa$family$percent        = 100 * t(t(subSQM$taxa$family$abund)       / subSQM$total_reads) #colSums(subSQM$taxa$family$abund))
    subSQM$taxa$genus$percent         = 100 * t(t(subSQM$taxa$genus$abund)        / subSQM$total_reads) #colSums(subSQM$taxa$genus$abund))
    subSQM$taxa$species$percent       = 100 * t(t(subSQM$taxa$species$abund)      / subSQM$total_reads) #colSums(subSQM$taxa$species$abund))

    if('KEGG' %in% names(subSQM$functions))
        {
        KEGG                          = aggregate.fun(subSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$KEGG$abund   = KEGG$abund
        subSQM$functions$KEGG$bases   = KEGG$bases
        subSQM$functions$KEGG$cov     = KEGG$cov
        }
    
    if('COG' %in% names(subSQM$functions))
        {
        COG                           = aggregate.fun(subSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$COG$abund    = COG$abund
        subSQM$functions$COG$bases    = COG$bases
        subSQM$functions$COG$cov      = COG$cov
        }


    if('PFAM' %in% names(subSQM$functions))
        {
        PFAM                          = aggregate.fun(subSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions$PFAM$abund   = PFAM$abund
        subSQM$functions$PFAM$bases   = PFAM$bases
        subSQM$functions$PFAM$cov     = PFAM$cov
        }

    ext_annots = list()
    for(method in subSQM$misc$ext_annot_sources)
        {
        ext_annots[[method]]          = aggregate.fun(subSQM, method, trusted_functions_only, ignore_unclassified_functions)
        subSQM$functions[[method]]$abund = ext_annots[[method]]$abund
        subSQM$functions[[method]]$bases = ext_annots[[method]]$bases
        subSQM$functions[[method]]$cov   = ext_annots[[method]]$cov
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
        #if('bins' in names(subSQM)
        #    {
        #    subSQM$bins$tpm           = 1000000 * t(t(subSQM$bins$tpm)   /colSums(subSQM$bins$tpm)   )
        #    }
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

    return(subSQM)
    }

