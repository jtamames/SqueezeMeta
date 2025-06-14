#' Combine several SQM objects
#'
#' Combine an arbitrary number of SQM objects into a single SQM object (if the input objects contain the same samples, i.e. they come from the same SqueezeMeta run) or a single SQMbunch object. For combining results from sqm_reads.pl or sqm_longreads.pl please check \code{\link{combineSQMlite}}. The parameters below (other than ...) will take only effect if the input objects contain the same samples. Otherwise the input objects will be taken as they are, with no recalculation of taxonomy, function or rescaling,
#' @param ... an arbitrary number of SQM objects. Alternatively, a single list containing an arbitrary number of SQM objects.
#' @param tax_source character. Features used for calculating aggregated abundances at the different taxonomic ranks. Either \code{"orfs"} or \code{"contigs"} (default \code{"orfs"}). If the objects being combined contain a subset of taxa or bins, this parameter can be set to \code{TRUE}.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{TRUE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers will represent the average copy number per function \emph{in the genomes of the selected bins or contigs}. If any SQM objects that are being combined contain a functional subset rather than a contig/bins subset, this parameter should be set to \code{FALSE}.
#' @param recalculate_bin_stats logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).
#' @return A SQM or SQMbunch object
#' @seealso \code{\link{subsetFun}}, \code{\link{subsetTax}}, \code{\link{combineSQMlite}}
#' @examples
#' data(Hadza)
#' # Select Carbohydrate metabolism ORFs in Bacteroidota,
#' #  and Amino acid metabolism ORFs in Proteobacteria
#' bact = subsetTax(Hadza, "phylum", "Bacteroidota")
#' bact.carb = subsetFun(bact, "Carbohydrate metabolism")
#' baci = subsetTax(Hadza, "phylum", "Bacillota")
#' baci.amins = subsetFun(baci, "Amino acid metabolism")
#' bact.carb_proteo.amins = combineSQM(bact.carb, baci.amins, rescale_copy_number=FALSE)
#' @export
combineSQM = function(..., tax_source = 'orfs', trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE, recalculate_bin_stats = TRUE)
    {
    inSQM = list(...)
    ### if there is only one argument and this argument is a list, treat it as a list containing SQM objects
    if(length(inSQM) == 1 & inherits(inSQM[[1]], 'list')) { inSQM = list(...)[[1]] }
    ### check whether the input objects come from the same project or not
    projNames = sapply(inSQM, FUN=function(proj) proj$misc$project_name)
    if(length(unique(projNames))==1)
        {
        # intermediate function so that we can pass extra args to combineSQM
        myFun = function(SQM1, SQM2) combineSQM_(SQM1, SQM2, tax_source, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number, recalculate_bin_stats)
        combSQM = Reduce(myFun, inSQM)
    } else
        {
	combSQM = combineSQMlite(inSQM)
	combSQM$projects = inSQM
	names(combSQM$projects) = projNames
        class(combSQM) = 'SQMbunch'
        }
    return(combSQM)
    }


#' @importFrom stats aggregate
combineSQM_ = function(SQM1, SQM2, tax_source = 'orfs', trusted_functions_only = FALSE, ignore_unclassified_functions = FALSE, rescale_tpm = TRUE, rescale_copy_number = TRUE, recalculate_bin_stats = TRUE)
    {

    if(!inherits(SQM1, 'SQM') | !inherits(SQM2, 'SQM')) { stop('This function only accepts SQM objects') }

    if (!identical(colnames(SQM1$orfs$table), colnames(SQM2$orfs$table)))
        {
        stop('The input objects do not seem to come from the same SQM project')
        }
    
    combSQM = SQM1
    ### ORFs
    extraORFs                          = setdiff(rownames(SQM2$orfs$table), rownames(SQM1$orfs$table))
    #    Table
    combSQM$orfs$table                 = rbind(combSQM$orfs$table, SQM2$orfs$table[extraORFs,,drop=FALSE])
    combSQM$orfs$table                 = combSQM$orfs$table[sort(rownames(combSQM$orfs$table)),,drop=FALSE]
    #    Abundances
    combSQM$orfs$abund                 = rbind(combSQM$orfs$abund, SQM2$orfs$abund[extraORFs,,drop=FALSE])
    combSQM$orfs$abund                 = combSQM$orfs$abund[rownames(combSQM$orfs$table),,drop=FALSE]
    combSQM$orfs$bases                 = rbind(combSQM$orfs$bases, SQM2$orfs$bases[extraORFs,,drop=FALSE])
    combSQM$orfs$bases                 = combSQM$orfs$bases[rownames(combSQM$orfs$table),,drop=FALSE]
    combSQM$orfs$cov                   = rbind(combSQM$orfs$cov, SQM2$orfs$cov[extraORFs,,drop=FALSE])
    combSQM$orfs$cov                   = combSQM$orfs$cov[rownames(combSQM$orfs$table),,drop=FALSE]
    combSQM$orfs$cpm                   = t(t(combSQM$orfs$cov) / (combSQM$total_reads / 1000000))
    combSQM$orfs$tpm                   = rbind(combSQM$orfs$tpm, SQM2$orfs$tpm[extraORFs,,drop=FALSE])
    combSQM$orfs$tpm                   = combSQM$orfs$tpm[rownames(combSQM$orfs$table),,drop=FALSE]
    #    Sequences
    if(!is.null(combSQM$orfs$seqs))
        {
        combSQM$orfs$seqs              = c(combSQM$orfs$seqs, SQM2$orfs$seqs[extraORFs])
        combSQM$orfs$seqs              = combSQM$orfs$seqs[rownames(combSQM$orfs$table)]
        }
    #    Taxonomy
    combSQM$orfs$tax                   = rbind(combSQM$orfs$tax, SQM2$orfs$tax[extraORFs,,drop=FALSE])
    combSQM$orfs$tax                   = combSQM$orfs$tax[rownames(combSQM$orfs$table),,drop=FALSE]
    if('tax16S' %in% names(combSQM$orfs))
        {
        combSQM$orfs$tax16S            = c(combSQM$orfs$tax16S, SQM2$orfs$tax16S[extraORFs])
	combSQM$orfs$tax16S            = combSQM$orfs$tax16S[rownames(combSQM$orfs$table)]
        }
    #    Markers
    if('markers' %in% names(combSQM$orfs))
        {
        combSQM$orfs$markers               = c(combSQM$orfs$markers, SQM2$orfs$markers[extraORFs])
        combSQM$orfs$markers               = combSQM$orfs$markers[rownames(combSQM$orfs$table)]
        }
    ### Contigs
    extraContigs                       = setdiff(rownames(SQM2$contigs$table), rownames(SQM1$contigs$table))
    #    Table
    combSQM$contigs$table              = rbind(combSQM$contigs$table, SQM2$contigs$table[extraContigs,,drop=FALSE])
    combSQM$contigs$table              = combSQM$contigs$table[sort(rownames(combSQM$contigs$table)),]
    #    Abundances
    combSQM$contigs$abund              = rbind(combSQM$contigs$abund, SQM2$contigs$abund[extraContigs,,drop=FALSE])
    combSQM$contigs$abund              = combSQM$contigs$abund[rownames(combSQM$contigs$table),,drop=FALSE]
    combSQM$contigs$bases              = rbind(combSQM$contigs$bases, SQM2$contigs$bases[extraContigs,,drop=FALSE])
    combSQM$contigs$bases              = combSQM$contigs$bases[rownames(combSQM$contigs$table),,drop=FALSE]
    combSQM$contigs$cov                = rbind(combSQM$contigs$cov, SQM2$contigs$cov[extraContigs,,drop=FALSE])
    combSQM$contigs$cov                = combSQM$contigs$cov[rownames(combSQM$contigs$table),,drop=FALSE]
    combSQM$contigs$cpm                = t(t(combSQM$contigs$cov) / (combSQM$total_reads / 1000000))
    combSQM$contigs$tpm                = rbind(combSQM$contigs$tpm, SQM2$contigs$tpm[extraContigs,,drop=FALSE])
    combSQM$contigs$tpm                = combSQM$contigs$tpm[rownames(combSQM$contigs$table),,drop=FALSE]
    #    Sequences
    if(!is.null(combSQM$contigs$seqs))
        {
        combSQM$contigs$seqs           = c(combSQM$contigs$seqs, SQM2$contigs$seqs[extraContigs])
        combSQM$contigs$seqs           = combSQM$contigs$seqs[rownames(combSQM$contigs$table)]
        }
    #    Taxonomy
    combSQM$contigs$tax                = rbind(combSQM$contigs$tax, SQM2$contigs$tax[extraContigs,,drop=FALSE])
    combSQM$contigs$tax                = combSQM$contigs$tax[rownames(combSQM$contigs$table),,drop=FALSE]
    #    Binning info
    if('bins' %in% names(combSQM))
        {
        combSQM$contigs$bins           = rbind(combSQM$contigs$bins, SQM2$contigs$bins[extraContigs,,drop=FALSE])
        combSQM$contigs$bins           = combSQM$contigs$bins[rownames(combSQM$contigs$table),,drop=FALSE]
	### Bins
        if(recalculate_bin_stats & (!('tax16S' %in% names(combSQM$orfs)) & 'markers' %in% names(combSQM$orfs)))
            {
            warning('You requested to recalculate bin stats but 16S or marker gene info are missing. Will not recalculate completeness/contamination')
            }
	# Table and Taxonomy
        extraBins                     = setdiff(rownames(SQM2$bins$table), rownames(SQM1$bins$table))
        if(recalculate_bin_stats & ('tax16S' %in% names(combSQM$orfs) & 'markers' %in% names(combSQM$orfs)))
            {
            bin_stats                 = get.bin.stats(combSQM)
            combSQM$bins$table        = bin_stats[['table']]
            combSQM$bins$tax          = bin_stats[['tax']] 
        }else
            {
            combSQM$bins$table        = rbind(combSQM$bins$table, SQM2$bins$table[extraBins,,drop=FALSE])
            combSQM$bins$table        = combSQM$bins$table[sort(rownames(combSQM$bins$table)),,drop=FALSE]
            combSQM$bins$tax          = rbind(combSQM$bins$tax, SQM2$bins$tax[extraBins,,drop=FALSE])
            combSQM$bins$tax          = combSQM$bins$tax[rownames(combSQM$bins$table),,drop=FALSE]
            }
        # Abundances
        if(recalculate_bin_stats)
            {
            bin_abunds                = get.bin.abunds(combSQM)
            combSQM$bins$abund        = bin_abunds[['abund']]
            combSQM$bins$percent      = bin_abunds[['percent']]
            combSQM$bins$bases        = bin_abunds[['bases']]
            combSQM$bins$length       = bin_abunds[['length']]
            combSQM$bins$cov          = bin_abunds[['cov']]
            combSQM$bins$cpm          = bin_abunds[['cpm']]
        } else
            {
            combSQM$bins$abund        = rbind(combSQM$bins$abund,   SQM2$bins$abund  [extraBins,,drop=FALSE])
            combSQM$bins$percent      = rbind(combSQM$bins$percent, SQM2$bins$percent[extraBins,,drop=FALSE])
            combSQM$bins$bases        = rbind(combSQM$bins$bases,   SQM2$bins$bases  [extraBins,,drop=FALSE])
            combSQM$bins$length       = c    (combSQM$bins$length,  SQM2$bins$length [extraBins            ])
            combSQM$bins$cov          = rbind(combSQM$bins$cov,     SQM2$bins$cov    [extraBins,,drop=FALSE])
            combSQM$bins$cpm          = rbind(combSQM$bins$cpm,     SQM2$bins$cpm    [extraBins,,drop=FALSE])
            }
        # GTDB-tax
        if(!is.null(SQM1$bins$tax_gtdb) | !is.null(SQM2$bins$tax_gtdb))
            {
            combSQM$bins$tax_gtdb     = rbind(combSQM$bins$tax_gtdb, SQM2$bins$tax_gtdb[extraBins,,drop=FALSE])
            combSQM$bins$tax_gtdb     = combSQM$bins$tax_gtdb[rownames(combSQM$bins$table),,drop=FALSE]
            }
        }

    ### Taxonomy   
    combSQM$taxa$superkingdom$abund   = aggregate_taxa(combSQM, 'superkingdom', tax_source)
    combSQM$taxa$phylum$abund         = aggregate_taxa(combSQM, 'phylum'      , tax_source)
    combSQM$taxa$class$abund          = aggregate_taxa(combSQM, 'class'       , tax_source)
    combSQM$taxa$order$abund          = aggregate_taxa(combSQM, 'order'       , tax_source)
    combSQM$taxa$family$abund         = aggregate_taxa(combSQM, 'family'      , tax_source)
    combSQM$taxa$genus$abund          = aggregate_taxa(combSQM, 'genus'       , tax_source)
    combSQM$taxa$species$abund        = aggregate_taxa(combSQM, 'species'     , tax_source)

    combSQM$taxa$superkingdom$percent = 100 * t(t(combSQM$taxa$superkingdom$abund) / combSQM$total_reads) #colSums(combSQM$taxa$superkingdom$abund))
    combSQM$taxa$phylum$percent       = 100 * t(t(combSQM$taxa$phylum$abund)       / combSQM$total_reads) #colSums(combSQM$taxa$phylum$abund))
    combSQM$taxa$class$percent        = 100 * t(t(combSQM$taxa$class$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$class$abund))
    combSQM$taxa$order$percent        = 100 * t(t(combSQM$taxa$order$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$order$abund))
    combSQM$taxa$family$percent       = 100 * t(t(combSQM$taxa$family$abund)       / combSQM$total_reads) #colSums(combSQM$taxa$family$abund))
    combSQM$taxa$genus$percent        = 100 * t(t(combSQM$taxa$genus$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$genus$abund))
    combSQM$taxa$species$percent      = 100 * t(t(combSQM$taxa$species$abund)      / combSQM$total_reads) #colSums(combSQM$taxa$species$abund))

    ### Functions
    if('KEGG' %in% names(combSQM$functions))
        {
        KEGG                          = aggregate_fun(combSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$KEGG$abund  = KEGG$abund
        combSQM$functions$KEGG$bases  = KEGG$bases
        combSQM$functions$KEGG$cov    = KEGG$cov
	combSQM$functions$KEGG$cpm    = t(t(combSQM$functions$KEGG$cov) / (combSQM$total_reads / 1000000))
        }

    if('COG' %in% names(combSQM$functions))
        {
        COG                           = aggregate_fun(combSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$COG$abund   = COG$abund
        combSQM$functions$COG$bases   = COG$bases
        combSQM$functions$COG$cov     = COG$cov
	combSQM$functions$COG$cpm     = t(t(combSQM$functions$COG$cov) / (combSQM$total_reads / 1000000))
        }


    if('PFAM' %in% names(combSQM$functions))
        {
        PFAM                          = aggregate_fun(combSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$PFAM$abund  = PFAM$abund
        combSQM$functions$PFAM$bases  = PFAM$bases
        combSQM$functions$PFAM$cov    = PFAM$cov
	combSQM$functions$PFAM$cpm    = t(t(combSQM$functions$PFAM$cov) / (combSQM$total_reads / 1000000))
        }

    ext_annots = list()
    for(method in combSQM$misc$ext_annot_sources)
        {
        ext_annots[[method]]              = aggregate_fun(combSQM, method, trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions[[method]]$abund = ext_annots[[method]]$abund
        combSQM$functions[[method]]$bases = ext_annots[[method]]$bases
        combSQM$functions[[method]]$cov   = ext_annots[[method]]$cov
	combSQM$functions[[method]]$cpm   = t(t(combSQM$functions[[method]]$cov) / (combSQM$total_reads / 1000000))
        }

    if(rescale_tpm)
        {
        if('KEGG' %in% names(combSQM$functions)) { combSQM$functions$KEGG$tpm = KEGG$tpm_rescaled }
        if('COG' %in% names(combSQM$functions))  { combSQM$functions$COG$tpm  = COG$tpm_rescaled  }
        if('PFAM' %in% names(combSQM$functions)) { combSQM$functions$PFAM$tpm = PFAM$tpm_rescaled }
         for(method in combSQM$misc$ext_annot_sources)
            {
            combSQM$functions[[method]]$tpm = ext_annots[[method]]$tpm_rescaled
            }

    }else
        {
        if('KEGG' %in% names(combSQM$functions)) { combSQM$functions$KEGG$tpm = KEGG$tpm }
        if('COG' %in% names(combSQM$functions))  { combSQM$functions$COG$tpm  = COG$tpm  }
        if('PFAM' %in% names(combSQM$functions)) { combSQM$functions$PFAM$tpm = PFAM$tpm }
        for(method in combSQM$misc$ext_annot_sources)
            {
            combSQM$functions[[method]]$tpm = ext_annots[[method]]$tpm
            }
        }
    
    ### COPY NUMBERS
    if(has_copy_numbers(SQM1) | has_copy_numbers(SQM2))
        {
        will_rescale_cg = FALSE
        if(rescale_copy_number)
            {
            scg = get_median_single_copy_cov(combSQM) # this will emit a warning and return NA if we can't reliably calculate single copy gene coverage
            if(!any(is.na(scg))) # we want to rescale AND can get get single copy coverage in all samples
                {
                will_rescale_cg = TRUE
                combSQM$misc$single_copy_cov = scg
                }
            }
        if(!will_rescale_cg)
            {
            warning('    Single copy gene coverage will not be recalculated. Instead, we will use the highest values present in the input objects')
            combSQM$misc$single_copy_cov = pmax(SQM1$misc$single_copy_cov, SQM2$misc$single_copy_cov)
            }

        for(method in names(combSQM$functions))
            {
            combSQM$functions[[method]]$copy_number = t(t(combSQM$functions[[method]]$cov) / combSQM$misc$single_copy_cov)
            }
        }
    return(combSQM)
    }

