#' Combine several SQM objects
#'
#' Combine an arbitrary number of SQM objects into a single SQM object. The input objects must be subsets of the same original SQM object (i.e. from the same SqueezeMeta run). For combining results from different runs please check \code{\link[combineSQMlite]{combineSQMlite}}.
#' @param ... an arbitrary number of SQM objects. Alternatively, a single list containing an arbitrary number of SQM objects.
#' @param tax_source character. Features used for calculating aggregated abundances at the different taxonomic ranks. Either \code{"orfs"} or \code{"contigs"} (default \code{"orfs"}). If the objects being combined contain a subset of taxa or bins, this parameter can be set to \code{TRUE}.
#' @param trusted_functions_only logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).
#' @param ignore_unclassified_functions logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).
#' @param rescale_tpm logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{TRUE}).
#' @param rescale_copy_number logical. If \code{TRUE}, copy numbers with be recalculated using the RecA/RadA coverages in the subset. Otherwise, RecA/RadA coverages will be taken from the parent object with the highest RecA/RadA coverages. By default it is set to \code{TRUE}, which means that the returned copy numbers will represent the average copy number per function \emph{in the genomes of the selected bins or contigs}. If any SQM objects that are being combined contain a functional subset rather than a contig/bins subset, this parameter should be set to \code{FALSE}.
#' @return A SQM object
#' @seealso \code{\link[subsetFun]{subsetFun}}, \code{\link[subsetTax]{subsetTax}}, \code{\link[combineSQMlite]{combineSQMlite}}
#' @examples
#' data(Hadza)
#' # Select Carbohydrate metabolism ORFs in Bacteroidetes, and Amino acid metabolism ORFs in Proteobacteria
#' bact = subsetTax(Hadza, "phylum", "Bacteroidetes")
#' bact.carb = subsetFun(bact, "Carbohydrate metabolism")
#' proteo = subsetTax(Hadza, "phylum", "Proteobacteria")
#' proteo.amins = subsetFun(proteo, "Amino acid metabolism")
#' bact.carb_proteo.amins = combineSQM(bact.carb, proteo.amins, rescale_copy_number=F)
#' @export
combineSQM = function(..., tax_source = 'orfs', trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = T, rescale_copy_number = T)
    {
    inSQM = list(...)
    # if there is only one argument and this argument is a list, treat it as a list containing SQM objects
    if(length(inSQM) == 1 & class(inSQM[[1]]) == 'list') { inSQM = list(...)[[1]] }
    # intermediate function so that we can pass extra args to combineSQM
    myFun = function(SQM1, SQM2) combineSQM_(SQM1, SQM2, tax_source, trusted_functions_only, ignore_unclassified_functions, rescale_tpm, rescale_copy_number)
    return(Reduce(myFun, inSQM))
    }


combineSQM_ = function(SQM1, SQM2, tax_source = 'orfs', trusted_functions_only = F, ignore_unclassified_functions = F, rescale_tpm = T, rescale_copy_number = T)
    {

    if(class(SQM1) != 'SQM' | class(SQM2) != 'SQM') { stop('This function only accepts SQM objects') }

    if (!identical(colnames(SQM1$orfs$table), colnames(SQM2$orfs$table)))
        {
        stop('The input objects do not seem to come from the same SQM project')
        }

    combSQM = SQM1
    ### ORFs
    extraORFs                          = setdiff(rownames(SQM2$orfs$table), rownames(SQM1$orfs$table))
    #    Table
    combSQM$orfs$table                 = rbind(combSQM$orfs$table, SQM2$orfs$table[extraORFs,,drop=F])
    combSQM$orfs$table                 = combSQM$orfs$table[sort(rownames(combSQM$orfs$table)),,drop=F]
    #    Abundances
    combSQM$orfs$abund                 = rbind(combSQM$orfs$abund, SQM2$orfs$abund[extraORFs,,drop=F])
    combSQM$orfs$abund                 = combSQM$orfs$abund[rownames(combSQM$orfs$table),,drop=F]
    combSQM$orfs$bases                 = rbind(combSQM$orfs$bases, SQM2$orfs$bases[extraORFs,,drop=F])
    combSQM$orfs$bases                 = combSQM$orfs$bases[rownames(combSQM$orfs$table),,drop=F]
    combSQM$orfs$cov                   = rbind(combSQM$orfs$cov, SQM2$orfs$cov[extraORFs,,drop=F])
    combSQM$orfs$cov                   = combSQM$orfs$cov[rownames(combSQM$orfs$table),,drop=F]
    combSQM$orfs$cpm                   = t(t(combSQM$orfs$cov) / (combSQM$total_reads / 1000000))
    combSQM$orfs$tpm                   = rbind(combSQM$orfs$tpm, SQM2$orfs$tpm[extraORFs,,drop=F])
    combSQM$orfs$tpm                   = combSQM$orfs$tpm[rownames(combSQM$orfs$table),,drop=F]
    #    Sequences
    combSQM$orfs$seqs                  = c(combSQM$orfs$seqs, SQM2$orfs$seqs[extraORFs])
    combSQM$orfs$seqs                  = combSQM$orfs$seqs[rownames(combSQM$orfs$table)]
    #    Taxonomy
    combSQM$orfs$tax                   = rbind(combSQM$orfs$tax, SQM2$orfs$tax[extraORFs,,drop=F])
    combSQM$orfs$tax                   = combSQM$orfs$tax[rownames(combSQM$orfs$table),,drop=F]
    
    ### Contigs
    extraContigs                       = setdiff(rownames(SQM2$contigs$table), rownames(SQM1$contigs$table))
    #    Table
    combSQM$contigs$table              = rbind(combSQM$contigs$table, SQM2$contigs$table[extraContigs,,drop=F])
    combSQM$contigs$table              = combSQM$contigs$table[sort(rownames(combSQM$contigs$table)),]
    #    Abundances
    combSQM$contigs$abund              = rbind(combSQM$contigs$abund, SQM2$contigs$abund[extraContigs,,drop=F])
    combSQM$contigs$abund              = combSQM$contigs$abund[rownames(combSQM$contigs$table),,drop=F]
    combSQM$contigs$bases              = rbind(combSQM$contigs$bases, SQM2$contigs$bases[extraContigs,,drop=F])
    combSQM$contigs$bases              = combSQM$contigs$bases[rownames(combSQM$contigs$table),,drop=F]
    combSQM$contigs$cov                = rbind(combSQM$contigs$cov, SQM2$contigs$cov[extraContigs,,drop=F])
    combSQM$contigs$cov                = combSQM$contigs$cov[rownames(combSQM$contigs$table),,drop=F]
    combSQM$contigs$cpm                = t(t(combSQM$contigs$cov) / (combSQM$total_reads / 1000000))
    combSQM$contigs$tpm                = rbind(combSQM$contigs$tpm, SQM2$contigs$tpm[extraContigs,,drop=F])
    combSQM$contigs$tpm                = combSQM$contigs$tpm[rownames(combSQM$contigs$table),,drop=F]
    #    Sequences
    combSQM$contigs$seqs               = c(combSQM$contigs$seqs, SQM2$contigs$seqs[extraContigs])
    combSQM$contigs$seqs               = combSQM$contigs$seqs[rownames(combSQM$contigs$table)]
    #    Taxonomy
    combSQM$contigs$tax                = rbind(combSQM$contigs$tax, SQM2$contigs$tax[extraContigs,,drop=F])
    combSQM$contigs$tax                = combSQM$contigs$tax[rownames(combSQM$contigs$table),,drop=F]
    #    Binning info
    if('bins' %in% names(combSQM))
        {
        combSQM$contigs$bins           = rbind(combSQM$contigs$bins, SQM2$contigs$bins[extraContigs,,drop=F])
        combSQM$contigs$bins           = combSQM$contigs$bins[rownames(combSQM$contigs$table),,drop=F]

    ### Bins

        extraBins                      = setdiff(rownames(SQM2$bins$table), rownames(SQM1$bins$table))
        #    Table
        combSQM$bins$table             = rbind(combSQM$bins$table, SQM2$bins$table[extraBins,,drop=F])
        combSQM$bins$table             = combSQM$bins$table[sort(rownames(combSQM$bins$table)),,drop=F]
        #    Abundances
	x = aggregate(combSQM$contigs$abund, by=list(combSQM$contigs$bins[,1]), FUN=sum)
        rownames(x)                    = x[,1]
        x = x[rownames(combSQM$bin$table),-1]
        nobin                          = colSums(combSQM$contigs$abund) - colSums(x)
        if(sum(nobin)>0)               { x['No_bin',] = nobin }

        combSQM$bins$abund             = as.matrix(x)
        combSQM$bins$percent           = 100 * t(t(combSQM$bins$abund) / combSQM$total_reads)

        x = aggregate(combSQM$contigs$bases, by=list(combSQM$contigs$bins[,1]), FUN=sum)
        rownames(x)                    = x[,1]
        x = x[rownames(combSQM$bin$table),-1]
        combSQM$bins$bases             = as.matrix(x)

        l = aggregate(combSQM$contigs$table$Length, by=list(combSQM$contigs$bins[,1]), FUN=sum)
        n = l[,1]; l = l[,-1]; names(l) = n
        l = l[rownames(combSQM$bin$table)]
        combSQM$bins$length            = l
        combSQM$bins$cov               = combSQM$bins$bases / combSQM$bins$length
	combSQM$bins$cpm               = t(t(combSQM$bins$cov) / (combSQM$total_reads / 1000000))

        #    Taxonomy
        combSQM$bins$tax               = rbind(combSQM$bins$tax, SQM2$bins$tax[extraBins,,drop=F])
        combSQM$bins$tax               = combSQM$bins$tax[rownames(combSQM$bins$table),,drop=F]
        }

    ### Taxonomy   
    combSQM$taxa$superkingdom$abund    = aggregate.taxa(combSQM, 'superkingdom', tax_source)
    combSQM$taxa$phylum$abund          = aggregate.taxa(combSQM, 'phylum'      , tax_source)
    combSQM$taxa$class$abund           = aggregate.taxa(combSQM, 'class'       , tax_source)
    combSQM$taxa$order$abund           = aggregate.taxa(combSQM, 'order'       , tax_source)
    combSQM$taxa$family$abund          = aggregate.taxa(combSQM, 'family'      , tax_source)
    combSQM$taxa$genus$abund           = aggregate.taxa(combSQM, 'genus'       , tax_source)
    combSQM$taxa$species$abund         = aggregate.taxa(combSQM, 'species'     , tax_source)

    combSQM$taxa$superkingdom$percent  = 100 * t(t(combSQM$taxa$superkingdom$abund) / combSQM$total_reads) #colSums(combSQM$taxa$superkingdom$abund))
    combSQM$taxa$phylum$percent        = 100 * t(t(combSQM$taxa$phylum$abund)       / combSQM$total_reads) #colSums(combSQM$taxa$phylum$abund))
    combSQM$taxa$class$percent         = 100 * t(t(combSQM$taxa$class$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$class$abund))
    combSQM$taxa$order$percent         = 100 * t(t(combSQM$taxa$order$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$order$abund))
    combSQM$taxa$family$percent        = 100 * t(t(combSQM$taxa$family$abund)       / combSQM$total_reads) #colSums(combSQM$taxa$family$abund))
    combSQM$taxa$genus$percent         = 100 * t(t(combSQM$taxa$genus$abund)        / combSQM$total_reads) #colSums(combSQM$taxa$genus$abund))
    combSQM$taxa$species$percent       = 100 * t(t(combSQM$taxa$species$abund)      / combSQM$total_reads) #colSums(combSQM$taxa$species$abund))

    ### Functions
    if('KEGG' %in% names(combSQM$functions))
        {
        KEGG                           = aggregate.fun(combSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$KEGG$abund   = KEGG$abund
        combSQM$functions$KEGG$bases   = KEGG$bases
        combSQM$functions$KEGG$cov     = KEGG$cov
	combSQM$functions$KEGG$cpm     = t(t(combSQM$functions$KEGG$cov) / (combSQM$total_reads / 1000000))
        }

    if('COG' %in% names(combSQM$functions))
        {
        COG                            = aggregate.fun(combSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$COG$abund    = COG$abund
        combSQM$functions$COG$bases    = COG$bases
        combSQM$functions$COG$cov      = COG$cov
	combSQM$functions$COG$cpm      = t(t(combSQM$functions$COG$cov) / (combSQM$total_reads / 1000000))
        }


    if('PFAM' %in% names(combSQM$functions))
        {
        PFAM                           = aggregate.fun(combSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
        combSQM$functions$PFAM$abund   = PFAM$abund
        combSQM$functions$PFAM$bases   = PFAM$bases
        combSQM$functions$PFAM$cov     = PFAM$cov
	combSQM$functions$PFAM$cpm     = t(t(combSQM$functions$PFAM$cov) / (combSQM$total_reads / 1000000))
        }

    ext_annots = list()
    for(method in combSQM$misc$ext_annot_sources)
        {
        ext_annots[[method]]              = aggregate.fun(combSQM, method, trusted_functions_only, ignore_unclassified_functions)
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


    if(!is.null(combSQM$misc$RecA_cov))
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
		    RecA = pmax(SQM1$misc$RecA_cov, SQM2$misc$RecA_cov) # use the largest and hope for the best.
                    }
            }else
                {
                warning('RecA is not present in this subset. Will not rescale copy numbers.')
                RecA = pmax(SQM1$misc$RecA_cov, SQM2$misc$RecA_cov) # use the largest and hope for the best.
                }
        }else
            {
            RecA = pmax(SQM1$misc$RecA_cov, SQM2$misc$RecA_cov) # use the largest and hope for the best.
            }
        if('KEGG' %in% names(combSQM$functions)) { combSQM$functions$KEGG$copy_number = t(t(KEGG$cov) / RecA) }
        if('COG' %in% names(combSQM$functions))  { combSQM$functions$COG$copy_number  = t(t(COG$cov ) / RecA) }
        if('PFAM' %in% names(combSQM$functions)) { combSQM$functions$PFAM$copy_number = t(t(PFAM$cov) / RecA) }
        for(method in combSQM$misc$ext_annot_sources)
            {
            combSQM$functions[[method]]$copy_number = t(t(ext_annots[[method]]$cov) / RecA)
            }

        combSQM$misc$RecA_cov              = RecA
        }

    ### Total reads
    #combSQM$total_reads               = colSums(combSQM$contigs$abund)

    return(combSQM)
    }


