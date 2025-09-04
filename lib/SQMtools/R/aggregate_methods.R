#' @importFrom stats aggregate
aggregate_bin = function(SQM)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object')}
    res = aggregate(SQM$contigs$abund, by=list(SQM$contigs$bin), FUN=sum)
    rownames(res) = res[,1]
    res = res[,-1]
    return(as.matrix(res))
    }


aggregate_taxa = function(SQM, tax_source = 'bins', allow_missing_annots = FALSE)
    {
    tk = get_tax_keys(SQM, tax_source)
    tax_source = tk[1]
    tax_key = tk[2]
    tax = SQM[[tax_source]][[tax_key]]
    abunds = SQM[[tax_source]][['abund']]
    ranks = colnames(tax)
    res = list()
    for(rank in ranks)
        {
        aggr = aggregate_taxa_rank(SQM, abunds, tax, rank, allow_missing_annots)
        percents = 100 * t(t(aggr) / SQM$total_reads)
        res[[rank]] = list(abund = aggr, percent = percents)
        }
    return(res)
    }


#' @importFrom stats aggregate
aggregate_taxa_rank = function(SQM, abunds, tax, rank, allow_missing_annots = FALSE)
    {
    missing_annots = setdiff(rownames(abunds), rownames(tax))
    missing_abunds = setdiff(rownames(tax), rownames(abunds))
    stopifnot(length(missing_abunds)==0)
    if(length(missing_annots))
        {
        # Some things like "Unclassified" or "No bin" may not be in our tax list
        # But we still want to aggregate and re-add them later
        stopifnot(allow_missing_annots)
        extra_abunds = abunds[missing_annots,,drop=F]
        abunds = abunds[!rownames(abunds) %in% missing_annots,,drop=F]
        }
    res = aggregate(data.table::data.table(abunds), by=list(tax[,rank]), FUN=sum)
    rownames(res) = res[,1] # SQM$misc$tax_names_long[[rank]][res[,1]]
    res = res[,-1,drop=FALSE]
    if(length(missing_annots)) { res = rbind(res, extra_abunds) }
    return(as.matrix(res))
    }


#' @importFrom stats aggregate
aggregate_fun = function(SQM, fun, trusted_functions_only, ignore_unclassified_functions)
    {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    stopifnot(fun %in% names(SQM$functions))
    stopifnot(is.logical(trusted_functions_only))
    stopifnot(is.logical(ignore_unclassified_functions))
    if(fun %in% c('KEGG', 'COG')) {funCol = sprintf('%s ID', fun)
    }else if(fun == 'PFAM') { funCol = 'PFAM'
    }else{funCol = fun}
    funs = SQM$orfs$table[,funCol]
    funs[funs==''] = 'Unclassified'
    if(trusted_functions_only & fun!='PFAM')
        {
        funs[!grepl('*', funs, fixed=TRUE)] = 'Unclassified'
        }
    funs = gsub('*', '', funs, fixed=TRUE)
    abund                      = aggregate(SQM$orfs$abund, by=list(funs), FUN=sum)
    rownames(abund)            = abund[,1]
    abund                      = as.matrix(abund[,-1,drop=FALSE])

    bases                      = aggregate(SQM$orfs$bases, by=list(funs), FUN=sum)
    rownames(bases)            = bases[,1]
    bases                      = as.matrix(bases[,-1,drop=FALSE])

    coverage                   = aggregate(SQM$orfs$cov  , by=list(funs), FUN=sum)
    rownames(coverage)         = coverage[,1]
    coverage                   = as.matrix(coverage[,-1,drop=FALSE])
    lengths                    = replicate(ncol(SQM$orfs$abund), SQM$orfs$table[,'Length NT'])
#    if(is.null(dim(lengths)))  { lengths = data.table::data.table(matrix(lengths, nrow=1)) } # The replicate function generates a vector if there is only one ORF. Avoid that.
    if(is.null(dim(lengths)))  { lengths = matrix(lengths, nrow=1) } # The replicate function generates a vector if there is only one ORF. Avoid that.
    dimnames(lengths)          = dimnames(SQM$orfs$abund)
    lengths[SQM$orfs$abund==0] = 0 #Don't count its length if it's missing fron that sample.
    lengths                    = aggregate(lengths       , by=list(funs), FUN=sum)
    rownames(lengths)          = lengths[,1]
    lengths                    = as.matrix(lengths[,-1,drop=FALSE])

    copies                     = aggregate(SQM$orfs$abund>0, by=list(funs), FUN=sum) # aggregate 1 if present, 0 if absent
    rownames(copies)           = copies[,1]
    copies                     = as.matrix(copies[,-1,drop=FALSE])

    tpm                        = aggregate(SQM$orfs$tpm  , by=list(funs), FUN=sum)
    rownames(tpm)              = tpm[,1]
    tpm                        = as.matrix(tpm[,-1,drop=FALSE])
    
    if(ignore_unclassified_functions)
        {
        abund    = abund[rownames(abund)      !='Unclassified',]
        bases    = bases[rownames(bases)      !='Unclassified',]
        tpm      = tpm  [rownames(tpm)        !='Unclassified',]
        coverage = coverage[rownames(coverage)!='Unclassified',]
        lengths  = lengths[rownames(lengths)  !='Unclassified',]
        copies   = copies[rownames(copies)    !='Unclassified',]
        }

    stopifnot(identical(rownames(abund), rownames(tpm)))
    stopifnot(identical(rownames(abund), rownames(bases)))
    stopifnot(identical(rownames(abund), rownames(coverage)))
    stopifnot(identical(rownames(abund), rownames(lengths)))
    stopifnot(identical(rownames(abund), rownames(copies)))

    if(fun %in% c('KEGG', 'COG', 'PFAM', 'COGonly'))
        {
	if(fun=='PFAM') { pattern = '];'
	} else { pattern = ';' }	
        multiFuns = rownames(abund)[grepl(pattern, rownames(abund), fixed=TRUE)]
        for(mf in multiFuns)
            {

            mfs = unlist(strsplit(mf, split=pattern))
	    if(fun=='PFAM') { mfs = c(paste(mfs[1:length(mfs)-1], ']', sep=''), mfs[length(mfs)]) }

            for(fun_ in mfs)
                {
                if(fun_ %in% rownames(abund))
                    {
                    abund[fun_,] = abund[fun_,] + abund[mf,] / length(mfs)
                } else
                    {
                    abund_r = abund[mf,,drop=FALSE] / length(mfs)
                    rownames(abund_r) = fun_
                    abund = rbind(abund, abund_r)
                    }

                if(fun_ %in% rownames(bases))
                    {
                    bases[fun_,] = bases[fun_,] + bases[mf,] / length(mfs)
                } else
                    {
                    bases_r = bases[mf,,drop=FALSE] / length(mfs)
                    rownames(bases_r) = fun_
                    bases = rbind(bases, bases_r)
                    }

                if(fun_ %in% rownames(coverage))
                    { coverage[fun_,] = coverage[fun_,] + coverage[mf,] / length(mfs)
                } else
                    {
                    coverage_r = coverage[mf,,drop=FALSE] / length(mfs)
                    rownames(coverage_r) = fun_
                    coverage = rbind(coverage, coverage_r)
                    }

                if(fun_ %in% rownames(lengths))
                    { lengths[fun_,] = lengths[fun_,] + lengths[mf,] / length(mfs)
                } else
                    {
                    lengths_r = lengths[mf,,drop=FALSE] / length(mfs)
                    rownames(lengths_r) = fun_
                    lengths = rbind(lengths, lengths_r)
                    }

                if(fun_ %in% rownames(copies))
                # We treat every fun in a multi-fun annotation as an individual smaller gene: less size, less reads, one copy.
                    { copies[fun_,] = copies[fun_,] + copies[mf,]
                } else
                    {
                    copies_r = copies[mf,,drop=FALSE] 
                    rownames(copies_r) = fun_
                    copies = rbind(copies, copies_r)
                    }

                if(fun_ %in% rownames(tpm))
                    { tpm[fun_,] = tpm[fun_,] + tpm[mf,] / length(mfs)
                } else
                    {
                    tpm_r = tpm[mf,,drop=FALSE] / length(mfs)
                    rownames(tpm_r) = fun_
                    tpm = rbind(tpm, tpm_r)
                    }
                }
            }

        abund        = abund   [!rownames(abund)    %in% multiFuns,,drop=FALSE]
        bases        = bases   [!rownames(bases)    %in% multiFuns,,drop=FALSE]
        coverage     = coverage[!rownames(coverage) %in% multiFuns,,drop=FALSE]
        lengths      = lengths [!rownames(lengths)  %in% multiFuns,,drop=FALSE]
        copies       = copies  [!rownames(copies)   %in% multiFuns,,drop=FALSE]
        tpm          = tpm     [!rownames(tpm)      %in% multiFuns,,drop=FALSE]
        }

    #avgLengths   = lengths / copies
    #rpk          = 1000 * abund/avgLengths
    #rpk[is.na(rpk)] = 0
    #tpm_rescaled = 1000000 * t(t(rpk)/colSums(rpk))

    tpm_rescaled = 1000000 * t(t(tpm)/colSums(tpm))

    abund        = abund   [sort(rownames(abund))           ,,drop=FALSE]
    bases        = bases   [sort(rownames(bases))           ,,drop=FALSE]
    coverage     = coverage[sort(rownames(coverage))        ,,drop=FALSE]
    tpm          = tpm     [sort(rownames(tpm))             ,,drop=FALSE]
    tpm_rescaled = tpm_rescaled[sort(rownames(tpm_rescaled)),,drop=FALSE]

    tpm          = t(t(tpm) * SQM$misc$coding_fraction[[fun]])

    return(list(abund=round(abund), bases=round(bases), cov=coverage, tpm=tpm, tpm_rescaled=tpm_rescaled))

    }
