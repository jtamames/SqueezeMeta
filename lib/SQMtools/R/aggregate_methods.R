aggregate.bin = function(SQM)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object')}
    res = aggregate(SQM$contigs$abund, by=list(SQM$contigs$bin), FUN=sum)
    rownames(res) = res[,1]
    res = res[,-1]
    return(as.matrix(res))
    }



aggregate.taxa = function(SQM, rank, tax_source)
    {
    stopifnot(tax_source %in% c('contigs', 'orfs'))
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object')}
    res = aggregate(SQM[[tax_source]][['abund']], by=list(SQM[[tax_source]][['tax']][,rank]), FUN=sum)
    rownames(res) = res[,1] # SQM$misc$tax_names_long[[rank]][res[,1]]
    res = res[,-1,drop=F]
    return(as.matrix(res))
    }



aggregate.fun = function(SQM, fun, trusted_functions_only, ignore_unclassified_functions)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    stopifnot(fun %in% c('KEGG', 'COG', 'PFAM'))
    stopifnot(is.logical(trusted_functions_only))
    stopifnot(is.logical(ignore_unclassified_functions))
    if(fun %in% c('KEGG', 'COG')) {funCol = sprintf('%s ID', fun)
    }else{funCol = 'PFAM'}
    
    funs = SQM$orfs$table[,funCol]
    funs[funs==''] = 'Unclassified'
    if(trusted_functions_only & fun!='PFAM')
        {
        funs[!grepl('*', funs, fixed=T)] = 'Unclassified'
        }
    funs = gsub('*', '', funs, fixed=T)

    abund                      = aggregate(SQM$orfs$abund, by=list(funs), FUN=sum)
    rownames(abund)            = abund[,1]
    abund                      = as.matrix(abund[,-1,drop=F])

    coverage                   = aggregate(SQM$orfs$cov  , by=list(funs), FUN=sum)
    rownames(coverage)         = coverage[,1]
    coverage                   = as.matrix(coverage[,-1,drop=F])

    lengths                    = replicate(ncol(SQM$orfs$abund), SQM$orfs$table[,'Length NT'])
    dimnames(lengths)          = dimnames(SQM$orfs$abund)
    lengths[SQM$orfs$abund==0] = 0 #Don't count its length if it's missing fron that sample.
    lengths                    = aggregate(lengths       , by=list(funs), FUN=sum)
    rownames(lengths)          = lengths[,1]
    lengths                    = as.matrix(lengths[,-1,drop=F])

    copies                     = aggregate(SQM$orfs$abund>0, by=list(funs), FUN=sum) # aggregate 1 if present, 0 if absent
    rownames(copies)           = copies[,1]
    copies                     = as.matrix(copies[,-1,drop=F])

    tpm                        = aggregate(SQM$orfs$tpm  , by=list(funs), FUN=sum)
    rownames(tpm)              = tpm[,1]
    tpm                        = as.matrix(tpm[,-1,drop=F])

    if(ignore_unclassified_functions)
        {
        abund    = abund[rownames(abund)      !='Unclassified',]
        tpm      = tpm  [rownames(tpm)        !='Unclassified',]
        coverage = coverage[rownames(coverage)!='Unclassified',]
        lengths  = lengths[rownames(lengths)  !='Unclassified',]
        copies   = copies[rownames(copies)    !='Unclassified',]
        }

    stopifnot(identical(rownames(abund), rownames(tpm)))
    stopifnot(identical(rownames(abund), rownames(coverage)))
    stopifnot(identical(rownames(abund), rownames(lengths)))
    stopifnot(identical(rownames(abund), rownames(copies)))

    multiFuns = rownames(abund)[grepl(';', rownames(abund), fixed=T)]
    for(mf in multiFuns)
        {
        mfs = unlist(strsplit(mf, split=';'))

        for(fun_ in mfs)
            {
            if(fun_ %in% rownames(abund))
                {
                abund[fun_,] = abund[fun_,] + abund[mf,] / length(mfs)
            } else
                {
                abund_r = abund[mf,,drop=F] / length(mfs)
                rownames(abund_r) = fun_
                abund = rbind(abund, abund_r)
                }

            if(fun_ %in% rownames(coverage))
                { coverage[fun_,] = coverage[fun_,] + coverage[mf,] / length(mfs)
            } else
                {
                coverage_r = coverage[mf,,drop=F] / length(mfs)
                rownames(coverage_r) = fun_
                coverage = rbind(coverage, coverage_r)
                }

            if(fun_ %in% rownames(lengths))
                { lengths[fun_,] = lengths[fun_,] + lengths[mf,] / length(mfs)
            } else
                {
                lengths_r = lengths[mf,,drop=F] / length(mfs)
                rownames(lengths_r) = fun_
                lengths = rbind(lengths, lengths_r)
                }

            if(fun_ %in% rownames(copies))
            # We treat every fun in a multi-fun annotation as an individual smaller gene: less size, less reads, one copy.
                { copies[fun_,] = copies[fun_,] + copies[mf,]
            } else
                {
                copies_r = copies[mf,,drop=F] 
                rownames(copies_r) = fun_
                copies = rbind(copies, copies_r)
                }

            if(fun_ %in% rownames(tpm))
                { tpm[fun_,] = tpm[fun_,] + tpm[mf,] / length(mfs)
            } else
                {
                tpm_r = tpm[mf,,drop=F] / length(mfs)
                rownames(tpm_r) = fun_
                tpm = rbind(tpm, tpm_r)
                }
            }
        }

    abund        = abund   [!rownames(abund)    %in% multiFuns,,drop=F]
    coverage     = coverage[!rownames(coverage) %in% multiFuns,,drop=F]
    lengths      = lengths [!rownames(lengths)  %in% multiFuns,,drop=F]
    copies       = copies  [!rownames(copies)   %in% multiFuns,,drop=F]
    tpm          = tpm     [!rownames(tpm)      %in% multiFuns,,drop=F]

    avgLengths   = lengths / copies
    rpk          = 1000 * abund/avgLengths
    rpk[is.na(rpk)] = 0
    tpm_rescaled = 1000000 * t(t(rpk)/colSums(rpk))

    abund        = abund   [sort(rownames(abund))           ,,drop=F]
    coverage     = coverage[sort(rownames(coverage))        ,,drop=F]
    tpm          = tpm     [sort(rownames(tpm))             ,,drop=F]
    tpm_rescaled = tpm_rescaled[sort(rownames(tpm_rescaled)),,drop=F]

    return(list(abund=round(abund), cov=coverage, tpm=tpm, tpm_rescaled=tpm_rescaled))

    }
