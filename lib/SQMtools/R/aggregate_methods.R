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

    abund              = aggregate(SQM$orfs$abund, by=list(funs), FUN=sum)
    rownames(abund)    = abund[,1]
    abund              = as.matrix(abund[,-1,drop=F])
    coverage           = aggregate(SQM$orfs$cov  , by=list(funs), FUN=sum)
    rownames(coverage) = coverage[,1]
    coverage           = as.matrix(coverage[,-1,drop=F])
    tpm                = aggregate(SQM$orfs$tpm  , by=list(funs), FUN=sum)
    rownames(tpm)      = tpm[,1]
    tpm                = as.matrix(tpm[,-1,drop=F])

    if(ignore_unclassified_functions)
        {
        abund    = abund[rownames(abund)      !='Unclassified',]
        tpm      = tpm  [rownames(tpm)        !='Unclassified',]
        coverage = coverage[rownames(coverage)!='Unclassified']
        }

    stopifnot(identical(rownames(abund), rownames(tpm)))
    stopifnot(identical(rownames(abund), rownames(coverage)))
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

    abund    = abund   [!rownames(abund)    %in% multiFuns,,drop=F]
    coverage = coverage[!rownames(coverage) %in% multiFuns,,drop=F]
    tpm      = tpm     [!rownames(tpm)      %in% multiFuns,,drop=F]
    abund    = abund   [sort(rownames(abund))   ,,drop=F]
    coverage = coverage[sort(rownames(coverage)),,drop=F]
    tpm      = tpm     [sort(rownames(tpm))     ,,drop=F]

    return(list(abund=round(abund), cov=coverage, tpm=tpm))

    }
