require(reshape2)

expand_proteo = function(phylum_abund, class_abund)
    {
    proteoall = phylum_abund['Proteobacteria',]
    proteoclass = colSums(class_abund[grep('proteobacteria',rownames(class_abund)),])
    proteodiff = proteoall-proteoclass
    rownames(proteodiff) = c('Unclassified Proteobacteria')

    proteo = class_abund[grep('proteobacteria',rownames(class_abund)),]
    proteo = rbind(proteo, proteodiff)

    res = phylum_abund[rownames(phylum_abund)!='Proteobacteria',]
    res = rbind(res, proteo)
    stopifnot(all(colSums(res)==colSums(class_abund)))
    return(res)
    }



read.namedvector = function(file)
    {
    ta = read.table(file, header=T, row.names=1, as.is=T)
    res = ta[,1]
    names(res) = rownames(ta)
    return(res)
    }



aggregate.taxa = function(SQM, rank)
    {
    res = aggregate(SQM$contigs$abund, by=list(SQM$contigs$tax[,rank]), FUN=sum)
    rownames(res) = SQM$misc$long_names[[rank]][res[,1]]
    res = res[,-1]
    return(as.matrix(res))
    }



aggregate.fun = function(SQM, fun, trusted_functions_only, ignore_unclassified_functions)
    {
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

    abund           = aggregate(SQM$orfs$abund, by=list(funs), FUN=sum)
    rownames(abund) = abund[,1]
    abund           = as.matrix(abund[,-1])
    tpm             = aggregate(SQM$orfs$tpm  , by=list(funs), FUN=sum)
    rownames(tpm)   = tpm[,1]
    tpm             = as.matrix(tpm[,-1])

    if(ignore_unclassified_functions)
        {
        abund = abund[rownames(abund)!='Unclassified',]
        tpm   = tpm  [rownames(tpm)  !='Unclassified',]
        }

    stopifnot(identical(rownames(abund),rownames(tpm)))
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

    abund = abund[!rownames(abund) %in% multiFuns,,drop=F]
    tpm   = tpm  [!rownames(tpm)   %in% multiFuns,,drop=F]

    abund = abund[sort(rownames(abund)),,drop=F]
    tpm   = tpm  [sort(rownames(tpm))  ,,drop=F]

    return(list(abund=round(abund), tpm=tpm))

    }



subsetFun = function(SQM, fun, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    goodORFs = rownames(SQM$orfs$table)[
                                        grepl(fun, SQM$orfs$table[,'KEGG ID'],  fixed=T) |
                                        grepl(fun, SQM$orfs$table[,'KEGGFUN'],  fixed=T) |
                                        grepl(fun, SQM$orfs$table[,'KEGGPATH'], fixed=T) |
                                        grepl(fun, SQM$orfs$table[,'COG ID'],   fixed=T) |
                                        grepl(fun, SQM$orfs$table[,'COGFUN'],   fixed=T) |
                                        grepl(fun, SQM$orfs$table[,'COGPATH'],  fixed=T) | 
                                        grepl(fun, SQM$orfs$table[,'PFAM'],     fixed=T)
                                       ]
    return ( subsetORFs(SQM, goodORFs) )
    }




subsetTax = function(SQM, rank, tax, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    goodContigs = rownames(SQM$contigs$tax)[SQM$contigs$tax[,rank] == tax]   
    return ( subsetContigs(SQM, goodContigs) )
    }




subsetBins = function(SQM, bins, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    goodContigs = names(SQM$contigs$bins)[SQM$contigs$bins %in% bins]
    return ( subsetContigs(SQM, goodContigs) )
    }




subsetContigs = function(SQM, contigs, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    goodORFs = rownames(SQM$orfs$table)[SQM$orfs$table[,'Contig ID'] %in% contigs]
    return ( subsetORFs(SQM, goodORFs) )
    }




subsetORFs = function(SQM, orfs, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    
    orfs    = rownames(SQM$orfs$table[orfs,]) # Make sure it will work if orfs is a bool vector too.
    contigs = unique(SQM$orfs$table[orfs,'Contig ID'])
    bins    = unique( unlist(SQM$contigs$bins[contigs,]) )
    bins    = bins[bins!='No_bin']
    
    subSQM = SQM

    subSQM$orfs$table                = SQM$orfs$table[orfs,]
    subSQM$orfs$abund                = SQM$orfs$abund[orfs,]
    subSQM$orfs$tpm                  = SQM$orfs$tpm[orfs,]
    subSQM$orfs$seqs                 = SQM$orfs$seqs[orfs]
    subSQM$orfs$tax                  = SQM$orfs$tax[orfs,]

    subSQM$contigs$table             = SQM$contigs$table[contigs,]
    subSQM$contigs$abund             = SQM$contigs$abund[contigs,]
    subSQM$contigs$tpm               = SQM$contigs$tpm[contigs,]
    subSQM$contigs$seqs              = SQM$contigs$seqs[contigs]
    subSQM$contigs$tax               = SQM$contigs$tax[contigs,]
    subSQM$contigs$bins              = SQM$contigs$bins[contigs,,drop=F]

    subSQM$bins$table                = subSQM$bins$table[bins,]
    subSQM$bins$tpm                  = subSQM$bins$tpm[bins,]
    subSQM$bins$tax                  = subSQM$bins$tax[bins,]

    subSQM$taxa$superkingdom$abund   = aggregate.taxa(subSQM, 'superkingdom')
    subSQM$taxa$phylum$abund         = aggregate.taxa(subSQM, 'phylum')
    subSQM$taxa$class$abund          = aggregate.taxa(subSQM, 'class')
    subSQM$taxa$order$abund          = aggregate.taxa(subSQM, 'order')
    subSQM$taxa$family$abund         = aggregate.taxa(subSQM, 'family')
    subSQM$taxa$genus$abund          = aggregate.taxa(subSQM, 'genus')
    subSQM$taxa$species$abund        = aggregate.taxa(subSQM, 'species')

    subSQM$taxa$superkingdom$percent = 100 * t(t(subSQM$taxa$superkingdom$abund) / colSums(subSQM$taxa$superkingdom$abund))
    subSQM$taxa$phylum$percent       = 100 * t(t(subSQM$taxa$phylum$abund)       / colSums(subSQM$taxa$phylum$abund))
    subSQM$taxa$class$percent        = 100 * t(t(subSQM$taxa$class$abund)        / colSums(subSQM$taxa$class$abund))
    subSQM$taxa$order$percent        = 100 * t(t(subSQM$taxa$order$abund)        / colSums(subSQM$taxa$order$abund))
    subSQM$taxa$family$percent       = 100 * t(t(subSQM$taxa$family$abund)       / colSums(subSQM$taxa$family$abund))
    subSQM$taxa$genus$percent        = 100 * t(t(subSQM$taxa$genus$abund)        / colSums(subSQM$taxa$genus$abund))
    subSQM$taxa$species$percent      = 100 * t(t(subSQM$taxa$species$abund)      / colSums(subSQM$taxa$species$abund))

    KEGG                             = aggregate.fun(subSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
    COG                              = aggregate.fun(subSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
    PFAM                             = aggregate.fun(subSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
    subSQM$functions$KEGG$abund      = KEGG$abund
    subSQM$functions$KEGG$tpm        = KEGG$tpm
    subSQM$functions$COG$abund       = COG$abund
    subSQM$functions$COG$tpm         = COG$tpm
    subSQM$functions$PFAM$abund      = PFAM$abund
    subSQM$functions$PFAM$tpm        = PFAM$tpm

    subSQM$total_reads		     = colSums(subSQM$contigs$abund)

    return(subSQM)

    }




combineSQM = function(..., trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    # intermediate function so that we can pass extra args to combineSQM
    myFun = function(SQM1, SQM2) combineSQM_(SQM1, SQM2, trusted_functions_only, ignore_unclassified_functions)
    return(Reduce(myFun, list(...)))
    }




combineSQM_ = function(SQM1, SQM2, trusted_functions_only = F, ignore_unclassified_functions = F)
    {
    stopifnot(identical(colnames(SQM1$orfs$table), colnames(SQM2$orfs$table)))
    combSQM = SQM1

    ### ORFs
    extraORFs                         = setdiff(rownames(SQM2$orfs$table), rownames(SQM1$orfs$table))
    #    Table
    combSQM$orfs$table                = rbind(combSQM$orfs$table, SQM2$orfs$table[extraORFs,])
    combSQM$orfs$table                = combSQM$orfs$table[sort(rownames(combSQM$orfs$table)),]
    #    Abundances
    combSQM$orfs$abund                = as.matrix(combSQM$orfs$table[,grepl('Raw read count', colnames(combSQM$orfs$table))])
    colnames(combSQM$orfs$abund)      = gsub('Raw read count ', '', colnames(combSQM$orfs$abund), fixed=T)
    combSQM$orfs$tpm                  = as.matrix(combSQM$orfs$table[,grepl('TPM', colnames(combSQM$orfs$table))])
    colnames(combSQM$orfs$tpm)        = gsub('TPM ', '', colnames(combSQM$orfs$tpm), fixed=T)
    #    Sequences
    combSQM$orfs$sequences            = c(combSQM$orfs$sequences, SQM2$orfs$sequences[extraORFs])
    combSQM$orfs$sequences            = combSQM$orfs$sequences[rownames(combSQM$orfs$table)]
    #    Taxonomy
    combSQM$orfs$tax                  = rbind(combSQM$orfs$tax, SQM2$orfs$tax[extraORFs,])
    combSQM$orfs$tax                  = combSQM$orfs$tax[rownames(combSQM$orfs$table),]
    
    ### Contigs
    extraContigs                      = setdiff(rownames(SQM2$contigs$table), rownames(SQM1$contigs$table))
    #    Table
    combSQM$contigs$table             = rbind(combSQM$contigs$table, SQM2$contigs$table[extraContigs,])
    combSQM$contigs$table             = combSQM$contigs$table[sort(rownames(combSQM$contigs$table)),]
    #    Abundances
    combSQM$contigs$abund             = as.matrix(combSQM$contigs$table[,grepl('Raw read count', colnames(combSQM$contigs$table))])
    colnames(combSQM$contigs$abund)   = gsub('Raw read count ', '', colnames(combSQM$contigs$abund), fixed=T)
    combSQM$contigs$tpm               = as.matrix(combSQM$contigs$table[,grepl('TPM', colnames(combSQM$contigs$table))])
    colnames(combSQM$contigs$tpm)     = gsub('TPM ', '', colnames(combSQM$contigs$tpm), fixed=T)
    #    Sequences
    combSQM$orfs$sequences            = c(combSQM$contigs$sequences, SQM2$contigs$sequences[extraContigs])
    combSQM$orfs$sequences            = combSQM$contigs$sequences[rownames(combSQM$contigs$table)]
    #    Taxonomy
    combSQM$contigs$tax               = rbind(combSQM$contigs$tax, SQM2$contigs$tax[extraContigs,])
    combSQM$contigs$tax               = combSQM$contigs$tax[rownames(combSQM$contigs$table),]
    #    Binning info
    combSQM$contigs$bins              = rbind(combSQM$contigs$bins, SQM2$contigs$bins[extraContigs,,drop=F])
    combSQM$contigs$bins              = combSQM$contigs$bins[rownames(combSQM$contigs$table),,drop=F]

    ### Bins
    extraBins                         = setdiff(rownames(SQM2$bins$table), rownames(SQM1$bins$table))
    #    Table
    combSQM$bins$table                = rbind(combSQM$bins$table, SQM2$bins$table[extraBins,])
    combSQM$bins$table                = combSQM$bins$table[sort(rownames(combSQM$bins$table)),]
    #    Abundances
    combSQM$bins$tpm                  = as.matrix(combSQM$bins$table[,grepl('TPM', colnames(combSQM$bins$table))])
    colnames(combSQM$bins$tpm)        = gsub('TPM ', '', colnames(combSQM$bins$tpm), fixed=T)
    #    Taxonomy
    combSQM$bins$tax                  = rbind(combSQM$bins$tax, SQM2$bins$tax[extraBins,])
    combSQM$bins$tax                  = combSQM$bins$tax[rownames(combSQM$bins$table),]

    ### Taxonomy   
    combSQM$taxa$superkingdom$abund   = aggregate.taxa(combSQM, 'superkingdom')
    combSQM$taxa$phylum$abund         = aggregate.taxa(combSQM, 'phylum')
    combSQM$taxa$class$abund          = aggregate.taxa(combSQM, 'class')
    combSQM$taxa$order$abund          = aggregate.taxa(combSQM, 'order')
    combSQM$taxa$family$abund         = aggregate.taxa(combSQM, 'family')
    combSQM$taxa$genus$abund          = aggregate.taxa(combSQM, 'genus')
    combSQM$taxa$species$abund        = aggregate.taxa(combSQM, 'species')

    combSQM$taxa$superkingdom$percent = 100 * t(t(combSQM$taxa$superkingdom$abund) / colSums(combSQM$taxa$superkingdom$abund))
    combSQM$taxa$phylum$percent       = 100 * t(t(combSQM$taxa$phylum$abund)       / colSums(combSQM$taxa$phylum$abund))
    combSQM$taxa$class$percent        = 100 * t(t(combSQM$taxa$class$abund)        / colSums(combSQM$taxa$class$abund))
    combSQM$taxa$order$percent        = 100 * t(t(combSQM$taxa$order$abund)        / colSums(combSQM$taxa$order$abund))
    combSQM$taxa$family$percent       = 100 * t(t(combSQM$taxa$family$abund)       / colSums(combSQM$taxa$family$abund))
    combSQM$taxa$genus$percent        = 100 * t(t(combSQM$taxa$genus$abund)        / colSums(combSQM$taxa$genus$abund))
    combSQM$taxa$species$percent      = 100 * t(t(combSQM$taxa$species$abund)      / colSums(combSQM$taxa$species$abund))

    ### Functions
    KEGG                              = aggregate.fun(combSQM, 'KEGG', trusted_functions_only, ignore_unclassified_functions)
    COG                               = aggregate.fun(combSQM, 'COG' , trusted_functions_only, ignore_unclassified_functions)
    PFAM                              = aggregate.fun(combSQM, 'PFAM', trusted_functions_only, ignore_unclassified_functions)
    combSQM$functions$KEGG$abund      = KEGG$abund
    combSQM$functions$KEGG$tpm        = KEGG$tpm
    combSQM$functions$COG$abund       = COG$abund
    combSQM$functions$COG$tpm         = COG$tpm
    combSQM$functions$PFAM$abund      = PFAM$abund
    combSQM$functions$PFAM$tpm        = PFAM$tpm

    ### Total reads
    combSQM$total_reads               = colSums(combSQM$contigs$abund)

    return(combSQM)
    }


exportTable = function(table, outName)
    {
    write.table(table, outName, col.names=NA, sep='\t', quote=F)
    }


loadSQM = function(project_path, tax_mode = 'allfilter')
    {
    
    if(!tax_mode %in% c('allfilter', 'prokfilter'))
        {
        stop('tax_mode must be either \'allfilter\' (apply minimum identity threshold for all taxa) or \'prokfilter\' (don\'t appy thresholds to Eukaryotes)')
        }

    # include my hierarchy somehow?? -> sqm2tables, aggregate KEGG hierarchy levels??
    
    SQM                           = list()
    
    project_name                  = tail(unlist(strsplit(project_path, split='/')), 1)
    SQM$misc                      = list()
    SQM$misc$projectName          = project_name

    cat('Loading orfs\n')
    SQM$orfs                      = list()

    cat('    table...\n')         # option to remove table from memory after getting abund & TPM?
    SQM$orfs$table                = read.table(sprintf('%s/results/13.%s.orftable', project_path, project_name),
                                               header=T, sep='\t', row.names=1, quote='', comment.char='', skip=1, as.is=TRUE, check.names=F)
    cat('    abundances...\n')
    SQM$orfs$abund                = as.matrix(SQM$orfs$table[,grepl('Raw read count', colnames(SQM$orfs$table))])
    colnames(SQM$orfs$abund)      = gsub('Raw.read.count ', '', colnames(SQM$orfs$abund), fixed=T)
    SQM$orfs$tpm                  = as.matrix(SQM$orfs$table[,grepl('TPM', colnames(SQM$orfs$table))])
    colnames(SQM$orfs$tpm)        = gsub('TPM ', '', colnames(SQM$orfs$tpm), fixed=T)
    SQM$misc$samples              = colnames(SQM$orfs$abund)
    
    cat('    sequences\n')    
    SQM$orfs$seqs                 = read.namedvector(sprintf('%s/results/tables/%s.orf.sequences.tsv', project_path, project_name))
    
    cat('    taxonomy...\n')
    SQM$orfs$tax                  = as.matrix(read.table(sprintf('%s/results/tables/%s.orf.tax.%s.tsv', project_path, project_name, tax_mode),
                                                         header=T, row.names=1, sep='\t'))
    # Remove orfs with no nt length (which should be fixed at some point). The tax table contains the correct number of orfs.
    SQM$orfs$table                = SQM$orfs$table[rownames(SQM$orfs$table) %in% rownames(SQM$orfs$tax),]
    SQM$orfs$abund                = SQM$orfs$abund[rownames(SQM$orfs$table),]
    SQM$orfs$tpm                  = SQM$orfs$tpm[rownames(SQM$orfs$table),]
    SQM$orfs$seqs                 = SQM$orfs$seqs[rownames(SQM$orfs$table)]


    cat('Loading contigs\n')
    SQM$contigs                   = list()

    cat('    table...\n')                    # option to remove table from memory after getting abund & TPM?
    SQM$contigs$table             = read.table(sprintf('%s/results/20.%s.contigtable', project_path, project_name),
                                                       header=T, sep='\t', row.names=1, quote='', comment.char='', skip=1, as.is=T, check.names=F)
    cat('    abundances...\n')
    SQM$contigs$abund             = as.matrix(SQM$contigs$table[,grepl('Raw read count', colnames(SQM$contigs$table))])
    colnames(SQM$contigs$abund)   = gsub('Raw read count ', '', colnames(SQM$contigs$abund), fixed=T)
    SQM$contigs$tpm               = as.matrix(SQM$contigs$table[,grepl('TPM', colnames(SQM$contigs$table))])
    colnames(SQM$contigs$tpm)     = gsub('TPM ', '', colnames(SQM$contigs$tpm), fixed=T)

    cat('    sequences...\n')                                                 
    SQM$contigs$seqs              = read.namedvector(sprintf('%s/results/tables/%s.contig.sequences.tsv', project_path, project_name))
    SQM$contigs$seqs              = SQM$contigs$seqs[rownames(SQM$contigs$table)]

    cat('    taxonomy...\n')
    SQM$contigs$tax               = as.matrix(read.table(sprintf('%s/results/tables/%s.contig.tax.tsv', project_path, project_name),
                                                         header=T, row.names=1, sep='\t'))
    SQM$contigs$tax               = SQM$contigs$tax[rownames(SQM$contigs$table),]

    cat('    binning info...\n')
    inBins                        = read.table(sprintf('%s/intermediate/19.%s.contigsinbins', project_path, project_name),
                                               header=T, sep='\t', quote='', comment.char='', skip=1, as.is=T)
    inBins                        = dcast(inBins, X..Contig~Method, value.var="Bin.ID")
    rownames(inBins)              = inBins[,1]
    SQM$contigs$bins              = as.matrix(inBins[,-1,drop=F])
    notInBins                     = setdiff(rownames(SQM$contigs$table), SQM$contigs$bins)
    notInBins                     = matrix(NA, nrow=length(notInBins), ncol=ncol(SQM$contigs$bins), dimnames=list(notInBins, colnames(SQM$contigs$bins)))
    SQM$contigs$bins              = rbind(SQM$contigs$bins, notInBins)
    SQM$contigs$bins              = SQM$contigs$bins[rownames(SQM$contigs$table),,drop=F]
    SQM$contigs$bins[is.na(SQM$contigs$bins)] = 'No_bin'


    cat('Loading bins\n')
    cat('    table...\n')
    SQM$bins                      = list()
    SQM$bins$table                = read.table(sprintf('%s/results/19.%s.bintable', project_path, project_name),
                                               header=T, sep='\t', row.names=1, quote='', comment.char='', skip=1, as.is=T, check.names=F)
    cat('    abundances...\n')
    SQM$bins$tpm                  = as.matrix(SQM$bins$table[,grepl('TPM', colnames(SQM$bins$table))])
    colnames(SQM$bins$tpm)        = gsub('TPM ', '', colnames(SQM$bins$tpm), fixed=T)
    cat('    taxonomy...\n')
    SQM$bins$taxonomy             = as.matrix(read.table(sprintf('%s/results/tables/%s.bin.tax.tsv', project_path, project_name),
                                                         header=T, row.names=1, sep='\t'))
    cat('Loading taxonomies\n')                                   
    SQM$taxa                      = list()
    SQM$taxa$superkingdom         = list()
    SQM$taxa$phylum               = list()
    SQM$taxa$class                = list()
    SQM$taxa$order                = list()
    SQM$taxa$family               = list()
    SQM$taxa$genus                = list()
    SQM$taxa$species              = list()
    

    SQM$taxa$superkingdom$abund   = as.matrix(read.table(sprintf('%s/results/tables/%s.superkingdom.%s.abund.tsv', project_path, project_name, tax_mode, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$phylum$abund         = as.matrix(read.table(sprintf('%s/results/tables/%s.phylum.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$class$abund          = as.matrix(read.table(sprintf('%s/results/tables/%s.class.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$order$abund          = as.matrix(read.table(sprintf('%s/results/tables/%s.order.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$family$abund         = as.matrix(read.table(sprintf('%s/results/tables/%s.family.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$genus$abund          = as.matrix(read.table(sprintf('%s/results/tables/%s.genus.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$taxa$species$abund        = as.matrix(read.table(sprintf('%s/results/tables/%s.species.%s.abund.tsv', project_path, project_name, tax_mode),
                                                         header=T, sep='\t', row.names=1, check.names=F))
                                                          
     
    SQM$taxa$superkingdom$percent = 100 * t(t(SQM$taxa$superkingdom$abund) / colSums(SQM$taxa$superkingdom$abund)) 
    SQM$taxa$phylum$percent       = 100 * t(t(SQM$taxa$phylum$abund)       / colSums(SQM$taxa$phylum$abund)) 
    SQM$taxa$class$percent        = 100 * t(t(SQM$taxa$class$abund)        / colSums(SQM$taxa$class$abund))  
    SQM$taxa$order$percent        = 100 * t(t(SQM$taxa$order$abund)        / colSums(SQM$taxa$order$abund))  
    SQM$taxa$family$percent       = 100 * t(t(SQM$taxa$family$abund)       / colSums(SQM$taxa$family$abund))
    SQM$taxa$genus$percent        = 100 * t(t(SQM$taxa$genus$abund)        / colSums(SQM$taxa$genus$abund))
    SQM$taxa$species$percent      = 100 * t(t(SQM$taxa$species$abund)      / colSums(SQM$taxa$species$abund))

    
    SQM$misc$long_names               = list()

    SQM$misc$long_names$superkingdom  = rownames(SQM$tax$superkingdom$abund)
    names(SQM$misc$long_names$superkingdom) = gsub('^k_', '', SQM$misc$long_names$superkingdom) # :'(

    SQM$misc$long_names$phylum        = rownames(SQM$tax$phylum$abund)
    names(SQM$misc$long_names$phylum) = sapply(strsplit(SQM$misc$long_names$phylum, split=';p_'), FUN = function(x) x[2])

    SQM$misc$long_names$class         = rownames(SQM$tax$class$abund)
    names(SQM$misc$long_names$class)  = sapply(strsplit(SQM$misc$long_names$class, split=';c_'), FUN = function(x) x[2])

    SQM$misc$long_names$order         = rownames(SQM$tax$order$abund)
    names(SQM$misc$long_names$order)  = sapply(strsplit(SQM$misc$long_names$order, split=';o_'), FUN = function(x) x[2])

    SQM$misc$long_names$family        = rownames(SQM$tax$family$abund)
    names(SQM$misc$long_names$family) = sapply(strsplit(SQM$misc$long_names$family, split=';f_'), FUN = function(x) x[2])

    SQM$misc$long_names$genus         = rownames(SQM$tax$genus$abund)
    names(SQM$misc$long_names$genus)  = sapply(strsplit(SQM$misc$long_names$genus, split=';g_'), FUN = function(x) x[2])

    SQM$misc$long_names$species       = rownames(SQM$tax$species$abund)
    names(SQM$misc$long_names$species)= sapply(strsplit(SQM$misc$long_names$species, split=';s_'), FUN = function(x) x[2]) 


    cat('Loading functions\n')
    SQM$functions                 = list()
    SQM$functions$KEGG            = list()
    SQM$functions$KEGG$abund      = as.matrix(read.table(sprintf('%s/results/tables/%s.KO.abund.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$functions$KEGG$tpm        = as.matrix(read.table(sprintf('%s/results/tables/%s.KO.tpm.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))
                                                             
    SQM$functions$COG             = list()
    SQM$functions$COG$abund       = as.matrix(read.table(sprintf('%s/results/tables/%s.COG.abund.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$functions$COG$tpm         = as.matrix(read.table(sprintf('%s/results/tables/%s.COG.tpm.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))

    SQM$functions$PFAM            = list()
    SQM$functions$PFAM$abund      = as.matrix(read.table(sprintf('%s/results/tables/%s.PFAM.abund.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))
    SQM$functions$PFAM$tpm        = as.matrix(read.table(sprintf('%s/results/tables/%s.PFAM.tpm.tsv', project_path, project_name),
                                                         header=T, sep='\t', row.names=1, check.names=F))

    cat('Loading total reads\n')
    SQM$total_reads               = as.matrix(
                                              read.table(sprintf('%s/results/10.%s.mappingstat', project_path, project_name), 
                                                         header=T, sep='\t', row.names=1, skip=1, comment.char='')
                                             )[,'Total.reads']

    return(SQM)

    }
