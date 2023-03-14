#' summary method for class SQM
#'
#' Computes different statistics of the data contained in the SQM object.
#' @param object SQM object to be summarized.
#' @param ... Additional parameters (ignored).
#' @return A list of summary statistics.
#' @export
summary.SQM = function(object, ...)
    {

    SQM = object # so that CRAN is happy

    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }

    res = list()

    res$project_name = SQM$misc$project_name

    res$reads_in_orfs        = colSums(SQM$orfs$abund)
    res$original_reads       = SQM$total_reads
    res$percent_in_orfs      = 100 * res$reads_in_orfs / res$original_reads

    res$contigs = list()
    res$contigs$nContigs     = nrow(SQM$contigs$table)
    res$contigs$total_length = sum(as.numeric(SQM$contigs$table$Length)) # as.numeric to avoid integer overflow.
    res$contigs$longest      = max(SQM$contigs$table$Length)
    res$contigs$shortest     = min(SQM$contigs$table$Length)
    res$contigs$N50          = Npercent(as.numeric(SQM$contigs$table$Length), 50) # as.numeric to avoid integer overflow.
    res$contigs$N90          = Npercent(as.numeric(SQM$contigs$table$Length), 90) # as.numeric to avoid integer overflow.

    sk = SQM$contigs$tax[,'superkingdom']
    sk = sk[!grepl('^Unclassified|^Unmapped|^No CDS', sk)] # Remove Unclassified
    sk = sk[!grepl('in NCBI)', sk)]      # Remove "virtual" taxa (coming for ranks missing in NCBI for some taxa).
    p  = SQM$contigs$tax[,'phylum']
    p  = p[!grepl('^Unclassified|^Unmapped|^No CDS', p)]
    p  = p[!grepl('in NCBI)', p)]
    c_ = SQM$contigs$tax[,'class']
    c_ = c_[!grepl('^Unclassified|^Unmapped|^No CDS', c_)]
    c_ = c_[!grepl('in NCBI)', c_)]
    o  = SQM$contigs$tax[,'order']
    o  = o[!grepl('^Unclassified|^Unmapped|^No CDS', o)]
    o  = o[!grepl('in NCBI)', o)]
    f  = SQM$contigs$tax[,'family']
    f  = f[!grepl('^Unclassified|^Unmapped|^No CDS', f)]
    f  = f[!grepl('in NCBI)', f)]
    g  = SQM$contigs$tax[,'genus']
    g  = g[!grepl('^Unclassified|^Unmapped|^No CDS', g)]
    g  = g[!grepl('in NCBI)', g)]
    s  = SQM$contigs$tax[,'species']
    s  = s[!grepl('^Unclassified|^Unmapped|^No CDS', s)]
    s  = s[!grepl('in NCBI)', s)]

    res$contigs$superkingdom = list(nContigs = length(sk), nTaxa = length(unique(sk)), most_abundant = most_abundant_row(SQM$tax$superkingdom$abund) )
    res$contigs$phylum       = list(nContigs = length(p ), nTaxa = length(unique(p )), most_abundant = most_abundant_row(SQM$tax$phylum$abund      ) )
    res$contigs$class        = list(nContigs = length(c_), nTaxa = length(unique(c_)), most_abundant = most_abundant_row(SQM$tax$class$abund       ) )
    res$contigs$order        = list(nContigs = length(o ), nTaxa = length(unique(o )), most_abundant = most_abundant_row(SQM$tax$order$abund       ) )
    res$contigs$family       = list(nContigs = length(f ), nTaxa = length(unique(f )), most_abundant = most_abundant_row(SQM$tax$family$abund      ) )
    res$contigs$genus        = list(nContigs = length(g ), nTaxa = length(unique(g )), most_abundant = most_abundant_row(SQM$tax$genus$abund       ) )
    res$contigs$species      = list(nContigs = length(s ), nTaxa = length(unique(s )), most_abundant = most_abundant_row(SQM$tax$species$abund     ) )

    res$orfs = list()
    res$orfs$total     = colSums(SQM$orfs$abund>0)
    res$orfs$KEGG      = colSums( (SQM$orfs$abund>0) * (SQM$orfs$table[,'KEGG ID'] != '') )
    res$orfs$KEGG_good = colSums( (SQM$orfs$abund>0) * grepl('*', SQM$orfs$table[,'KEGG ID'], fixed=T) )
    res$orfs$COG       = colSums( (SQM$orfs$abund>0) * (SQM$orfs$table[,'COG ID' ] != '') )
    res$orfs$COG_good  = colSums( (SQM$orfs$abund>0) * grepl('*', SQM$orfs$table[,'COG ID' ], fixed=T) )
    res$orfs$PFAM      = colSums( (SQM$orfs$abund>0) * (SQM$orfs$table[,'PFAM'   ] != '') )
    for(method in SQM$misc$ext_annot_sources)
        {
        res$orfs[[method]]                     = colSums( (SQM$orfs$abund>0) * (SQM$orfs$table[,method] != '') )
        res$orfs[[sprintf('%s_good', method)]] = colSums( (SQM$orfs$abund>0) * grepl('*', SQM$orfs$table[,method], fixed=T) )
        }

    res$samples           = SQM$misc$samples
    res$ext_annot_sources = SQM$misc$ext_annot_sources

    class(res) = 'summary.SQM'

    return(res)
    }


#' @export
#' @noRd
print.summary.SQM = function(x, ...)
    {
    summ = x # so that CRAN is happy
    cat('\n')
    cat( sprintf('\tBASE PROJECT NAME: %s\n', summ$project_name) )
    cat('\n\t----------------------------------------------------------\n\n')
    cat('\tREADS:\n\n')
    cat( sprintf('\t\t%s\n'               , paste(summ$samples                  , collapse='\t')) )
    cat( sprintf('\tMapping to ORFs\t%s\n', paste(summ$reads_in_orfs            , collapse='\t')) )
    cat( sprintf('\tInput reads\t%s\n'    , paste(summ$original_reads           , collapse='\t')) )
    cat( sprintf('\tPercent\t%s\n'        , paste(round(summ$percent_in_orfs, 1), collapse='\t')) )
    cat('\n\t----------------------------------------------------------\n\n')
    cat('\tCONTIGS:\n')
    cat( sprintf('\tNumber of contigs:\t%s\n'    , summ$contigs$nContigs    ) ) 
    cat( sprintf('\tTotal length:\t%s bases\n'   , summ$contigs$total_length) )
    cat( sprintf('\tLongest contig:\t%s bases\n' , summ$contigs$longest     ) )
    cat( sprintf('\tShortest contig:\t%s bases\n', summ$contigs$shortest    ) )
    cat( sprintf('\tN50:\t%s bases\n'            , summ$contigs$N50         ) )
    cat( sprintf('\tN90:\t%s bases\n'            , summ$contigs$N90         ) )
    cat('\n')
    cat( sprintf('\tSuperkingdom:\t %s (%s%%) classified contigs in %s taxa\n',
                 summ$contigs$superkingdom$nContigs, round(100*summ$contigs$superkingdom$nContigs/summ$contigs$nContigs, 1), summ$contigs$superkingdom$nTaxa ))
    cat( sprintf('\tPhylum:\t %s (%s%%) classified contigs in %s taxa\n'      ,
                 summ$contigs$phylum$nContigs      , round(100*summ$contigs$phylum$nContigs/summ$contigs$nContigs,       1), summ$contigs$phylum$nTaxa       ))
    cat( sprintf('\tClass:\t %s (%s%%) classified contigs in %s taxa\n'       ,
                 summ$contigs$class$nContigs       , round(100*summ$contigs$class$nContigs/summ$contigs$nContigs,        1), summ$contigs$class$nTaxa        ))
    cat( sprintf('\tOrder:\t %s (%s%%) classified contigs in %s taxa\n'       , 
                 summ$contigs$order$nContigs       , round(100*summ$contigs$order$nContigs/summ$contigs$nContigs,        1), summ$contigs$order$nTaxa        ))
    cat( sprintf('\tFamily:\t %s (%s%%) classified contigs in %s taxa\n'      ,
                 summ$contigs$family$nContigs      , round(100*summ$contigs$family$nContigs/summ$contigs$nContigs,       1), summ$contigs$family$nTaxa       ))
    cat( sprintf('\tGenus:\t %s (%s%%) classified contigs in %s taxa\n'       ,
                 summ$contigs$genus$nContigs       , round(100*summ$contigs$genus$nContigs/summ$contigs$nContigs,        1), summ$contigs$genus$nTaxa        ))
    cat( sprintf('\tSpecies:\t %s (%s%%) classified contigs in %s taxa\n'     ,
                 summ$contigs$species$nContigs     , round(100*summ$contigs$species$nContigs/summ$contigs$nContigs,      1), summ$contigs$species$nTaxa      ))
    cat('\n')
    cat('\tMost abundant taxa (ignoring Unclassified):\n')
    cat( sprintf('\t\t%s\n'            , paste(summ$samples                           , collapse='\t')) )
    cat( sprintf('\tSuperkingdom\t%s\n', paste(summ$contigs$superkingdom$most_abundant, collapse='\t')) )
    cat( sprintf('\tPhylum\t%s\n'      , paste(summ$contigs$phylum$most_abundant      , collapse='\t')) )
    cat( sprintf('\tClass\t%s\n'       , paste(summ$contigs$class$most_abundant       , collapse='\t')) )
    cat( sprintf('\tOrder\t%s\n'       , paste(summ$contigs$order$most_abundant       , collapse='\t')) )
    cat( sprintf('\tFamily\t%s\n'      , paste(summ$contigs$family$most_abundant      , collapse='\t')) )
    cat( sprintf('\tGenus\t%s\n'       , paste(summ$contigs$genus$most_abundant       , collapse='\t')) )
    cat( sprintf('\tSpecies\t%s\n'     , paste(summ$contigs$species$most_abundant     , collapse='\t')) )
    cat('\n\t----------------------------------------------------------\n\n')
    cat('\tORFs:\n')
    cat( sprintf('\t\t%s\n'                     , paste(summ$samples       , collapse='\t')) )
    cat( sprintf('\tTotal\t%s\n'                , paste(summ$orfs$total    , collapse='\t')) )
    cat( sprintf('\tWith KEGG\t%s\n'            , paste(summ$orfs$KEGG     , collapse='\t')) )
    cat( sprintf('\tWith KEGG (best aver)\t%s\n', paste(summ$orfs$KEGG_good, collapse='\t')) )
    cat( sprintf('\tWith COG\t%s\n'             , paste(summ$orfs$COG      , collapse='\t')) )
    cat( sprintf('\tWith COG (best aver)\t%s\n' , paste(summ$orfs$COG_good , collapse='\t')) )
    cat( sprintf('\tWith PFAM\t%s\n'            , paste(summ$orfs$PFAM     , collapse='\t')) )
    for(method in summ$ext_annot_sources)
        {
        cat( sprintf('\tWith %s\t%s\n'            , method, paste(summ$orfs[[method]]                      , collapse='\t')) )
        cat( sprintf('\tWith %s (best aver)\t%s\n', method, paste(summ$orfs[[sprintf('%s_good', method)]] , collapse='\t')) )
        }
    cat('\n')
    }



Npercent  = function(len, percent)
    {
    # Taken from https://gist.github.com/shujishigenobu/1858458
    len_sorted = rev(sort(len))
    Npercent = len_sorted[cumsum(len_sorted) >= sum(len_sorted)*percent/100][1]
    return(Npercent)
    }


#' Return the name of the most abundant row for all the colums of a taxa.
#' @noRd
most_abundant_row = function(table, ignore_unclassified=T)
    {
    if(ignore_unclassified) { table = table[!grepl('Unclassified|Unmapped|^No CDS', rownames(table)),,drop=F] }
    colMaxsIdx = apply(table, 2, which.max)
    res = rownames(table)[colMaxsIdx]
    if(is.null(res)) { res = rep('Unclassified', ncol(table)) }
    names(res) = colnames(table)
    return(res)
    }
