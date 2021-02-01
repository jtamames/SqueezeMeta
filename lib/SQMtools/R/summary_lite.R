#' summary method for class SQMlite
#'
#' Computes different statistics of the data contained in the SQMlite object.
#' @return A list of summary statistics.
#' @export
summary.SQMlite = function(SQM)
    {
    if(!class(SQM)=='SQMlite') { stop('The first argument must be a SQMlite object') }

    res = list()
    res$project_name = SQM$misc$project_name

    res$nReads = colSums(SQM$taxa$superkingdom$abund)
    res$nORFs  = colSums(SQM$functions[[1]]$abund)
    res$taxa   = list()

    sk = rownames(SQM$taxa$superkingdom$abund)
    sk = sk[!grepl('^Unclassified', sk)] # Remove Unclassified
    sk = sk[!grepl('in NCBI)', sk)]      # Remove "virtual" taxa (coming for ranks missing in NCBI for some taxa).
    p  = rownames(SQM$taxa$phylum$abund)
    p  = p[!grepl('^Unclassified', p)]
    p  = p[!grepl('in NCBI)', p)]
    c_ = rownames(SQM$taxa$class$abund)
    c_ = c_[!grepl('^Unclassified', c_)]
    c_ = c_[!grepl('in NCBI)', c_)]
    o  = rownames(SQM$taxa$order$abund)
    o  = o[!grepl('^Unclassified', o)]
    o  = o[!grepl('in NCBI)', o)]
    f  = rownames(SQM$taxa$family$abund)
    f  = f[!grepl('^Unclassified', f)]
    f  = f[!grepl('in NCBI)', f)]
    g  = rownames(SQM$taxa$genus$abund)
    g  = g[!grepl('^Unclassified', g)]
    g  = g[!grepl('in NCBI)', g)]
    s  = rownames(SQM$taxa$species$abund)
    s  = s[!grepl('^Unclassified', s)]
    s  = s[!grepl('in NCBI)', s)]

    res$taxa$superkingdom = list(
                                 nReads = colSums(SQM$taxa$superkingdom$abund[sk,]), 
                                 nTaxa  = colSums(SQM$taxa$superkingdom$abund[sk,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$superkingdom$abund)
                                )
    res$taxa$phylum       = list(
                                 nReads = colSums(SQM$taxa$phylum$abund[p,]),
                                 nTaxa  = colSums(SQM$taxa$phylum$abund[p,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$phylum$abund)
                                )
    res$taxa$class        = list(
                                 nReads = colSums(SQM$taxa$class$abund[c_,]),
                                 nTaxa  = colSums(SQM$taxa$class$abund[c_,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$class$abund)
                                )
    res$taxa$order        = list(
                                 nReads = colSums(SQM$taxa$order$abund[o,]),
                                 nTaxa  = colSums(SQM$taxa$order$abund[o,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$order$abund)
                                )
    res$taxa$family       = list(
                                 nReads = colSums(SQM$taxa$family$abund[f,]),
                                 nTaxa  = colSums(SQM$taxa$family$abund[f,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$family$abund)
                                )
    res$taxa$genus        = list(
                                 nReads = colSums(SQM$taxa$genus$abund[g,]),
                                 nTaxa  = colSums(SQM$taxa$genus$abund[g,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$genus$abund)
                                )
    res$taxa$species      = list(
                                 nReads = colSums(SQM$taxa$species$abund[s,]),
                                 nTaxa  = colSums(SQM$taxa$species$abund[s,]>0),
                                 most_abundant = most_abundant_row(SQM$taxa$species$abund)
                                )


    res$functions         = list()
    res$functions$KEGG    = colSums( SQM$functions$KEGG$abund[rownames(SQM$functions$KEGG$abund) != 'Unclassified',] )
    res$functions$COG     = colSums( SQM$functions$COG$abund [rownames(SQM$functions$COG$abund ) != 'Unclassified',] )
    for(method in SQM$misc$ext_annot_sources)
        {
        res$functions[[method]]    = colSums( SQM$functions[[method]]$abund[rownames(SQM$functions[[method]]$abund) != 'Unclassified',] )
        }

    res$samples           = SQM$misc$samples
    res$ext_annot_sources = SQM$misc$ext_annot_sources

    class(res) = 'summary.SQMlite'

    return(res)
    }


#' @export
#' @noRd
print.summary.SQMlite = function(summ)
    {
    cat('\n')
    cat( sprintf('\tBASE PROJECT NAME: %s\n', summ$project_name) )
    cat('\n')
    cat( sprintf('\t\t%s\n'            , paste(summ$samples                        , collapse='\t')) )
    cat( sprintf('\tTOTAL READS\t%s\n' , paste(summ$nReads                         , collapse='\t')) )    
    cat( sprintf('\tTOTAL ORFs\t%s\n' , paste(summ$nORFs                           , collapse='\t')) )
    cat('\n\t----------------------------------------------------------\n\n')
    cat('\tTAXONOMY:\n\n')
    cat('\tClassified reads:\n')
    cat( sprintf('\t\t%s\n'            , paste(summ$samples                        , collapse='\t')) )
    cat( sprintf('\tSuperkingdom\t%s\n', paste(summ$taxa$superkingdom$nReads       , collapse='\t')) )
    cat( sprintf('\tPhylum\t%s\n'      , paste(summ$taxa$phylum$nReads             , collapse='\t')) )
    cat( sprintf('\tClass\t%s\n'       , paste(summ$taxa$class$nReads              , collapse='\t')) )
    cat( sprintf('\tOrder\t%s\n'       , paste(summ$taxa$order$nReads              , collapse='\t')) )
    cat( sprintf('\tFamily\t%s\n'      , paste(summ$taxa$family$nReads             , collapse='\t')) )
    cat( sprintf('\tGenus\t%s\n'       , paste(summ$taxa$genus$nReads              , collapse='\t')) )
    cat( sprintf('\tSpecies\t%s\n'     , paste(summ$taxa$species$nReads            , collapse='\t')) )
    cat('\n')
    cat('\tMost abundant taxa (ignoring Unclassified):\n')
    cat( sprintf('\t\t%s\n'            , paste(summ$samples                        , collapse='\t')) )
    cat( sprintf('\tSuperkingdom\t%s\n', paste(summ$taxa$superkingdom$most_abundant, collapse='\t')) )
    cat( sprintf('\tPhylum\t%s\n'      , paste(summ$taxa$phylum$most_abundant      , collapse='\t')) )
    cat( sprintf('\tClass\t%s\n'       , paste(summ$taxa$class$most_abundant       , collapse='\t')) )
    cat( sprintf('\tOrder\t%s\n'       , paste(summ$taxa$order$most_abundant       , collapse='\t')) )
    cat( sprintf('\tFamily\t%s\n'      , paste(summ$taxa$family$most_abundant      , collapse='\t')) )
    cat( sprintf('\tGenus\t%s\n'       , paste(summ$taxa$genus$most_abundant       , collapse='\t')) )
    cat( sprintf('\tSpecies\t%s\n'     , paste(summ$taxa$species$most_abundant     , collapse='\t')) )
    cat('\n\t----------------------------------------------------------\n\n')
    cat('\tFUNCTIONS:\n')
    cat('\tClassified ORFs:\n')
    cat( sprintf('\t\t%s\n'            , paste(summ$samples                        , collapse='\t')) )
    for(method in names(summ$functions))
        {
        cat( sprintf('\t%s\t%s\n'      , method, paste(summ$functions[[method]]    , collapse='\t')) )
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
