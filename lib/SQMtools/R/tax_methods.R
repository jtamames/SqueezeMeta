get_tax_abund = function(SQM, tax_source, rank, count)
    {
    tk = get_tax_keys(SQM, tax_source, taxa_abunds = TRUE)
    tax_source = tk[1]
    tax_key = tk[2]
    taxa_abund = SQM[[tax_source]][[tax_key]]
    stopifnot(!is.null(taxa_abund[[rank]]))
    stopifnot(!is.null(taxa_abund[[rank]][[count]]))
    return(taxa_abund[[rank]][[count]])
    }

get_tax_table = function(SQM, tax_source)
    {
    tk = get_tax_keys(SQM, tax_source)
    tax_source = tk[1]
    tax_key = tk[2]
    return(SQM[[tax_source]][[tax_key]])
    }


get_tax_keys = function(SQM, tax_source, taxa_abunds = FALSE)
    {
    stopifnot(tax_source %in% c('orfs', 'contigs', 'bins', 'bins_sqm', 'bins_gtdb'))
    # For bins, return GTDB tax when available, unless asked explicitly for another taxonomy
    tax_key = 'tax'
    if(tax_source == 'bins' & !is.null(SQM$bins$tax_gtdb)) { tax_key = 'tax_gtdb'
    } else if (tax_source == 'bins_gtdb') { tax_source = 'bins'; tax_key = 'tax_gtdb' 
    } else if (tax_source == 'bins_sqm')  { tax_source = 'bins' }
    if(taxa_abunds) { tax_key = gsub('tax', 'tax_abund', tax_key) }
    if(is.null(SQM[[tax_source]][[tax_key]]))
        {
        stop(sprintf('Tax source "%s$%s" is not present in this object', tax_source, tax_key))
        }
    return(c(tax_source, tax_key))
    }


get_preferred_tax = function(SQM)
    {
    tk = get_tax_keys(SQM, SQM$misc$tax_source, taxa_abunds = TRUE)
    tax_source = tk[1]
    tax_key = tk[2]
    return(SQM[[tax_source]][[tax_key]])
    }


fix_gtdbtk_tax = function(tax_table)
    {
    for(bin in rownames(tax_table))
        {
        taxvec = tax_table[bin,]
        newtaxvec = c()
        last_tax = NA
        prefixes = c('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
        for(i in seq_along(prefixes))
            {
            tax = gsub(prefixes[i], '', taxvec[i])
            if(i==1) { last_tax = gsub('^Unclassified', '', tax) }
            if(nchar(tax) == 0)
                {
                tax = sprintf('Unclassified %s', last_tax)
            } else
                {
                last_tax = tax
                }
            newtaxvec = c(newtaxvec, tax)
            }
        tax_table[bin,] = newtaxvec
        }
    return(tax_table)
    }
