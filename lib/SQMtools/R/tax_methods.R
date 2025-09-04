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
