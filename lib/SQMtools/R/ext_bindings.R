#' Convert a SQM object into a microtable object from the \emph{microeco} package
#'
#' This function will convert the selected features from a SQM object into an object of the \code{\link[microeco]{microtable}} class from the \href{https://chiliubio.github.io/microeco/}{microeco} package. When possible, it will also include the taxonomy of the included features (for functional classifications, the taxonomy table will instead include the description of each feature ID). Optionally, it accepts a meta table that will be passed as provided to \code{microtable$new}.  
#'  
#'
#' @param SQM A SQM, SQMbunch or SQMlite object.
#' @param features character. Either \code{"orfs"}, \code{"contigs"}, \code{"bins"}, any taxonomic rank included in \code{SQM$taxa} or any functional classication included in \code{SQM$functions} (default \code{"tax"}). Note that a given feature type might not be available in this objects (e.g. \code{"contigs"} in SQMlite objects originating from a SQM reads project).
#' @param count character. Either \code{"abund"} for raw abundances, \code{"percent"} for percentages, \code{"bases"} for raw base counts, \code{"cov"} for coverages, \code{"cpm"} for coverages per million reads, \code{"tpm"} for TPM normalized values or \code{"copy_number"} for copy numbers (default \code{"abund"}). Note that a given count type might not available in this object (e.g. TPM or copy number in SQMlite objects originating from a SQM reads project).
#' @param md data.frame. A optional data.frame containing metadata for the samples in the SQM object.
#' @param nocds character. Either \code{"treat_separately"} to treat features annotated as No CDS separately, \code{"treat_as_unclassified"} to treat them as Unclassified or \code{"ignore"} to ignore them in the output (default \code{"treat_separately"}).
#' @param no_partial_classifications logical. When \code{features} is a taxonomic rank, treat features not fully classified at the requested level (e.g. "Unclassified bacteroidota" at the class level or below) as fully unclassified. This takes effect before \code{ignore_unclassified}, so if both are \code{TRUE} the plot will only contain features that were fully classified at the requested level (default \code{FALSE}).
#' @param ignore_unclassified logical. When \code{features} is a taxonomic rank or functional category, don't include unclassified reads in the output (default \code{FALSE}).
#' @param ignore_unmapped logical. Don't include unmapped reads in the output (default \code{FALSE}).
#' @param bin_tax_source character. Source of taxonomy when \code{features = "bins"}, either \code{"SQM"} of \code{"gtdb"} (default  \code{"gtdb"}).
#' @param include_seqs logical. Whether to include sequences or not if creating a microtable from contigs (default \code{FALSE}).
#' @return A \code{\link[microeco]{microtable}}. 
#' @seealso \code{\link{SQM_to_phyloseq}} for exporting a SQM/SQMlite/SQM object as a phyloseq object.
#' @export
SQM_to_microeco = function(SQM, features = 'genus', count = 'abund', md = NULL,
                           nocds = 'treat_separately', no_partial_classifications = FALSE,
                           ignore_unclassified = FALSE, ignore_unmapped = FALSE,
                           bin_tax_source = 'SQM', include_seqs = FALSE)
    {
    if(!requireNamespace("microeco")) { stop('The `microeco` package must be installed in order to use this function') }
    # We let argument checking be done by prepare_export_tables
    res = prepare_export_tables(SQM, features, count,
                                nocds, no_partial_classifications, ignore_unclassified, ignore_unmapped,
                                bin_tax_source, include_seqs)
    if(!is.null(res[['tax']])) { res[['tax']] = as.data.frame(res[['tax']]) }
    if(features=='orfs') { res[['seqs']] = NULL } # since microtable only accepts DNAStringSets
    mt = microeco::microtable$new(otu_table = as.data.frame(res[['counts']]),
		                  tax_table = res[['tax']],
		                  rep_fasta = res[['seqs']],
                                  sample_table = md)
    return(mt)
    }


#' Convert a SQM object into a phyloseq object from the \emph{phyloseq} package
#'
#' This function will convert the selected features from a SQM object into a phyloseq object from the \href{https://joey711.github.io/phyloseq/}{phyloseq} package. When possible, it will also include the taxonomy of the included features  (for functional classifications, the taxonomy table will instead include the description of each feature ID). Optionally, it accepts a meta table that will be passed as provided to the \code{phyloseq} object constructor.
#'
#'
#' @param SQM A SQM, SQMbunch or SQMlite object.
#' @param features character. Either \code{"orfs"}, \code{"contigs"}, \code{"bins"}, any taxonomic rank included in \code{SQM$taxa} or any functional classication included in \code{SQM$functions} (default \code{"tax"}). Note that a given feature type might not be available in this objects (e.g. \code{"contigs"} in SQMlite objects originating from a SQM reads project).
#' @param count character. Either \code{"abund"} for raw abundances, \code{"percent"} for percentages, \code{"bases"} for raw base counts, \code{"cov"} for coverages, \code{"cpm"} for coverages per million reads, \code{"tpm"} for TPM normalized values or \code{"copy_number"} for copy numbers (default \code{"abund"}). Note that a given count type might not available in this object (e.g. TPM or copy number in SQMlite objects originating from a SQM reads project).
#' @param md data.frame. A optional data.frame containing metadata for the samples in the SQM object.
#' @param nocds character. Either \code{"treat_separately"} to treat features annotated as No CDS separately, \code{"treat_as_unclassified"} to treat them as Unclassified or \code{"ignore"} to ignore them in the output (default \code{"treat_separately"}).
#' @param no_partial_classifications logical. When \code{features} is a taxonomic rank, treat features not fully classified at the requested level (e.g. "Unclassified bacteroidota" at the class level or below) as fully unclassified. This takes effect before \code{ignore_unclassified}, so if both are \code{TRUE} the plot will only contain features that were fully classified at the requested level (default \code{FALSE}).
#' @param ignore_unclassified logical. When \code{features} is a taxonomic rank or functional category, don't include unclassified reads in the output (default \code{FALSE}).
#' @param ignore_unmapped logical. Don't include unmapped reads in the output (default \code{FALSE}).
#' @param bin_tax_source character. Source of taxonomy when \code{features = "bins"}, either \code{"SQM"} of \code{"gtdb"} (default  \code{"gtdb"}).
#' @param include_seqs logical. Whether to include sequences or not if creating a microtable from ORFs or contigs (default \code{FALSE}).
#' @return A phyloseq object.
#' @seealso \code{\link{SQM_to_microeco}} for exporting a SQM/SQMlite/SQM object as a microtable object.
#' @export
SQM_to_phyloseq = function(SQM, features = 'genus', count = 'abund', md = NULL,                                                                                 nocds = 'treat_separately', no_partial_classifications = FALSE,
                           ignore_unclassified = FALSE, ignore_unmapped = FALSE,
                           bin_tax_source = 'SQM', include_seqs = FALSE)
    {
    if(!requireNamespace("phyloseq")) { stop('The `phyloseq` package must be installed in order to use this function') }
    # We let argument checking be done by prepare_export_tables
    res = prepare_export_tables(SQM, features, count,
                                nocds, no_partial_classifications, ignore_unclassified, ignore_unmapped,
                                bin_tax_source, include_seqs)
    ps = phyloseq::phyloseq(phyloseq::otu_table(res[['counts']], taxa_are_rows = TRUE))
    if(!is.null(res[['tax' ]]))
        {
        tt = phyloseq::tax_table(as.matrix(res[['tax']]))
        rownames(tt) = rownames(res[['tax']]) # because phyloseq is dumb
        ps = phyloseq::merge_phyloseq(ps, tt)
        }
    if(!is.null(res[['seqs']])) { ps = phyloseq::merge_phyloseq(ps, res[['seqs']]          ) }
    if(!is.null(md)           ) { ps = phyloseq::merge_phyloseq(ps, phyloseq::sample_data(md)        ) }
    return(ps)
    }

prepare_export_tables = function(SQM, features, count,
                                 nocds, no_partial_classifications, ignore_unclassified, ignore_unmapped,
                                 bin_tax_source, include_seqs)
    {
    if(!inherits(SQM, c('SQM', 'SQMlite', 'SQMbunch'))) { stop('The first argument must be a SQM or SQMbunch object') }
    if(!bin_tax_source %in% c('SQM', 'gtdb')) { stop('bin_tax_source mush be either "SQM" or "gtdb"') }
    if(!nocds %in% c('treat_separately', 'treat_as_unclassified', 'ignore'))
        {
        stop('nocds must be "treat_separately", "treat_as_unclassified" or "ignore"')
        }
    if(features %in% c('orfs', 'contigs', 'bins'))
        {
        if(!inherits(SQM, c('SQM')))
            {
            stop(sprintf('features="%s" can only be requested with SQM objects, not with SQMlite or SQMbunch objects', features))
            }
        if(features == 'bins' & is.null(SQM$bins))
            {
            stop('This SQM object contains no binning results')
            }
        counts = get_counts(SQM[[features]], features, count)
        tax = SQM[[features]][['tax']]
        if(features == 'bins' & bin_tax_source == "gtdb")
            {
            if('tax_gtdb' %in% names(SQM[[features]]))
                {
                tax = SQM[[features]][['tax_gtdb']]
            } else
                {
                warning('GTDB taxonomy is not present in this SQM object, falling back to SQM taxonomy for bins')
                }
            } 
        seqs = NULL
        if(include_seqs & features %in% c('contigs', 'orfs'))
            {
            seqs = SQM[[features]][['seqs']]
            }
    } else if(features %in% names(SQM$taxa))
        {
        counts = get_counts(SQM$taxa[[features]], features, count)
        tax = as.data.frame(t(data.frame(strsplit(SQM$misc$tax_names_long[[features]], ';'), check.names = FALSE)))
        seqs = NULL
    } else if(features %in% names(SQM$functions))
        {
        counts = get_counts(SQM$functions[[features]], features, count)
        tax = NULL
        feature_names = paste0(features, "_names")
        feature_paths = paste0(features, "_paths")
        if(feature_names %in% names(SQM$misc))
            {
            tax = data.frame(SQM$misc[[feature_names]])[rownames(counts),,drop=FALSE]
            colnames(tax) = features
            if(feature_paths %in% names(SQM$misc))
                {
                tax[,feature_paths] = SQM$misc[[feature_paths]][rownames(counts)]
                tax = tax[,c(feature_paths, features),drop=FALSE]
                }
            }
        seqs = NULL
    } else
        {
        stop('features must be "orfs", "contigs", "bins", any taxonomic rank included in SQM$taxa or any functional classication included in SQM$functions')
        }

    # Handle No CDS
    if(nocds == 'treat_as_unclassified' & 'No CDS' %in% rownames(counts))
        {
        if(!'Unclassified' %in% rownames(counts))
            {
            # Create from scratch
            counts = rbind(counts, Unclassified = counts['No CDS',,drop=FALSE])
        } else
            {
            # Or add Unclassified and No CDS
            counts['Unclassified',] = counts['Unclassified',,drop=FALSE] + counts['No CDS',,drop=FALSE]
            }
        }
    if(nocds %in% c('ignore', 'treat_as_unclassified')) # we remove them from the table in both cases
        {
        counts = counts[rownames(counts)!='No CDS',,drop=FALSE]
        }

    # Handle partial classifications
    if(no_partial_classifications & features %in% names(SQM$taxa))
        {
        partials = grepl('^Unclassified', rownames(counts)) & rownames(counts) != 'Unclassified'
        partials_abund = colSums(counts[partials,,drop=FALSE])
        if(!'Unclassified' %in% rownames(counts))
            {
            counts = rbind(counts, Unclassified = partials_abund)
        } else 
            {
            counts['Unclassified',] = counts['Unclassified',,drop=FALSE] + partials_abund
            }
        counts = counts[!partials,,drop=FALSE]
        }

    # Remove Unclassified if needed
    if(ignore_unclassified)
        {
        counts = counts[rownames(counts) != 'Unclassified',,drop=FALSE]
        }    

    # Remove Unmapped if needed
    if(ignore_unmapped)
        {
        counts = counts[rownames(counts)!='Unmapped',,drop=FALSE]
        }

    # Fix tax table
    if(!is.null(tax))
        {
        if(features %in% c('bins', names(SQM$taxa)))
            {
            # Change rank names to fit microeco defaults
            colnames(tax) = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')[1:ncol(tax)]
            }
        # Add extra rows for things like "Unmapped" and "No_bin"
        extra_rows = setdiff(rownames(counts), rownames(tax))
        for(rn in extra_rows)
            {
            tax = rbind(tax, matrix(rn, nrow=1, ncol=ncol(tax), dimnames = list(c(rn), colnames(tax))))
            }
        tax = tax[rownames(counts),,drop=FALSE]
        }

    return( list(counts = counts, tax = tax, seqs = seqs) )
    }


get_counts = function(SQMlist, features, count)
    {
    good_counts = intersect(names(SQMlist), c('abund', 'bases', 'percent', 'cov', 'cpm', 'tpm', 'copy_number'))
    if(!count %in% good_counts)
        {
        stop(sprintf('count="%s" is not available for %s in this object. Available counts are: %s',
                     count, features, paste0(good_counts, collapse = ', ')
                    )
            )
        }
    counts = SQMlist[[count]]
    if(count %in% c('abund', 'bases'))
        {
        # Some may not be integers (if dealing with multi-function KEGG/COG), so force rounding
        counts = round(counts)
        }
    return(counts)
    }
