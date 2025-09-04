#' Get the N most abundant rows (or columns) from a numeric table
#'
#' Return a subset of an input matrix or data frame, containing only the N most abundant rows (or columns), sorted. Alternatively, a custom set of rows can be returned.
#' @param data numeric matrix or data frame
#' @param N integer Number of rows to return (default \code{10}).
#' @param items character vector. Custom row names to return. If provided, it will override \code{N} and \code{extra_items} (default \code{NULL}).
#' @param extra_items character vector. Extra row names to return on top of the N most abundant (default \code{NULL})
#' @param ignore character. Custom row names to drop before abundance calculation.
#' @param others logical. If \code{TRUE}, an extra row will be returned containing the aggregated abundances of the elements not selected with \code{N} or \code{items} (default \code{FALSE}).
#' @param rescale logical. Scale result to percentages column-wise (default \code{FALSE}).
#' @param bycol logical. Operate on columns instead of rows (default \code{FALSE}).
#' @return A matrix or data frame (same as input) with the selected rows (or columns).
#' @examples
#' data(Hadza)
#' Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
#' # Which are the 20 most abundant KEGG functions in the ORFs related to carbohydrate metabolism?
#' topCarb = mostAbundant(Hadza.carb$functions$KEGG$tpm, N=20)
#' # Now print them with nice names.
#' rownames(topCarb) = paste(rownames(topCarb),
#'                           Hadza.carb$misc$KEGG_names[rownames(topCarb)], sep="; ")
#' topCarb
#' # We can pass this to any R function.
#' heatmap(topCarb)
#' # But for convenience we provide wrappers for plotting ggplot2 heatmaps and barplots.
#' plotHeatmap(topCarb, label_y="TPM")
#' plotBars(topCarb, label_y="TPM")
#' @export
mostAbundant = function(data, N = 10, items = NULL, extra_items = NULL,
                        ignore = NULL, others = FALSE, rescale = FALSE, bycol = FALSE)
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
    type = typeof(data)

    if(bycol) { data = t(data) }

    if(!is.null(ignore)) # User ignores some rows
        {
        data = data[!rownames(data) %in% ignore,,drop=FALSE]
        }
    if(!is.null(items))  # User selects custom data.
        {
        # Check that items selection is possible and user have not asked for unknown things!
        if(bycol) { s = 'columns'; tgt = colnames(data) } else { s = 'rows'; tgt = rownames(data) }
        if(any(!items %in% tgt))
            {
            stop(sprintf('At least one of your custom items is not in the %s', s))
            }
    } else
        {
        extra_items = intersect(extra_items, rownames(data))
        total_items = nrow(data)
        if (N <= total_items) # Do we have at least N items?
            { # Do we have at least N items?
            if (N + length(extra_items) <= 0)
                {
                stop('N<=0 and no vector of items items was supplied. There is nothing to return')
                }
            } else
                { # User asks for sth impossible
                warning(sprintf('N=%s but only %s items exist. Returning %s items', N, total_items, total_items))
                N = total_items - length(extra_items)
                others = FALSE
            }
        data2 = data[!rownames(data) %in% extra_items,]
        items = names(sort(rowSums(data2), decreasing = TRUE)[1:N])
        items = c(items, extra_items)
        }
    other_items = colSums(data[!rownames(data) %in% items,, drop = FALSE])
    data = data[items, ,drop = FALSE]
    # Sum the abundances of the non-selected taxa
    if (others)
        {
        data = as.data.frame(rbind('Other' = other_items, data))
        }
    if (rescale)
        {
        data = 100 * t(t(data) / colSums(data))
        }
    if(type!='list')
        {
        data = as.matrix(data)
    }else
        {
        data = as.data.frame(data)
        }

    if(bycol) { data = t(data) }

    return(data)
    }

