#' Get the N most variable rows (or columns) from a numeric table
#'
#' Return a subset of an input matrix or data frame, containing only the N most variable rows (or columns), sorted. Variability is calculated as the Coefficient of Variation (sd/mean).
#' @param data numeric matrix or data frame
#' @param N integer Number of rows to return (default \code{10}).
#' @param bycol logical. Operate on columns instead of rows (default \code{FALSE}).
#' @return A matrix or data frame (same as input) with the selected rows or columns.
#' @examples
#' data(Hadza)
#' Hadza.carb = subsetFun(Hadza, "Carbohydrate metabolism")
#' # Which are the 20 most variable KEGG functions in the ORFs related to carbohydrate metabolism?
#' topCarb = mostVariable(Hadza.carb$functions$KEGG$tpm, N=20)
#' # Now print them with nice names
#' rownames(topCarb) = paste(rownames(topCarb),
#'                           Hadza.carb$misc$KEGG_names[rownames(topCarb)], sep="; ")
#' topCarb
#' # We can pass this to any R function
#' heatmap(topCarb)
#' # But for convenience we provide wrappers for plotting ggplot2 heatmaps and barplots
#' plotHeatmap(topCarb, label_y="TPM")
#' plotBars(topCarb, label_y="TPM")
#' @importFrom stats sd
#' @export
mostVariable = function(data, N = 10, bycol = FALSE)
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }

    if(bycol) { data = t(data) }

    total_items = nrow(data)

    if (N <= total_items) # Do we have at least N items?
        { # Do we have at least N items?
        if (N <= 0)
            {
            stop('N<=0. There is nothing to return')
            }
        } else
            { # User asks for sth impossible
            warning(sprintf('N=%s but only %s items exist. Returning %s items', N, total_items, total_items))
            N = total_items
        }
        items = names(sort(apply(data,1,sd)/apply(data,1,mean), decreasing = TRUE)[1:N])
        

    data = data[items,,drop=FALSE]

    if(bycol) { data = t(data) }

    return(data)
    }

