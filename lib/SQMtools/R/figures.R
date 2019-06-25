library(ggplot2)
library(reshape2)

#' Plot a heatmap using ggplot2
#'
#' Plot a ggplot2 heatmap from a matrix or data frame. The data should be in tabular format (e.g. features in rows and samples in columns).
#' @param data numeric matrix or data frame.
#' @param label_x character Label for the x axis (default \code{"Samples"}).
#' @param label_y character Label for the y axis (default \code{"Features"}).
#' @param label_fill character Label for color scale (default \code{"Abundance"}).
#' @param gradient_col A vector of two colors representing the low and high ends of the color gradient (default \code{c("ghostwhite", "dodgerblue4")}).
#' @return A ggplot2 plot object.
#' @seealso \code{\link[plotFunctions]{plotFunctions}} for plotting the top functional categories of a SQM object; \code{\link[plotBars]{plotBars}} for plotting a barplot with arbitrary data; \code{\link[mostAbundant]{mostAbundant}} for selecting the most abundant rows in a dataframe or matrix.
#' @examples
#' data(Hadza)
#' topPFAM = mostAbundant(Hadza$functions$PFAM$tpm)
#' topPFAM = topPFAM[rownames(topPFAM) != 'Unclassified',] # Take out the Unclassified ORFs.
#' plotHeatmap(topPFAM, label_x = 'Samples', label_y = 'PFAMs', label_fill = 'TPM')
#' @export
#' @examples
#' data(Hadza)
#' phyla_percent = Hadza$taxa$phylum$percent
#' plotHeatmap(phyla_percent, label_y = 'Phylum', label_fill = 'Percentage')
plotHeatmap = function(data, label_x = 'Samples', label_y = 'Features', label_fill = 'Abundance', gradient_col = c('ghostwhite', 'dodgerblue4'))
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
    if(length(gradient_col) < 2)
        {
        stop('gradient_col must be a vector with two colors representing the low and high ends of the color gradient')
        }
    data=t(data)
    # data = data[, order(colSums(data), decreasing = T), drop = F] # Order functions according to their abundances
    data_melt = melt(as.matrix(data), value.name = 'abun')
    colnames(data_melt) = c('sample', 'item', 'abun')
    data_melt$abun = as.numeric(data_melt$abun)
    #PLOT DATA
    if(is.null(label_x)) { label_x = '' }
    if(is.null(label_y)) { label_y = '' }
    if(is.na(label_x)  ) { label_x = '' }
    if(is.na(label_y)  ) { label_y = '' }

    if(nchar(label_x)==0) { theme_x = element_blank()
    }else{ theme_x = element_text(size = 12) }
    if(nchar(label_y)==0) { theme_y = element_blank()
    }else{ theme_y = element_text(size = 12) }
        
    p = ggplot(data_melt, aes(x = sample, y = item, fill = abun))
    p = p + geom_tile()
    p = p + scale_fill_gradient(low = gradient_col[1], high = gradient_col[2])
    p = p + theme_light()
    p = p + theme(axis.title.x = theme_x, axis.title.y = theme_y)
    p = p + labs(x = label_x, y = label_y, fill = label_fill)
    return(p)
    }



#' Plot a barplot using ggplot2
#'
#' Plot a ggplot2 barplot from a matrix or data frame. The data should be in tabular format (e.g. features in rows and samples in columns).
#' @param data Numeric matrix or data frame.
#' @param label_x character Label for the x axis (default \code{"Samples"}).
#' @param label_y character Label for the y axis (default \code{"Abundances"}).
#' @param label_fill character Label for color categories (default \code{"Features"}).
#' @param color Vector with custom colors for the different features. If empty, the default ggplot2 palette will be used (default \code{NULL}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotTaxonomy]{plotTaxonomy}} for plotting the most abundant taxa of a SQM object; \code{\link[plotBars]{plotHeatmap}} for plotting a heatmap with arbitrary data; \code{\link[mostAbundant]{mostAbundant}} for selecting the most abundant rows in a dataframe or matrix.
#' data(Hadza)
#' sk = Hadza$taxa$superkingdom$abund
#' plotBars(sk, label_y = 'Raw reads', label_fill = 'Superkingdom')
#' @export
plotBars = function(data, label_x = 'Samples', label_y = 'Abundances', label_fill = 'Features', color = NULL)
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
    data=t(data)
    data_melt = melt(as.matrix(data), value.name = 'abun')
    colnames(data_melt) = c('sample', 'item', 'abun')
    data_melt$abun = as.numeric(data_melt$abun)
    #PLOT DATA
    p = ggplot(data_melt, aes(x = sample, y = abun, fill = item))
    p = p + geom_col()
    # There are two types of bar charts: geom_bar() and geom_col().
    # geom_bar() makes the height of the bar proportional to the number of cases in each group (or if the weight aesthetic is supplied, the sum of the weights).
    # If you want the heights of the bars to represent values in the data, use geom_col() instead.
    # geom_bar() uses stat_count() by default: it counts the number of cases at each x position.
    # geom_col() uses stat_identity(): it leaves the data as is.
    # geom_col == geom_bar(stat = 'identity')

    if(is.null(label_x)) { label_x = '' }
    if(is.null(label_y)) { label_y = '' }
    if(is.na(label_x)  ) { label_x = '' }
    if(is.na(label_y)  ) { label_y = '' }

    if(nchar(label_x)==0) { theme_x = element_blank()
    }else{ theme_x = element_text(size = 12) }
    if(nchar(label_y)==0) { theme_y = element_blank()
    }else{ theme_y = element_text(size = 12) }

    p = p + theme_light()
    p = p + theme(axis.title.x = theme_x, axis.title.y = theme_y)
    p = p + labs(x = label_x, y = label_y, fill = label_fill)
    if (!is.null(color) & length(color) >= length(unique(as.character(data_melt$item))))
        {
        p = p + scale_fill_manual( values = setNames(color, unique(as.character(data_melt$item))) )
        }
    return(p)
    }



#' Heatmap of the most abundant functions in a SQM object
#'
#' This function selects the most abundant functions across all samples in a SQM object and represents their abundances in a heatmap. Alternatively, a custom set of functions can be represented.
#' @param SQM A SQM object.
#' @param fun_level character. Either \code{"KEGG"}, \code{"COG"} or \code{"PFAM"} (default \code{"KEGG"}).
#' @param count character. Either \code{"tpm"} for TPM normalized values, \code{"abund"} for raw abundances or \code{"copy_number"} for copy numbers (default \code{"tpm"}).
#' @param N integer Plot the \code{N} most abundant functions (default \code{25}).
#' @param fun character. Custom functions to plot. If provided, it will override \code{N} (default \code{NULL}).
#' @param gradient_col A vector of two colors representing the low and high ends of the color gradient (default \code{c("ghostwhite", "dodgerblue4")}).
#' @param ignore_unclassified logical. Don't include unclassified ORFs in the plot (default \code{TRUE}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotTaxonomy]{plotTaxonomy}} for plotting the most abundant taxa of a SQM object; \code{\link[plotBars]{plotBars}} and \code{\link[plotBars]{plotHeatmap}} for plotting barplots or heatmaps with arbitrary data.
#' @examples
#' data(Hadza)
#' plotFunctions(Hadza)
#' @export
plotFunctions = function(SQM, fun_level = 'KEGG', count = 'tpm', N = 25, fun = c(), gradient_col = c('ghostwhite', 'dodgerblue4'), ignore_unclassified = T)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    if ((fun_level != 'KEGG' & fun_level != 'COG' & fun_level != 'PFAM'))
        {
        stop('Select function category among "KEGG", "COG" or "PFAM"')
        }
    if(!count %in% c('abund', 'tpm', 'copy_number'))
        {
        stop('count must be "abund", "tpm" or "copy_number"')
        }
    if( count == 'copy_number' & is.null(SQM[['functions']][[fun_level]][[count]]) )
        {
        stop('There are no copy number tables in your project, possibly because RecA was not present in the metagenome')
        }
    if (is.null(fun) & N <= 0)
        {
        warning(sprintf('We can\'t plot N=%s functions. Continuing with default values', N))
        N = 25
        }
    # Work with samples in rows (like vegan). Tranposition converts a df into list again, need to cast it to df.
    data = as.data.frame(SQM[['functions']][[fun_level]][[count]])
    data = mostAbundant(data, N = N, items = fun)
    # remove unclassified functions (only possible when not custom items)
    # if N include unclassified function, add one more function
    if (ignore_unclassified & is.null(fun) & 'Unclassified' %in% rownames(data) & N != 0)
        {
        data = as.data.frame(SQM[['functions']][[fun_level]][[count]])
        data = mostAbundant(data, N = (N + 1), items = fun)
        data = data[rownames(data) != 'Unclassified', , drop = F]
        }
    # Add KEGG & COG names
    if (fun_level == 'KEGG' | fun_level == 'COG')
        {
        item_name = SQM$misc[[paste0(fun_level, '_names')]][rownames(data)]
        item_name[is.na(item_name)] = 'no_name'
        item_name = paste0(rownames(data) , '; ', item_name)
        item_name = sapply(item_name, FUN = function(x) paste(strwrap(x, width = 50), collapse = '\n'))
        rownames(data) = item_name
        }
    if('Unclassified' %in% rownames(data)) # Put unclassified at the bottom.
        {
        niceOrder = c(rownames(data)[rownames(data)!='Unclassified'], 'Unclassified')
        data = data[niceOrder,]
        }
    nice_label = c(abund='Raw abundance', tpm='TPM', copy_number='Copy number')[count]
    p = plotHeatmap(data, label_y = fun_level, label_fill = nice_label, gradient_col = gradient_col)
    return(p)
    }  



#' Barplot of the most abundant taxa in a SQM object
#'
#' This function selects the most abundant taxa across all samples in a SQM object and represents their abundances in a barplot. Alternatively, a custom set of taxa can be represented.
#' @param SQM A SQM object.
#' @param rank Taxonomic rank to plot (default \code{phylum}).
#' @param count character. Either \code{"percent"} for percentages, or \code{"abund"} for raw abundances (default \code{"percent"}).
#' @param N integer Plot the \code{N} most abundant taxa (default \code{15}).
#' @param tax character. Custom taxa to plot. If provided, it will override \code{N} (default \code{NULL}).
#' @param others logical. Collapse the abundances of least abundant taxa, and include the result in the plot (default \code{TRUE}).
#' @param color Vector with custom colors for the different features. If empty, we will use our own hand-picked pallete if N<=15, and the default ggplot2 palette otherwise (default \code{NULL}).
#' @param rescale logical. Re-scale results to percentages (default \code{FALSE}).
#' @param ignore_unclassified logical. Don't include unclassified contigs in the plot (default \code{FALSE}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotFunctions]{plotFunctions}} for plotting the most abundant functions of a SQM object; \code{\link[plotBars]{plotBars}} and \code{\link[plotBars]{plotHeatmap}} for plotting barplots or heatmaps with arbitrary data.
#' @examples
#' data(Hadza)
#' Hadza.amin = subsetFun(Hadza, 'Amino acid metabolism')
#' # Taxonomic distribution of amino acid metabolism ORFs at the family level.
#' plotTaxonomy(Hadza.amin, 'family')
#' @export
plotTaxonomy = function(SQM, rank = 'phylum', count = 'percent', N = 15, tax = NULL, others = T, color = NULL, rescale = F, ignore_unclassified = F)
    {
    if(!class(SQM)=='SQM') { stop('The first argument must be a SQM object') }
    if (!rank %in% c('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'))
        {
        stop('Select rank among "superkingdom", "phylum", "class", "order", "family", "genus" or "species". and count between \'percent\' or \'abund\'')
        }
    if (!count %in% c('abund', 'percent'))
        {
        stop('count must be either "abund" or "percent"')
        }
    if ('Other' %in% rownames(SQM[['taxa']][[rank]][[count]]))
        {
        stop('One of your taxa is called "Other", please change its name')
        }
    if (is.null(tax) & N <= 0)
        {
        warning(sprintf('We can\'t plot N = %s? taxa. Continuing with default values', N))
        N = 15
        }
    # Work with samples in rows (like vegan). Tranposition converts a df into list again, need to cast it to df.
    data = as.data.frame(SQM[['taxa']][[rank]][[count]])
    data = mostAbundant(data, N = N, items = tax, others = others, rescale = rescale)
    # remove unclassified taxa (only possible when not custom items)
    # if N include unclassified taxa, add one more taxa
    if (ignore_unclassified & is.null(tax) & 'Unclassified' %in% rownames(data) & N != 0)
        { # We overwrite data from scratch!
        data = as.data.frame(SQM[['taxa']][[rank]][[count]]) # Pick one more taxa
        data = mostAbundant(data, N = N + 1, items = tax, others = others, rescale = F) # Create the data table again
        data = data[rownames(data) != 'Unclassified', , drop = F] # Remove 'Unclassified'
        data = mostAbundant(data, items = rownames(data), others = F, rescale = rescale) # Renormalize/Others
    }
    
    #print(data)
    defaultColors = c(
        '#97B065', '#B93838', '#83AAF0', '#E2AA48', '#A9A3D2',
        '#797900', '#F4BCE6', '#84F4F6', '#9AD63B', '#A285F5',
        '#D2B48C', '#4EA24E', '#465569', '#9C669C', '#6495ED')
    # Colors to plot. Checks.
    if (!is.null(color))
        {
        # Try to use user colors
        # Check if the user is trying to trick us
        if(!(is.null(tax)) & length(color) == length(tax) )
            {
            # User passes colors for tax
            color = color
            } else if(is.null(tax) & length(color) == N)
            {
            # User passes colors for the N most abundant taxa, use them
            color = color
            } else
            {
            # User passes less/more colors than taxa, use default colors.
            warning('You passed less/more colors than taxa. Using default colors')
            color = defaultColors
            }
    } else
        {
        # Use default colors from the beginning. User does not care about colors
        if(N <= length(defaultColors)) { color = defaultColors[1:nrow(data[!rownames(data) %in% c('Other', 'Unclassified'),,drop=F])]
        }else{ color = NULL }
        }
 
    # Add others color
    if (others & !is.null(color))
        {
        color = c('#F5DEB3', color)
        }
    # Add unclassified color and put Unclassified at the bottom.
    if('Unclassified' %in% rownames(data))
        {
        if(!is.null(color)) { color = c(color, 'azure3') }
        niceOrder = c(rownames(data)[rownames(data)!='Unclassified'], 'Unclassified')
        data = data[niceOrder,,drop=F]
        }
    nice_label = c(abund='Raw abundance', percent='Percentage')[count]
    nice_rank  = paste0(toupper(substr(rank,1,1)), substr(rank,2,nchar(rank)))
    p = plotBars(data, label_y = nice_label, color = color, label_fill = nice_rank)
    return(p)
    }

