require(ggplot2)
require(reshape2)

#' Plot a heatmap using ggplot2
#'
#' Plot a ggplot2 heatmap from a matrix or data frame. The data should be in tabular format (e.g. features in rows and samples in columns).
#' @param data numeric matrix or data frame.
#' @param label_x character Label for the x axis (default \code{"Samples"}).
#' @param label_y character Label for the y axis (default \code{"Features"}).
#' @param label_fill character Label for color scale (default \code{"Abundance"}).
#' @param gradient_col A vector of two colors representing the low and high ends of the color gradient (default \code{c("ghostwhite", "dodgerblue4")}).
#' @param base_size numeric. Base font size (default \code{11}).
#' @param metadata_groups list. Split the plot into groups defined by the user: list('G1' = c('sample1', sample2'), 'G2' = c('sample3', 'sample4')) default \code{NULL}).
#' @return A ggplot2 plot object.
#' @seealso \code{\link[plotFunctions]{plotFunctions}} for plotting the top functional categories of a SQM object; \code{\link[plotBars]{plotBars}} for plotting a barplot with arbitrary data; \code{\link[mostAbundant]{mostAbundant}} for selecting the most abundant rows in a dataframe or matrix.
#' @examples
#' data(Hadza)
#' topPFAM = mostAbundant(Hadza$functions$PFAM$tpm)
#' topPFAM = topPFAM[rownames(topPFAM) != "Unclassified",] # Take out the Unclassified ORFs.
#' plotHeatmap(topPFAM, label_x = "Samples", label_y = "PFAMs", label_fill = "TPM")
#' @export
plotHeatmap = function(data, label_x = 'Samples', label_y = 'Features', label_fill = 'Abundance', gradient_col = c('ghostwhite', 'dodgerblue4'), base_size = 11, metadata_groups = NULL)
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
    if(length(gradient_col) < 2)
        {
        stop('gradient_col must be a vector with two colors representing the low and high ends of the color gradient')
        }
    data=t(data)
    # data = data[, order(colSums(data), decreasing = T), drop = F] # Order functions according to their abundances
    data_melt = reshape2::melt(as.matrix(data), value.name = 'abun')
    colnames(data_melt) = c('sample', 'item', 'abun')
    data_melt$sample = as.factor(data_melt$sample)
    data_melt$item = as.factor(data_melt$item)
    data_melt$abun = as.numeric(data_melt$abun)
    if (!is.null(metadata_groups))
    	{
            if(sum(sapply(metadata_groups, length)) != length(unique(data_melt$sample)))
                {stop('metadata_groups must contain the same samples that data')}
        data_melt$group = apply(data_melt, 1, function(s) unlist(lapply(names(metadata_groups), function(x) if( s['sample'] %in% metadata_groups[[x]]){x})))
    	}
    #PLOT DATA
    if(is.null(label_x)) { label_x = '' }
    if(is.null(label_y)) { label_y = '' }
    if(is.na(label_x)  ) { label_x = '' }
    if(is.na(label_y)  ) { label_y = '' }

    if(nchar(label_x)==0) { theme_x = ggplot2::element_blank()
    }else{ theme_x = ggplot2::element_text() }
    if(nchar(label_y)==0) { theme_y = ggplot2::element_blank()
    }else{ theme_y = ggplot2::element_text() }
        
    p = ggplot2::ggplot(data_melt, ggplot2::aes(x = sample, y = item, fill = abun))
    p = p + ggplot2::geom_tile()
    p = p + ggplot2::scale_fill_gradient(low = gradient_col[1], high = gradient_col[2])
    p = p + ggplot2::theme_light(base_size = base_size)
    p = p + ggplot2::theme(axis.title.x = theme_x, axis.title.y = theme_y)
    p = p + ggplot2::labs(x = label_x, y = label_y, fill = label_fill)
    if (!is.null(metadata_groups))
    	{
        p = p + ggplot2::facet_grid(.~group, scale = 'free')
    	}	

    if (length(unique(data_melt$sample)) > 10)
    	{
        p = p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
    	}
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
#' @param base_size numeric. Base font size (default \code{11}).
#' @param max_scale_value numeric. Maximum value to include in the y axis. By default it is handled automatically by ggplot2 (default \code{NULL}).
#' @param metadata_groups list. Split the plot into groups defined by the user: list('G1' = c('sample1', sample2'), 'G2' = c('sample3', 'sample4')) default \code{NULL}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotTaxonomy]{plotTaxonomy}} for plotting the most abundant taxa of a SQM object; \code{\link[plotBars]{plotHeatmap}} for plotting a heatmap with arbitrary data; \code{\link[mostAbundant]{mostAbundant}} for selecting the most abundant rows in a dataframe or matrix.
#' @examples
#' data(Hadza)
#' sk = Hadza$taxa$superkingdom$abund
#' plotBars(sk, label_y = "Raw reads", label_fill = "Superkingdom")
#' @export
plotBars = function(data, label_x = 'Samples', label_y = 'Abundances', label_fill = 'Features', color = NULL, base_size = 11, max_scale_value = NULL, metadata_groups = NULL)
    {
    if (!is.data.frame(data) & !is.matrix(data)) { stop('The first argument must be a matrix or a data frame') }
    if(!is.null(max_scale_value) & !is.numeric(max_scale_value)) { stop('max_scale_value must be numeric') }
    data=t(data)
    data_melt = reshape2::melt(as.matrix(data), value.name = 'abun')
    colnames(data_melt) = c('sample', 'item', 'abun')
    data_melt$sample = as.factor(data_melt$sample)
    data_melt$item = as.factor(data_melt$item)
    data_melt$abun = as.numeric(data_melt$abun)
    
     if (!is.null(metadata_groups)) 
    	{
	    if(sum(sapply(metadata_groups, length)) != length(unique(data_melt$sample)))
	        {stop('metadata_groups must contain the same samples that data')}	    
        data_melt$group = apply(data_melt, 1, function(s) unlist(lapply(names(metadata_groups), function(x) if( s['sample'] %in% metadata_groups[[x]]){x})))
    	}
	
    #PLOT DATA
    p = ggplot2::ggplot(data_melt, ggplot2::aes(x = sample, y = abun, fill = item))
    p = p + ggplot2::geom_col()
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

    if(nchar(label_x)==0) { theme_x = ggplot2::element_blank()
    }else{ theme_x = ggplot2::element_text() }
    if(nchar(label_y)==0) { theme_y = ggplot2::element_blank()
    }else{ theme_y = ggplot2::element_text() }

    p = p + ggplot2::theme_light(base_size = base_size)
    p = p + ggplot2::theme(axis.title.x = theme_x, axis.title.y = theme_y)
    p = p + ggplot2::labs(x = label_x, y = label_y, fill = label_fill)
    if (!is.null(color) & length(color) >= length(unique(as.character(data_melt$item))))
        {
        p = p + ggplot2::scale_fill_manual( values = setNames(color, unique(as.character(data_melt$item))) )
        }
    if (!is.null(max_scale_value))
        {
        p = p + ggplot2::ylim(0, max_scale_value)
        }
    if (!is.null(metadata_groups))
     	{
        p = p + ggplot2::facet_grid(.~group, scale = 'free')
    	}

    if (length(unique(data_melt$sample)) > 10)
    	{
        p = p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
    	}	
    return(p)
    }



#' Heatmap of the most abundant functions in a SQM object
#'
#' This function selects the most abundant functions across all samples in a SQM object and represents their abundances in a heatmap. Alternatively, a custom set of functions can be represented.
#' @param SQM A SQM or SQMlite object.
#' @param fun_level character. Either \code{"KEGG"}, \code{"COG"}, \code{"PFAM"} or any other custom database used for annotation (default \code{"KEGG"}).
#' @param count character. Either \code{"abund"} for raw abundances, \code{"percent"} for percentages, \code{"bases"} for raw base counts, \code{"tpm"} for TPM normalized values or \code{"copy_number"} for copy numbers (default \code{"tpm"}). Note that a given count type might not available in this object (e.g. TPM or copy number in SQMlite objects originating from a SQM reads project).
#' @param N integer Plot the \code{N} most abundant functions (default \code{25}).
#' @param fun character. Custom functions to plot. If provided, it will override \code{N} (default \code{NULL}).
#' @param samples character. Character vector with the names of the samples to include in the plot. Can also be used to plot the samples in a custom order. If not provided, all samples will be plotted (default \code{NULL}).
#' @param ignore_unmapped logical. Don't include unmapped ORFs in the plot (default \code{TRUE}).
#' @param ignore_unclassified logical. Don't include unclassified ORFs in the plot (default \code{TRUE}).
#' @param gradient_col A vector of two colors representing the low and high ends of the color gradient (default \code{c("ghostwhite", "dodgerblue4")}).
#' @param base_size numeric. Base font size (default \code{11}).
#' @param metadata_groups list. Split the plot into groups defined by the user: list('G1' = c('sample1', sample2'), 'G2' = c('sample3', 'sample4')) default \code{NULL}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotTaxonomy]{plotTaxonomy}} for plotting the most abundant taxa of a SQM object; \code{\link[plotBars]{plotBars}} and \code{\link[plotBars]{plotHeatmap}} for plotting barplots or heatmaps with arbitrary data.
#' @examples
#' data(Hadza)
#' plotFunctions(Hadza)
#' @export
plotFunctions = function(SQM, fun_level = 'KEGG', count = 'tpm', N = 25, fun = NULL, samples = NULL, ignore_unmapped = T, ignore_unclassified = T, gradient_col = c('ghostwhite', 'dodgerblue4'), base_size = 11, metadata_groups = NULL)
    {
    if(!class(SQM) %in% c('SQM', 'SQMlite')) { stop('The first argument must be a SQM or a SQMlite object') }
    if (!fun_level %in% names(SQM$functions))
        {
        intro = 'Select function category among'
        opts  = paste(names(SQM$functions), collapse='" "')
        msg   = sprintf('%s "%s"', intro, opts)
        stop(msg)
        }
    if(!count %in% c('abund', 'percent', 'bases', 'tpm', 'copy_number'))
        {
        stop('count must be "abund", "percent", "bases", "tpm" or "copy_number"')
        }
    if( count == 'bases' & is.null(SQM[['functions']][[fun_level]][[count]]) )
        {
        stop('There are no tables containing aggregated bases per function in your project, possibly because it was it was generated by sqm_reads.pl. Try with "count=\'abund\'".')
        }
    if( count == 'tpm' & is.null(SQM[['functions']][[fun_level]][[count]]) )
        {
        stop('There are no tables containing TPM per function in your project, possibly because it was it was generated by sqm_reads.pl. Try with "count=\'percent\'".')
        }
    if( count == 'copy_number' & is.null(SQM[['functions']][[fun_level]][[count]]) )
        {
        stop('There are no copy number tables in your project, possibly because RecA was not present in the metagenome or it was generated by sqm_reads.pl. Try with "count=\'percent\'".')
        }
    if (is.null(fun) & N <= 0)
        {
        warning(sprintf('We can\'t plot N=%s functions. Continuing with default values', N))
        N = 25
        }

    check.samples(SQM, samples)

    # Work with samples in rows (like vegan). Tranposition converts a df into list again, need to cast it to df.
    if(count == 'percent')
        {
        percents = 100 * t(t(SQM[['functions']][[fun_level]][['abund']]) / colSums(SQM[['functions']][[fun_level]][['abund']]))
        data = as.data.frame(percents)
        }else { data = as.data.frame(SQM[['functions']][[fun_level]][[count]]) }
    data = mostAbundant(data, N = N, items = fun)

    # Verify whether there are Unclassified or Unmapped
    if(ignore_unmapped & !('Unmapped' %in% rownames(data))) { ignore_unmapped = F }
    if(ignore_unclassified & !('Unclassified' %in% rownames(data))) { ignore_unclassified = F }
    
    # remove unmapped functions (only possible when not custom items)
    # if N include unmapped function, add one more function
    if (ignore_unmapped & !(ignore_unclassified) & is.null(fun) & N != 0)
        {
        if(count == 'percent') { data = as.data.frame(percents)
        } else { data = as.data.frame(SQM[['functions']][[fun_level]][[count]]) }
        data = mostAbundant(data, N = (N + 1), items = fun)
        data = data[rownames(data) != 'Unmapped', , drop = F] # Remove 'Unmapped'
        }

    # remove unclassified functions (only possible when not custom items)
    # if N include unclassified function, add one more function
    if (ignore_unclassified & !(ignore_unmapped) & is.null(fun) & N != 0)
        {
        if(count == 'percent') { data = as.data.frame(percents)
        }else { data = as.data.frame(SQM[['functions']][[fun_level]][[count]]) }
        data = mostAbundant(data, N = (N + 1), items = fun)
        data = data[rownames(data) != 'Unclassified', , drop = F] # Remove 'Unclassified'
        }

    # remove unmapped and unclassified functions (only possible when not custom items)
    # if N include unmapped and unclassified functions, add two more functions
    if (ignore_unclassified & ignore_unmapped & is.null(fun) & N != 0)
        {
         if(count == 'percent') { data = as.data.frame(percents)
         }else { data = as.data.frame(SQM[['functions']][[fun_level]][[count]]) }
         data = mostAbundant(data, N = (N + 2), items = fun)
         data = data[rownames(data) != 'Unmapped', , drop = F]
         data = data[rownames(data) != 'Unclassified', , drop = F]
        }

    # Add KEGG & COG names
    if (fun_level %in% c('COG', 'KEGG'))
        {
        item_name = SQM$misc[[paste0(fun_level, '_names')]][rownames(data)]
        item_name[is.na(item_name)] = 'no_name'
        item_name = paste0(rownames(data) , '; ', item_name)
        item_name = sapply(item_name, FUN = function(x) paste(strwrap(x, width = 50), collapse = '\n'))
        rownames(data) = item_name
        }
    # Put unclassified and unmapped at the bottom
    if('Unclassified' %in% rownames(data)) # Put unclassified at the bottom.
        {
        niceOrder = c(rownames(data)[rownames(data)!='Unclassified'], 'Unclassified')
        data = data[niceOrder,]
        }
    if('Unmapped' %in% rownames(data)) # Put unmapped at the bottom.
        {
        niceOrder = c(rownames(data)[rownames(data)!='Unmapped'], 'Unmapped')
        data = data[niceOrder,]
        }
    nice_label = c(abund='Raw abundance', percent = 'Percentage', bases='Bases', tpm='TPM', copy_number='Copy number')[count]

    # If requested, plot only the selected samples
    if(!is.null(samples)) { data = data[,samples,drop=F] }
    if (!is.null(metadata_groups))
        {
            if(sum(sapply(metadata_groups, length)) != ncol(data))
               {stop('metadata_groups must contain the same samples that your SQM object. If you want to select some samples, you can use the flag samples and a metadata_groups list that includes all the selected samples')}
        }
    # Plot
    p = plotHeatmap(data, label_y = fun_level, label_fill = nice_label, gradient_col = gradient_col, base_size = base_size, metadata_groups = metadata_groups)
    return(p)
    }  



#' Barplot of the most abundant taxa in a SQM object
#'
#' This function selects the most abundant taxa across all samples in a SQM object and represents their abundances in a barplot. Alternatively, a custom set of taxa can be represented.
#' @param SQM A SQM or a SQMlite object.
#' @param rank Taxonomic rank to plot (default \code{phylum}).
#' @param count character. Either \code{"percent"} for percentages, or \code{"abund"} for raw abundances (default \code{"percent"}).
#' @param N integer Plot the \code{N} most abundant taxa (default \code{15}).
#' @param tax character. Custom taxa to plot. If provided, it will override \code{N} (default \code{NULL}).
#' @param others logical. Collapse the abundances of least abundant taxa, and include the result in the plot (default \code{TRUE}).
#' @param nocds character. Either \code{"treat_separately"} to treat reads annotated as No CDS separately, \code{"treat_as_unclassified"} to treat them as Unclassified or \code{"ignore"} to ignore them in the plot (default \code{"treat_separately"}).
#' @param ignore_unmapped logical. Don't include unmapped reads in the plot (default \code{FALSE}).
#' @param ignore_unclassified logical. Don't include unclassified reads in the plot (default \code{FALSE}).
#' @param samples character. Character vector with the names of the samples to include in the plot. Can also be used to plot the samples in a custom order. If not provided, all samples will be plotted (default \code{NULL}).
#' @param no_partial_classifications logical. Treat reads not fully classified at the requested level (e.g. "Unclassified bacteroidetes" at the class level or below) as fully unclassified. This takes effect before \code{ignore_unclassified}, so if both are \code{TRUE} the plot will only contain fully classified contigs (default \code{FALSE}).
#' @param rescale logical. Re-scale results to percentages (default \code{FALSE}).
#' @param color Vector with custom colors for the different features. If empty, we will use our own hand-picked pallete if N<=15, and the default ggplot2 palette otherwise (default \code{NULL}).
#' @param base_size numeric. Base font size (default \code{11}).
#' @param max_scale_value numeric. Maximum value to include in the y axis. By default it is handled automatically by ggplot2 (default \code{NULL}).
#' @return a ggplot2 plot object.
#' @seealso \code{\link[plotFunctions]{plotFunctions}} for plotting the most abundant functions of a SQM object; \code{\link[plotBars]{plotBars}} and \code{\link[plotBars]{plotHeatmap}} for plotting barplots or heatmaps with arbitrary data.
#' @examples
#' data(Hadza)
#' Hadza.amin = subsetFun(Hadza, "Amino acid metabolism")
#' # Taxonomic distribution of amino acid metabolism ORFs at the family level.
#' plotTaxonomy(Hadza.amin, "family")
#' @export
plotTaxonomy = function(SQM, rank = 'phylum', count = 'percent', N = 15, tax = NULL, others = T, samples = NULL, nocds = 'treat_separately', ignore_unmapped = F, ignore_unclassified = F, no_partial_classifications = F, rescale = F, color = NULL, base_size = 11, max_scale_value = NULL, metadata_groups = NULL)
    {
    if(!class(SQM) %in% c('SQM', 'SQMlite')) { stop('The first argument must be a SQM or a SQMlite object') }
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
    if (!nocds %in% c('treat_separately', 'treat_as_unclassified', 'ignore'))
        {
        stop('Select nocds among "treat_separately", "treat_as_unclassified", "ignore"')
        }
    if(!is.null(max_scale_value) & !is.numeric(max_scale_value)) { stop('max_scale_value must be numeric') }

    #if(!is.null(ignore_taxa)) 
    #    {
    #    remove_taxa = rep(T, length(ignore_taxa))
    #    names(remove_taxa) = ignore_taxa
    #    }
    check.samples(SQM, samples)

    data0 = SQM[['taxa']][[rank]][[count]]
    # First collapse partial classifications if required
    if(no_partial_classifications)
        {
        unclassified = grepl('[Uu]nclassified', rownames(data0))
        unclassified_counts = colSums(data0[unclassified,,drop=F])
	data0 = rbind(data0[!unclassified,,drop=F], 'Unclassified' = unclassified_counts)
        }

    # Work with samples in rows (like vegan). Tranposition converts a df into list again, need to cast it to df
    data = as.data.frame(data0)
    data = mostAbundant(data, N = N, items = tax, others = others, rescale = rescale)

    # Verify whether there are Unclassified or Unmapped or No CDS



    if(ignore_unmapped & !('Unmapped' %in% rownames(data))) {ignore_unmapped = F}
    if((nocds == 'ignore'| nocds == 'treat_as_unclassified') & !('No CDS' %in% rownames(data))) { nocds = 'treat_separately' }
    # Convert No CDS options in logical values
    if (nocds == 'treat_separately') 
        { 
        ignore_nocds = F 
        } else if (nocds == 'ignore') 
        { 
        ignore_nocds = T 
        } else if (nocds == 'treat_as_unclassified') 
        { 
        ignore_nocds = T # at the end we ignore this class 'cause they're included in ignore_unclassified 
        }
    if(ignore_unclassified & !('Unclassified' %in% rownames(data)) & nocds != 'treat_as_unclassified') {ignore_unclassified = F}

    ignore_cases = c('Unmapped' = ignore_unmapped, 'Unclassified' = ignore_unclassified, 'No CDS' = ignore_nocds)
    #if (!is.null(ignore_taxa)) {ignore_cases = c(ignore_cases, remove_taxa)}
    
    # any(ignore_cases) => return T if at least one of the cases is T and we should ignore it
    # remove unmapped/noCDS/Unclassified taxa (only possible when not custom items)
    # Add as many taxa as needed to recover N taxa excluding any of the special cases
    if ( any(ignore_cases) & is.null(tax)  & N != 0 )
        { # We overwrite data from scratch!
        data = as.data.frame(data0) # Pick more taxa to complete N
        # how many taxa should we replace?
        rr = sum(ignore_cases) # count T cases
        data = mostAbundant(data, N = N + rr, items = tax, others = others, rescale = F) # Create the data table again
        # If nocds == 'treat_as_unclassified' => Unclassified = Unclassified + nocds
        if (nocds == 'treat_as_unclassified') 
            {
            if (('Unclassified' %in% rownames(data)))
                {
                data['Unclassified', ] = data['Unclassified', ] + data['No CDS', ]
                } else {data['Unclassified', ] = data['No CDS', ]}
            }
	# Remove ignore_cases = T
        for (i in 1:length(ignore_cases))
            {
            if (ignore_cases[i])
                { data = data[rownames(data) != names(ignore_cases)[i], , drop = F] } # Remove  ignore cases
            }
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
                }else if(is.null(tax) & length(color) == N)
                {
                # User passes colors for the N most abundant taxa, use them
                color = color
                }else
                {
                # User passes less/more colors than taxa, use default colors.
                warning('You passed less/more colors than taxa. Using default colors')
                color = defaultColors
                }
        }else
        {
            # Use default colors from the beginning. User does not care about colors
            if(N <= length(defaultColors)) { color = defaultColors[1:nrow(data[!rownames(data) %in% c('Other', 'Unclassified', 'Unmapped', 'No CDS'), , drop=F])]
            }else{ color = NULL }
        }

    # Add others color
    if (others & !is.null(color))
        {
        color = c('#F5DEB3', color)
        }

    # Add NoCDS color and put No CDS at the bottom
    if('No CDS' %in% rownames(data))
        {
        if(!is.null(color)) { color = c(color, 'azure2') }
        niceOrder = c(rownames(data)[rownames(data)!='No CDS'], 'No CDS')
        data = data[niceOrder,,drop=F]
        }
    
    # Add unclassified color and put Unclassified at the bottom
    if('Unclassified' %in% rownames(data))
        {
        if(!is.null(color)) { color = c(color, 'azure3') }
        niceOrder = c(rownames(data)[rownames(data)!='Unclassified'], 'Unclassified')
        data = data[niceOrder,,drop=F]
        }

    # Add unmapped color and put Unmapped at the bottom
    if('Unmapped' %in% rownames(data))
        {
        if(!is.null(color)) { color = c(color, 'azure4') }
        niceOrder = c(rownames(data)[rownames(data)!='Unmapped'], 'Unmapped')
        data = data[niceOrder,,drop=F]
        }




    # If requested, plot only the selected samples
    if(!is.null(samples)) { data = data[,samples,drop=F] }
    if (!is.null(metadata_groups))
        {
            if(sum(sapply(metadata_groups, length)) != ncol(data))
               {stop('metadata_groups must contain the same samples that your SQM object. If you want to select some samples, you can use the flag samples and a metadata_groups list that includes all the selected samples')}
        }
    # Plot
    nice_label = c(abund='Raw abundance', percent='Percentage')[count]
    nice_rank  = paste0(toupper(substr(rank,1,1)), substr(rank,2,nchar(rank)))
    p = plotBars(data, label_y = nice_label, color = color, label_fill = nice_rank, base_size = base_size, max_scale_value = max_scale_value, metadata_groups = metadata_groups)
    return(p)
    }
