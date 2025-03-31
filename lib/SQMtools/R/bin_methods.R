utils::globalVariables(c('CheckMProkaryote'))

#' Find redundant contigs within a bin
#'
#' Find contigs with overlapping marker genes in a given bin, and suggest contigs to be removed in order to reduce contamination without affecting completeness. Note that this can give a quick idea of the contigs that are sources of contamination within a bin, but is not a replacement for proper bin refininement with other tools such as anvi\'o.
#'
#' @param SQM A SQM object.
#' @param bin character. Name of the bin to be created.
#' @param minimum_overlap_for_removal numeric. Fraction of marker genes in the contigs present in another contig needed to suggest it for removal. If set to \code{1} (default), contigs will only suggested for removal if their markers fully overlap with those in another contig (and thus completeness will not change after removing them). Smaller values will result in more contigs being suggested for removal, which will further reduce contamination at the expense of some completeness.
#' @return A character vector with the contigs deemed to be redundant. A heatmap showing how marker genes overlap over different contigs will also be produced. 
#' @seealso \code{\link{create_bin}}, \code{\link{remove_contigs_from_bin}}
#' @examples
#' data(Hadza)
#' bin_name = "Hadza2merged.concoct.28.fa.contigs"
#' # Get redundant contigs that could be removed from our bin
#' candidates_for_removal = find_redundant_contigs(Hadza, bin_name)
#' # We can now remove them from the bin
#' Hadza.new.1 = remove_contigs_from_bin(Hadza, bin_name, candidates_for_removal)
#' # Or we can create a new bin out of them
#' #  which will also remove them from the original bin
#' Hadza.new.2 = create_bin(Hadza, "new_bin_name", candidates_for_removal)
#' @export
#' @importFrom utils combn
find_redundant_contigs = function(SQM, bin, minimum_overlap_for_removal = 1) {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(!('tax16S' %in% names(SQM$orfs)) & 'markers' %in% names(SQM$orfs)) {
        stop('16S or marker gene info are missing from this dataset, can not compute marker gene overlaps')
        }
    #bin = 'Hadza2merged.concoct.28.fa.contigs'
    warning('While this can give a quick idea of which contigs are sources of contamination, it is not a replacement for proper bin refinenent e.g. with anvi\'o')
    contigs = rownames(SQM$contigs$bins)[SQM$contigs$bins==bin]
    markers = matrix(0, nrow=length(contigs), ncol=length(unique(unlist(CheckMProkaryote))), dimnames=list(contigs, unique(unlist(CheckMProkaryote))))
    orfs = rownames(SQM$orfs$table)[SQM$orfs$table$Contig %in% contigs]
    for(orf in orfs) {
        ms = intersect(SQM$orfs$markers[[orf]], colnames(markers))
        if(!length(ms)) { next }
        ctg = SQM$orfs$table[orf, 'Contig ID']
        for(m in ms) {
            markers[ctg,m] = markers[ctg,m] + 1
            }
        }

    # Sort contigs based on how many markers they have
    #  and markers based on in how many contigs they appear
    markers = markers[names(sort(rowSums(markers), decreasing=T)),]
    markers = markers[,names(sort(colSums(markers), decreasing=T))]

    # Get only contigs with markers
    markers = markers[rowSums(markers)>0,]
    # Get overlaps
    overlaps = matrix(1, nrow=nrow(markers), ncol=nrow(markers), dimnames = list(rownames(markers), rownames(markers)))
    combis = combn(rownames(markers),2)
    for(i in 1:ncol(combis)) {
        c1 = combis[1,i]
        c2 = combis[2,i]
        m1 = markers[c1,] == 1
        m2 = markers[c2,] == 1
        ol = sum(m1&m2)
        overlaps[c1,c2] = ol/sum(m1)
        overlaps[c2,c1] = ol/sum(m2)
        }

    # Plot overlaps
    plotHeatmap(overlaps, label_x='Target contig', label_y='Source contig', label_fill='Fraction of markers from source present in target')

    # Get candidates for removal
    removal_candidates = c()
    for(ctg in rownames(overlaps)) {
        if(ctg %in% removal_candidates) { next }
        potential_targets = setdiff(colnames(overlaps), c(ctg, removal_candidates))
        removal_candidates = c(removal_candidates, potential_targets[overlaps[potential_targets,ctg]>=1])
        }

    return(removal_candidates)
    }


#' Create a bin from a vector of contigs
#'
#' @param SQM A SQM object.
#' @param bin character. Name of the bin to be created.
#' @param contigs character. Vector with the names of the contigs that will be included in the new bin.
#' @param delete_overlapping_bins logical. If \code{TRUE}, bins that originally contained any of the provided contigs will be removed from the object. Default \code{FALSE}.
#' @return SQM object with the new binning information, including recalculated bin statistics if possible.
#' @seealso \code{\link{find_redundant_contigs}}, \code{\link{remove_contigs_from_bin}}
#' @export
create_bin = function(SQM, bin, contigs, delete_overlapping_bins = FALSE) {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(!('tax16S' %in% names(SQM$orfs)) & 'markers' %in% names(SQM$orfs)) {
        stop('16S or marker gene info are missing from this dataset, can not recalculate bin stats')
        }
    newSQM = SQM
    oldbins = unique(newSQM$contigs$bins[contigs,1])
    newSQM$contigs$bins[contigs,1] = bin
    if(delete_overlapping_bins) {
        newSQM$contigs$bins[newSQM$contigs$bins[,1] %in% oldbins,1] = 'No_bin'
        bins_in_new_table = bin
    } else { bins_in_new_table = c(oldbins, bin) }

    new_stats = get.bin.stats(newSQM, bins_in_new_table)
    newSQM$bins$table[bins_in_new_table,] = new_stats[['table']][bins_in_new_table,]
    newSQM$bins$table[bins_in_new_table,'Method'] = 'Custom'
    newSQM$bins$tax = as.data.frame(newSQM$bins$tax)
    newSQM$bins$tax  [bins_in_new_table,] = new_stats[['tax'  ]][bins_in_new_table,]
    newSQM$bins$tax = as.matrix(newSQM$bins$tax)
    if(delete_overlapping_bins) {
        newSQM$bins$table = newSQM$bins$table[!rownames(newSQM$bins$table) %in% oldbins,]
        newSQM$bins$tax   = newSQM$bins$tax  [!rownames(newSQM$bins$tax)   %in% oldbins,]
        }
    bin_abunds                          = get.bin.abunds(newSQM)
    newSQM$bins$abund                   = bin_abunds[['abund']]
    newSQM$bins$percent                 = bin_abunds[['percent']]
    newSQM$bins$bases                   = bin_abunds[['bases']]
    newSQM$bins$length                  = bin_abunds[['length']]
    newSQM$bins$cov                     = bin_abunds[['cov']]
    newSQM$bins$cpm                     = bin_abunds[['cpm']]
    newSQM$bins$table                   = newSQM$bins$table[names(newSQM$bins$length),]
    warning('Bin stats were only recalculated for the modified bins')
    return(newSQM)
    }


#' Remove contigs from a given bin
#'
#' @param SQM A SQM object.
#' @param bin character. Name of the bin from which the contigs will be removed.
#' @param contigs character. Vector with the names of the contigs that will be removed from the new bin.
#' @return SQM object with the new binning information, including recalculated bin statistics if possible.
#' @seealso \code{\link{find_redundant_contigs}}, \code{\link{create_bin}}
#' @export
remove_contigs_from_bin = function(SQM, bin, contigs) {
    if(!inherits(SQM, 'SQM')) { stop('The first argument must be a SQM object') }
    if(!('tax16S' %in% names(SQM$orfs)) & 'markers' %in% names(SQM$orfs)) {
        stop('16S or marker gene info are missing from this dataset, can not recalculate bin stats')
        }
    all_contigs = rownames(SQM$contigs$bins)[SQM$contigs$bins[,1]==bin]
    wrong_contigs = setdiff(contigs, all_contigs)
    if(length(wrong_contigs)>0) {
        stop(sprintf('The following contig does not belong to bin %s: %s\n', bin, contigs))
        }
    newSQM = SQM
    newSQM$contigs$bins[contigs,1] = 'No_bin'
    new_stats = get.bin.stats(newSQM, bin)
    newSQM$bins$table[bin,] = new_stats[['table']][bin,]
    newSQM$bins$table[bin,'Method'] = 'Custom'
    newSQM$bins$tax  [bin,] = new_stats[['tax'  ]][bin,]
    bin_abunds                          = get.bin.abunds(newSQM)
    newSQM$bins$abund                   = bin_abunds[['abund']]
    newSQM$bins$percent                 = bin_abunds[['percent']]
    newSQM$bins$bases                   = bin_abunds[['bases']]
    newSQM$bins$length                  = bin_abunds[['length']]
    newSQM$bins$cov                     = bin_abunds[['cov']]
    newSQM$bins$cpm                     = bin_abunds[['cpm']]
    newSQM$bins$table                   = newSQM$bins$table[names(newSQM$bins$length),]
    warning('Bin stats were only recalculated for the modified bins')
    return(newSQM)
    }


get.bin.stats = function(SQM, bins = NULL)
    {
    if(is.null(bins)) { bins = unique(SQM$contigs$bins[,1]) }
    bins = sort(bins[bins!='No_bin'])
    GC            = c()
    nContigs      = c()
    disparity     = c()
    completeness  = c()
    contamination = c()
    tax16S        = c()
    taxonomy      = matrix(NA, nrow = length(bins), ncol = ncol(SQM$contigs$tax),
                   dimnames = list(bins, colnames(SQM$contigs$tax)))
    if(!'markers' %in% names(SQM$orfs))
        {
        warning('There are no universal marker annotations in your project. Skipping completeness/contamination calculations')
    } else
        {
        warning('Bin completeness / contamination will be recalculated using the CheckM1 algorithm on root markers')
        }
    for(b in bins) {
        contigs  = rownames(SQM$contigs$bins)[SQM$contigs$bins[,1]==b]
        nContigs = c(nContigs, length(contigs))
        bl       = sum(SQM$contigs$table[contigs,'Length'])
        gc       = sum(SQM$contigs$table[contigs,'GC perc'] * SQM$contigs$table[contigs,'Length']) / bl
        GC       = c(GC, round(gc,2))
        orfs     = rownames(SQM$orfs$table)[SQM$orfs$table$`Contig ID` %in% contigs]
        orfs16S  = orfs[orfs %in% names(SQM$orfs$tax16S)]
        tax16S   = c(tax16S, paste(SQM$orfs$tax16S[orfs16S], collapse='|'))

        if(!'markers' %in% names(SQM$orfs)) {
            comp = NA
            cont = NA
        } else {

            comp = 0
            cont = 0

            marker_count = table(unlist(SQM$orfs$markers[orfs]))

            for(collocated_marker_set in CheckMProkaryote)
                {
                present_markers = marker_count[intersect(collocated_marker_set, names(marker_count))]
                comp = comp + sum(present_markers>0)/length(collocated_marker_set)
                extra_markers = present_markers - 1
                extra_markers[extra_markers<0] = 0
                cont = cont + sum(extra_markers)/length(collocated_marker_set)
                }

            comp  = 100 * comp / length(CheckMProkaryote)
            cont  = 100 * cont / length(CheckMProkaryote)

            }

        completeness  = c(completeness,  round(comp,2))
        contamination = c(contamination, round(cont,2))

        MINCONSPERC_ASIG16  = 0.6
        MINCONSPERC_TOTAL16 = 0.3
        tax = rep('Unclassified',ncol(SQM$contigs$tax))
        good_ranks = c()
        for(rank in colnames(SQM$contigs$tax))
            {
            taxa_hits = sort(table(SQM$contigs$tax[contigs,rank]),decreasing=T)
            total = sum(taxa_hits)
            class_taxa_hits = taxa_hits[!grepl('^Unclassified|in NCBI|No CDS|bacterium$', names(taxa_hits))]
            total_class = sum(class_taxa_hits)
            # Skip this level if we don't fulfill the conditions
            if(!length(class_taxa_hits))                              { break }
            if(!class_taxa_hits[1]/total_class >= MINCONSPERC_ASIG16) { break }
            if(!class_taxa_hits[1]/total >= MINCONSPERC_TOTAL16)      { break }
            # Otherwise we try to use this for the taxonomy
            good_ranks = c(good_ranks, rank)
            LCA_tax = names(class_taxa_hits)[1]
            taxstr = SQM$misc$tax_names_long[[rank]][LCA_tax]
            tax = gsub('^k_|^p_|^c_|^o_|^f_|^g_|^s_', '', unlist(strsplit(taxstr, split=';')))
            missing_ranks = ncol(SQM$contigs$tax) - length(tax)
            if(missing_ranks)
                {
                unclass_str = sprintf('Unclassified %s', LCA_tax)
                for(i in 1:missing_ranks) { tax = c(tax, unclass_str) }
                }
            }
        taxonomy[b,] = tax

        if(!length(good_ranks))
            {
            dis = 0
        } else if(length(contigs)==1)
            {
            dis = 0
        } else
            {
            ### And get the disparity
            # Get the last classified rank
            last_rank = good_ranks[length(good_ranks)]
            # Get the consensus taxonomy at that rank
            last_tax  = taxonomy[b,last_rank]
            # Get the contigs taxonomy at that rank
            contigs_tax = SQM$contigs$tax[contigs,last_rank]
            contigs_tax = contigs_tax[!grepl('^Unclassified|in NCBI|No CDS', contigs_tax)]
            dis = sum(contigs_tax != last_tax) / length(contigs_tax)
            }
        disparity = c(disparity, round(dis,3))
        }
    resDF = data.frame(Method = SQM$bins$table[bins,'Method'], `Num contigs` = nContigs, `GC perc` = GC, `Tax 16S` = tax16S,
                       Disparity = disparity, Completeness = completeness, Contamination = contamination, row.names = bins, check.names = F)
    return(list(table = resDF, tax = taxonomy))
    }


get.bin.abunds = function(SQM, track_unmapped=FALSE)
    {
    res = list()
    x                = aggregate(SQM$contigs$abund, by=list(SQM$contigs$bins[,1]), FUN=sum)
    rownames(x)      = x[,1]
    x                = x[rownames(SQM$bin$table),-1,drop=F]
    nobin            = colSums(SQM$contigs$abund) - colSums(x)
    if(sum(nobin)>0) { x['No_bin',] = nobin }
    if(track_unmapped) { x['Unmapped',] = SQM$total_reads - colSums(x) }

    res[['abund']]   = as.matrix(x)
    res[['percent']] = 100 * t(t(res[['abund']]) / SQM$total_reads)

    x = aggregate(SQM$contigs$bases, by=list(SQM$contigs$bins[,1]), FUN=sum)
    rownames(x)      = x[,1]
    x = x[rownames(SQM$bin$table),-1,drop=F]
    res[['bases']]   = as.matrix(x)

    l = aggregate(SQM$contigs$table$Length, by=list(SQM$contigs$bins[,1]), FUN=sum)
    n = l[,1]; l = l[,-1]; names(l) = n
    l = l[rownames(SQM$bin$table)]
    res[['length']]  = l
    res[['cov']]     = res[['bases']] / res[['length']]
    res[['cpm']]     = round(t(t(res[['cov']]) / (SQM$total_reads / 1000000)), 6)
    res[['cov']]     = round(res[['cov']], 2)
    return(res)
    }

