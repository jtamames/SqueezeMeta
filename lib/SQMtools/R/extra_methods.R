#' @importFrom utils read.table
read.namedvector = function(file, engine = 'data.frame')
    { 
    if(!engine %in% c('data.frame', 'data.table')) { stop('Engine must be "data.frame" or "data.table"') }
    if(engine == 'data.frame')
        {
        ta = read.table(file, header=TRUE, row.names=1, as.is=TRUE)
        res = ta[,1]
        names(res) = rownames(ta)
    } else if (engine == 'data.table')
        {
        ta = data.table::fread(file, sep='\t')
        res = unlist(ta[,2])
        names(res) = unlist(ta[,1])
        }
    return(res)
    }

#' @importFrom utils tail unzip
read.namedvector.zip = function(project_path, file_path, engine = 'data.frame')
    {
    zipmode = endsWith(project_path, '.zip')
    if(!zipmode)
        {
        res = read.namedvector(sprintf('%s/%s', project_path, file_path), engine = engine)
    } else # since data.table::fread can't read from connections we need to do uncompress the file before reading it
        {
        unzip(project_path, file_path, exdir = tempdir(), junkpaths = TRUE) # junkpaths=TRUE so the file is extracted directly into tempdir()
        f = tail(unlist(strsplit(file_path, split = '/')), 1)               #  instead of creating "results" or "intermediate" directories
        f = sprintf('%s/%s', tempdir(), f)                                  #  mostly so we can remove it easily after using it
        res = read.namedvector(f, engine = engine)
        unlink(f)
        }
    return(res)
    }


#' Return a vector with the row-wise minima of a matrix or dataframe.
#' @export
#' @param table matrix or dataframe.
#' @return a vector with the row-wise minima.
rowMins = function(table)
    {
    res = do.call(pmin, as.data.frame(table))
    names(res) = rownames(table)
    return(res)
    }


#' Return a vector with the row-wise maxima of a matrix or dataframe.
#' @export
#' @param table matrix or dataframe.
#' @return a vector with the row-wise maxima.
rowMaxs = function(table)
    {
    res = do.call(pmax, as.data.frame(table))
    names(res) = rownames(table)
    return(res)
    }


merge.numeric.matrices = function(m1, m2)
    {
    notIn1 = setdiff(rownames(m2), rownames(m1))
    m1 = rbind(m1, matrix(0, nrow=length(notIn1), ncol=ncol(m1), dimnames=list(notIn1, colnames(m1))))
    notIn2 = setdiff(rownames(m1), rownames(m2))
    m2 = rbind(m2, matrix(0, nrow=length(notIn2), ncol=ncol(m2), dimnames=list(notIn2, colnames(m2))))
    allRows = sort(rownames(m1))
    return(cbind(m1[allRows,,drop=FALSE], m2[allRows,,drop=FALSE]))
    }


named.unique = function(v)
    {
    return(v[!duplicated(v)])
    }


SQMtoSQMlite = function(SQM) # untested and unused
    {
    if(!inherits(SQM, c('SQM', 'SQMlite'))) { stop('This function only accepts SQM objects') }
    SQMlite             = list()
    SQMlite$taxa        = SQM$taxa
    SQMlite$functions   = SQM$functions
    SQMlite$total_reads = SQM$total_reads
    SQMlite$misc        = SQM$misc
    class(SQMlite)      = 'SQMlite'
    return(SQMlite)
    }


check.samples = function(SQM, samples)
    {
    if(!is.null(samples))
        {
        missing_samples = setdiff(samples, SQM$misc$samples)
        if(length(missing_samples) > 0)
            {
            str = paste(missing_samples, collapse = '", "')
            stop(sprintf('Samples "%s" are not present in this SQM object', str))
            }
        }
    }


get.bin.stats = function(SQM)
    {
    bins = unique(SQM$contigs$bins[,'DAS'])
    bins = sort(bins[bins!='No_bin'])
    GC            = c()
    nContigs      = c()
    disparity     = c()
    completeness  = c()
    contamination = c()
    tax16S        = c()
    taxonomy      = matrix(NA, nrow = length(bins), ncol = ncol(SQM$contigs$tax),
                   dimnames = list(bins, colnames(SQM$contigs$tax)))
    for(b in bins) {
        contigs  = rownames(SQM$contigs$bins)[SQM$contigs$bins[,'DAS']==b]
        nContigs = c(nContigs, length(contigs))
        GC       = c(GC, mean(SQM$contigs$table[contigs,'GC perc']))
        orfs     = rownames(SQM$orfs$table)[SQM$orfs$table$`Contig ID` %in% contigs]
        orfs16S  = orfs[orfs %in% names(SQM$orfs$tax16S)]
        tax16S   = c(tax16S, paste(SQM$orfs$tax16S[orfs16S], collapse='|'))

        marker_count = table(unlist(SQM$orfs$markers[orfs]))

        comp = 0
        cont = 0
        
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
            if(!length(class_taxa_hits))                             { break }
            if(!class_taxa_hits[1]/total_class > MINCONSPERC_ASIG16) { break }
            if(!class_taxa_hits[1]/total > MINCONSPERC_TOTAL16)      { break }
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
	} else if(length(contigs)>10000)
	    {
            dis = NA
            warning(sprintf('Bin %s has more than 10000 contigs, skipping disparity calculation.', b))
        } else
            {
            ### And get the disparity
            # Get a table with the classified ranks for this bin
            tax_table = SQM$contigs$tax[contigs,good_ranks,drop=F]
            # Mark unclassified as NA so they won't be used for comparison
            tax_table[grepl('^Unclassified|in NCBI|No CDS|bacterium$', tax_table)] = NA
            # Get all possible combinations of contigs
            combis = combn(contigs,2)
            chims = matrix(FALSE, nrow = length(good_ranks), ncol = ncol(combis), dimnames = list(good_ranks, 1:ncol(combis)))
    
            for(i in 1:nrow(chims))
                {
                t1 = tax_table[combis[1,], i]
                t2 = tax_table[combis[2,], i]
                chims[i,] = t1 != t2
                # A pair is chimeric (TRUE) if it was chimeric in a higher rank
                # Also it will be unclassified (NA) if it was unclassified in a higher rank
                # Both are achieved with the `|` operator
                if(i>1) { chims[i,] = chims[i,] | chims[i-1,] }
                }
            diffs = chims[nrow(chims),]
            diffs = diffs[!is.na(diffs)] # ignore Unclassified
            dis = sum(diffs) / length(diffs)
            }
        disparity = c(disparity, round(dis,3))
        }
    resDF = data.frame(Method = SQM$bins$table[bins,'Method'], `Num contigs` = nContigs, `GC perc` = GC, `Tax 16S` = tax16S,
                       Disparity = disparity, Completeness = completeness, Contamination = contamination, `Strain het` = NA, row.names = bins, check.names = F)
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


file.exists.zip = function(project_path, file_path)
    {
    zipmode = endsWith(project_path, '.zip')
    file_exists = TRUE
    if(zipmode) { f = unz(project_path, file_path) } else { f = sprintf('%s/%s', project_path, file_path) }
    if(zipmode)
        {
        res=try(suppressWarnings(scan(f, what = 'character', quiet = TRUE, n = 1)), silent=TRUE)
        if(inherits(res, 'try-error')) { file_exists = FALSE }
        close(f)
        }
    else
        {
        if(!file.exists(f)) { file_exists = FALSE }
        }
    return(file_exists)
    }


open.conn.zip = function(project_path, file_path)
    # remember to open/close conns explicitly later after calling this function if needed!
    #  some functions (eg read.table) will open and close the connection by themselves, so there's nothing to be done
    #  others (eg scan) seem to open it but not close it
    #  have fun exploring this amazingly consistent behaviour!
    {
    zipmode = endsWith(project_path, '.zip')
    if(zipmode) { conn = unz(project_path, file_path) } else { conn = file(sprintf('%s/%s', project_path, file_path)) }
    }


#' @importFrom utils unzip
list.files.zip = function(project_path, dir_path)
    {
    zipmode = endsWith(project_path, '.zip')
    if(!zipmode)
        {
        res = list.files(sprintf('%s/%s', project_path, dir_path))
        }
    else
        {
        all_files = unzip(project_path, list = TRUE)[[1]]
        res = all_files[startsWith(all_files, dir_path)]
        res = gsub(sprintf('^%s', dir_path), '', res)
        res = gsub('^/', '', res)
        res = res[res!='']
        res = unique(sapply(strsplit(res, '/'), `[`, 1))
        }
    return(res)
    }
