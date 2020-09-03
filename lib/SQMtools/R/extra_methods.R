read.namedvector = function(file, engine = 'data.frame')
    {
    if(!engine %in% c('data.frame', 'data.table')) { stop('Engine must be "data.frame" or "data.table"') }
    if(engine == 'data.frame')
        {
        ta = read.table(file, header=T, row.names=1, as.is=T)
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


#' Return a vector with the row-wise minima of a matrix or dataframe.
#' @export
rowMins = function(table)
    {
    res = do.call(pmin, as.data.frame(table))
    names(res) = rownames(table)
    return(res)
    }


#' Return a vector with the row-wise maxima of a matrix or dataframe.
#' @export
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
    return(cbind(m1[allRows,,drop=F], m2[allRows,,drop=F]))
    }


SQMtoSQMlite = function(SQM) # untested and unused
    {
    if(!class(SQM) %in% c('SQM', 'SQMlite')) { stop('This function only accepts SQM objects') }
    SQMlite             = list()
    SQMlite$taxa        = SQM$taxa
    SQMlite$functions   = SQM$functions
    SQMlite$total_reads = SQM$total_reads
    SQMlite$misc        = SQM$misc
    class(SQMlite)      = 'SQMlite'
    return(SQMlite)
    }
