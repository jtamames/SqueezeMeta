read.namedvector = function(file)
    {
    ta = read.table(file, header=T, row.names=1, as.is=T)
    res = ta[,1]
    names(res) = rownames(ta)
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
