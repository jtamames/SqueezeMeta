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
