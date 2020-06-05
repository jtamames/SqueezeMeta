### A generic table works internally using the data.frame or data.table engines, but has the same interface as a data.frame

generic.table = function(x)
    {
    if('data.table' %in% class(x))
        {
        setcolorder(x, c(colnames(x)[colnames(x) != colnames(x)[1]], colnames(x)[1])) # Put the first column (row names) at the end note that this happens in place!
        colnames(x)[ncol(x)] = 'rn'
        rownames(x) = x$rn
        class(x) = c('generic.data.table', class(x))
    }else if ('data,frame' %in% class(x))
        {
        class(x) = c('generic.data.frame', class(x))
        }
    return(x)
    }


as.data.table.generic.data.table = function(x)
    {
    res = x
    class(res) = class(res)[class(res) != 'generic.data.table']
    return(res)
    }


`[.generic.data.table` = function(x,i,j, drop = T)
    {
    # add support para indexados negativos

    class(x) = class(x)[class(x) != 'generic.data.table']

    if(missing(i))
        {
        i = '' # just so we can print it in tests
        si = ''
        names = rownames(x)
    }else
        {
        if(typeof(i) == 'character') { i = chmatch(i, x$rn) } # turn into indices and keep the order of i
        si = 'i'
        names = x[i,'rn']
        }

    if(missing(j)) { j = 1:(ncol(x)) }
    else if(typeof(j) %in% c('numeric', 'integer', 'double'))
        {
        if(max(abs(j)) > (ncol(x)-1)) { stop('Undefined columns selected') }
        if(max(j) > 0) { j = c(j, ncol(x)) }
        }
    else if(typeof(j) == 'character') { j = c(j, 'rn') }
    else if(typeof(j) == 'logical') { j = c(j, TRUE) }

    sj = 'j'


    expression = sprintf('x[%s,%s, with=FALSE]', si, sj)
    res = eval(parse(text=expression))
    if (ncol(res) == 2 & drop == TRUE) { res = unlist(res[,1]); names(res) = NULL
    } else { rownames(res) = names; class(res) = c('generic.data.table', class(res));  }
    return (res)
    }



print.generic.data.table = function(x)
    {
    `[` = data.table:::`[.data.table` # is this needed?
    res = x[,1:ncol(x),drop=F]
    rownames(res) = rownames(x)
    print.data.frame(res)
    }


dim.generic.data.table = function(x)
    {
    return(data.table:::dim.data.table(x) - c(0,1))
    }



dimnames.generic.data.table = function(x)
    {
    return(list(row.names(x), names(x)[-(ncol(x)+1)]))
    }


as.matrix.generic.data.table = function(x)
    {
    `[` = data.table:::`[.data.table`
    res = data.table:::as.matrix.data.table(x[,-(ncol(x)+1)])
    rownames(res) = rownames(x)
    return(res)
    }


colnames = function(x) UseMethod("colnames") # For making a non-generic R function generic and add a custom class method without breaking stuff
colnames.generic.data.table = function(x) dimnames(x)[[2]] # https://gist.github.com/datalove/88f5a24758b2d8356e32
colnames.default = base::colnames


