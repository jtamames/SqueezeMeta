### A generic table works internally using the data.frame or data.table engines, but has the same interface as a data.frame

library(data.table)

#' @importFrom utils read.table tail
#' @importFrom zip unzip
read.generic.table.zip = function(project_path, file_path, engine = 'data.frame', ...)
    {
    extra_args = list(...)
    if(engine == 'data.frame')
        {
        extra_args = extra_args[names(extra_args) %in% names(formals(read.table))]
        extra_args$file = open.conn.zip(project_path, file_path) # see open.conn.zip in extra_methods.R
        res = do.call(read.table, extra_args)
    } else if(engine == 'data.table')
        {
        extra_args = list(...)
        extra_args = extra_args[names(extra_args) %in% names(formals(data.table::fread))]
        zipmode = endsWith(project_path, '.zip')
        if(!zipmode)
            {
            extra_args$file = sprintf('%s/%s', project_path, file_path)
            res = do.call(data.table::fread, extra_args)
        } else
            {
            unzip(project_path, file_path, exdir = tempdir(), junkpaths = TRUE) # junkpaths=TRUE so the file is extracted directly into tempdir()
            f = tail(unlist(strsplit(file_path, split = '/')), 1)               #  instead of creating "results" or "intermediate" directories
            f = sprintf('%s/%s', tempdir(), f)                                  #  mostly so we can remove it easily after using it
            extra_args$file = f
            res = do.call(data.table::fread, extra_args)
            unlink(f)
            }
        }
    return(generic.table(res))
    }


generic.table = function(x)
    {
    if('data.table' %in% class(x))
        {
        rowns = x[[1]]
        x = x[,-1]
        rownames(x) = rowns
        class(x) = c('generic.data.table', class(x))
    }else if ('data.frame' %in% class(x))
        {
        class(x) = c('generic.data.frame', class(x))
        }
    return(x)
    }


#' @export
#' @noRd
#' @importFrom data.table data.table
`[.generic.data.table` = function(x,i,j, drop = TRUE)
    {
    
    class(x) = class(x)[class(x) != 'generic.data.table']

    if(missing(i))
        {
        i = '' # just so we can print it in tests
        si = ''
        names = rownames(x)
    }else
        {
        if(typeof(i) == 'character') { names = i; i = data.table::chmatch(i, rownames(x)); # turn into indices and keep the order of i
        } else { names = rownames(x)[i] }
        si = 'i'
        }

    if(missing(j)) { j = 1:(ncol(x)) }


    else if(typeof(j) %in% c('numeric', 'integer', 'double'))
        {
        if(max(abs(j)) > (ncol(x))) { stop('Undefined columns selected') }
        }

    sj = 'j'

    expression = sprintf('x[%s,%s, with=FALSE]', si, sj)
    res = eval(parse(text=expression))
    if(!is.null(dim(res))){ rownames(res) = names; class(res) = c('generic.data.table', class(res));  }
    if(drop & ncol(res)==1) { res = res[[1]] }
    return (res)
    }


#' @export
#' @noRd
print.generic.data.table = function(x, ...)
    {
    print.data.frame(x)
    }


#' @export
#' @noRd
dimnames.generic.data.table = function(x)
    {
    return(list(row.names(x), names(x)))
    }



# We don't export it since having an "as.data.table" method breaks it during install (works with any other name). Anyways the native "as.data.table" from the data.table package should work.
as.data.table.generic.data.table = function(x)
    {
    res = x
    class(res) = class(res)[class(res) != 'generic.data.table']
    return(res)
    }


#' @export
#' @noRd
#' @importFrom utils getFromNamespace
as.data.frame.generic.data.table = function(x, ...)
    {
    res = getFromNamespace('as.data.frame.data.table', 'data.table')(x)
    rownames(res) = rownames(x)
    return(res)
    }


#' @export
#' @noRd
#' @importFrom utils getFromNamespace
as.matrix.generic.data.table = function(x, ...)
    {
    res = getFromNamespace('as.matrix.data.table', 'data.table')(x)
    rownames(res) = rownames(x)
    return(res)
    }


#' @export
#' @noRd
rbind.generic.data.table = function(...)
    {
    res = data.table::rbindlist(list(...))
    rownames(res) = unlist(sapply(list(...), rownames))
    class(res) = c('generic.data.table', class(res))
    return(res)
    }


#cbind.generic.data.table = function(...)
#    {
#    inter = Reduce(intersect, lapply(list(...), rownames))
#    uni = Reduce(union, lapply(list(...), rownames))
#    if(!identical(inter, uni)) { stop('The input tables do not have the same row names') }
#    res = data.table:::cbind(...) # This will not work until data.table reintroduces cbind or adds cbindlist!
#    rownames(res) = rownames(list(...)[[1]])
#    class(res) = c('generic.data.table', class(res))
#    return(res)
#    }
