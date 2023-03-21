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
