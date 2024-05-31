#' summary method for class SQMbunch
#'
#' Computes different statistics of the data contained in the SQMbunch object.
#' @param object SQMbunch object to be summarized.
#' @param ... Additional parameters (ignored).
#' @return A list of summary statistics.
#' @export
summary.SQMbunch = function(object, ...)
    {

    SQM = object # so that CRAN is happy

    if(!inherits(SQM, 'SQMbunch')) { stop('The first argument must be a SQM object') }

    res = summary.SQMlite(SQM)

    return(res)
    }


#' @export
#' @noRd
print.summary.SQMbunch = function(x, ...)
    {
    print.summary.SQMlite(x)
    }
