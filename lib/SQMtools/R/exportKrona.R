#' Export the taxonomy of a SQM object into a krona Chart
#'
#' This function selects the most abundant functions across all samples in a SQM object and represents their abundances in a heatmap. Alternatively, a custom set of functions can be represented.
#'
#' Original code was kindly provided by Giuseppe D'Auria (dauria_giu@gva.es).
#'
#' @param SQM A SQM or SQMlite object.
#' @param output_name character. Name of the output file containing the Krona charts in html format (default \code{"<project_name>.krona.html")}.
#' @seealso \code{\link[plotTaxonomy]{plotTaxonomy}} for plotting the most abundant taxa of a SQM object.
#' @examples
#' data(Hadza)
#' exportKrona(Hadza)
#' @export
exportKrona = function(SQM, output_name = NA)
    {

    # Check params.
    if(!class(SQM) %in% c('SQM', 'SQMlite')) { stop('The first argument must be a SQM or a SQMlite object') }


    # Check that kronatools is present.
    ecode = system('ktImportText', ignore.stdout = T, ignore.stderr = T)
    if(ecode!=0) { stop('ktImportText (from KronaTools) is not present in your PATH. KronaTools can be downloaded from https://github.com/marbl/Krona') }

    # Should we use the default output name?
    if(is.na(output_name))
        {
        output_name = sprintf('%s.krona.html', SQM$misc$project_name)
	tempDir = 'kronaTemp'
	}
    tempDir = sprintf('%s/SQMtools.krona.temp', tempdir())
    system(sprintf('rm -r %s', tempDir)   , ignore.stdout = T, ignore.stderr = T)
    system(sprintf('mkdir %s', tempDir), ignore.stdout = T, ignore.stderr = T)

    # Prepare data.
    ta = SQM$taxa$species$abund
    taxFields = strsplit(SQM$misc$tax_names_long$species[rownames(ta)], split=';', fixed = T )
    taxFields = lapply(taxFields, function(x) sapply(strsplit(x, split='_', fixed=T), function(y) y[2]))
    taxFields = t(data.frame(taxFields))

    for(sample in SQM$misc$samples)
        {
        idx = which(ta[,sample]>0) # Taxa that are present in this sample.
	kronaTemp = cbind(ta[idx,sample], taxFields[idx,])
        write.table(kronaTemp, file = sprintf('%s/%s.tsv', tempDir, sample),
                    sep='\t', quote=F, col.names=F, row.names=F) # Write temp file for this sample
        }

    # Create krona report.
    system(sprintf('ktImportText -o %s %s/*.tsv', output_name, tempDir), ignore.stdout = T, ignore.stderr = T)
    }
