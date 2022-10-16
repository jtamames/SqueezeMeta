#' Export results in tabular format
#'
#' This function is a wrapper for R's write.table function.
#' @param table vector, matrix or data.frame. The table to be written.
#' @param output_name character. Name of the output file.
#' @examples
#' data(Hadza)
#' Hadza.iron = subsetFun(Hadza, "iron")
#' # Write the taxonomic distribution at the genus level of all the genes related to iron.
#' exportTable(Hadza.iron$taxa$genus$percent, "Hadza.ironGenes.genus.tsv")
#' # Now write the distribution of the different iron-related COGs (Clusters of Orthologous Groups) across samples.
#' exportTable(Hadza.iron$functions$COG$tpm, "Hadza.ironGenes.COG.tsv")
#' # Now write all the information contained in the ORF table.
#' exportTable(Hadza.iron$orfs$table, "Hadza.ironGenes.orftable.tsv")
#' @export
exportTable = function(table, output_name)
    {
    write.table(table, output_name, col.names=NA, sep='\t', quote=F)
    }

