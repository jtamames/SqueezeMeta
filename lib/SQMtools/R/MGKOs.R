#' Single Copy Phylogenetic Marker Genes from Sunagawa's group (KOs)
#'
#' Lists of Single Copy Phylogenetic Marker Genes.
#' These are useful for transforming coverages or tpms into copy numbers.
#' This is an alternative way of normalizing data in order to be able to
#' compare functional profiles in samples with different sequencing depths.
#'
#' @docType data
#'
#' @usage data(MGKOs)
#'
#' @format Character vector with the KEGG identifiers for 10 Single Copy Phylogenetic Marker Genes.
#'
#' @keywords datasets
#'
#' @references Salazar, G \emph{et al.} (2019). 
#' Gene Expression Changes and Community Turnover Differentially Shape the Global Ocean Metatranscriptome
#' \emph{Cell} \bold{179}:1068-1083.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31730850/}{PubMed}).
#'
#' @seealso \code{\link[MGKOs]{MGOGs}} for an equivalent list using OGs instead of KOs; \code{\link[USiCGs.R]{USiCGs}} for an alternative set of single copy genes, and for examples on how to generate copy numbers.
"MGKOs"
