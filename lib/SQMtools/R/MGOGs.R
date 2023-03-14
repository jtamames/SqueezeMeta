#' Single Copy Phylogenetic Marker Genes from Sunagawa's group (OGs)
#'
#' Lists of Single Copy Phylogenetic Marker Genes.
#' These are useful for transforming coverages or tpms into copy numbers.
#' This is an alternative way of normalizing data in order to be able to
#' compare functional profiles in samples with different sequencing depths.
#'
#' @docType data
#'
#' @usage data(MGOGs)
#'
#' @format Character vector with the COG identifiers for 10 Single Copy Phylogenetic Marker Genes.
#'
#' @keywords datasets
#'
#' @references Salazar, G \emph{et al.} (2019). 
#' Gene Expression Changes and Community Turnover Differentially Shape the Global Ocean Metatranscriptome
#' \emph{Cell} \bold{179}:1068-1083.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/31730850/}{PubMed}).
#'
#' @seealso \code{\link{MGKOs}} for an equivalent list using KOs instead of OGs; \code{\link{USiCGs}} for an alternative set of single copy genes, and for examples on how to generate copy numbers.
"MGOGs"
