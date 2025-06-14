% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combineSQM.R
\name{combineSQM}
\alias{combineSQM}
\title{Combine several SQM objects}
\usage{
combineSQM(
  ...,
  tax_source = "orfs",
  trusted_functions_only = FALSE,
  ignore_unclassified_functions = FALSE,
  rescale_tpm = TRUE,
  rescale_copy_number = TRUE,
  recalculate_bin_stats = TRUE
)
}
\arguments{
\item{...}{an arbitrary number of SQM objects. Alternatively, a single list containing an arbitrary number of SQM objects.}

\item{tax_source}{character. Features used for calculating aggregated abundances at the different taxonomic ranks. Either \code{"orfs"} or \code{"contigs"} (default \code{"orfs"}). If the objects being combined contain a subset of taxa or bins, this parameter can be set to \code{TRUE}.}

\item{trusted_functions_only}{logical. If \code{TRUE}, only highly trusted functional annotations (best hit + best average) will be considered when generating aggregated function tables. If \code{FALSE}, best hit annotations will be used (default \code{FALSE}).}

\item{ignore_unclassified_functions}{logical. If \code{FALSE}, ORFs with no functional classification will be aggregated together into an "Unclassified" category. If \code{TRUE}, they will be ignored (default \code{FALSE}).}

\item{rescale_tpm}{logical. If \code{TRUE}, TPMs for KEGGs, COGs, and PFAMs will be recalculated (so that the TPMs in the subset actually add up to 1 million). Otherwise, per-function TPMs will be calculated by aggregating the TPMs of the ORFs annotated with that function, and will thus keep the scaling present in the parent object (default \code{TRUE}).}

\item{rescale_copy_number}{logical. If \code{TRUE}, copy numbers with be recalculated using the median single-copy gene coverages in the subset. Otherwise, single-copy gene coverages will be taken from the parent object. By default it is set to \code{TRUE}, which means that the returned copy numbers will represent the average copy number per function \emph{in the genomes of the selected bins or contigs}. If any SQM objects that are being combined contain a functional subset rather than a contig/bins subset, this parameter should be set to \code{FALSE}.}

\item{recalculate_bin_stats}{logical. If \code{TRUE}, bin abundance, quality and taxonomy are recalculated based on the contigs present in the subsetted object (default \code{TRUE}).}
}
\value{
A SQM or SQMbunch object
}
\description{
Combine an arbitrary number of SQM objects into a single SQM object (if the input objects contain the same samples, i.e. they come from the same SqueezeMeta run) or a single SQMbunch object. For combining results from sqm_reads.pl or sqm_longreads.pl please check \code{\link{combineSQMlite}}. The parameters below (other than ...) will take only effect if the input objects contain the same samples. Otherwise the input objects will be taken as they are, with no recalculation of taxonomy, function or rescaling,
}
\examples{
data(Hadza)
# Select Carbohydrate metabolism ORFs in Bacteroidota,
#  and Amino acid metabolism ORFs in Proteobacteria
bact = subsetTax(Hadza, "phylum", "Bacteroidota")
bact.carb = subsetFun(bact, "Carbohydrate metabolism")
baci = subsetTax(Hadza, "phylum", "Bacillota")
baci.amins = subsetFun(baci, "Amino acid metabolism")
bact.carb_proteo.amins = combineSQM(bact.carb, baci.amins, rescale_copy_number=FALSE)
}
\seealso{
\code{\link{subsetFun}}, \code{\link{subsetTax}}, \code{\link{combineSQMlite}}
}
