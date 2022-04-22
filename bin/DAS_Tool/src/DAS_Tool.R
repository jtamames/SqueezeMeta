#
# DAS Tool for genome-resolved metagenomics
# by Christian MK Sieber (csieber@lbl.gov)
#
# Please use DAS_Tool.sh in your installation folder to run DAS Tool.
# DAS Tool is available at https://github.com/cmks/DAS_Tool
# Please cite https://doi.org/10.1101/107789
#

library(methods) # FPS: 2-04-19 Explicitly load this bc Rscript does not load it on startup (but R console does). What were the odds?
packages <- list('doMC','data.table','ggplot2','DASTool')
foo <- lapply(packages,function(x){suppressPackageStartupMessages(library(x,character.only=TRUE))})

#disable warnings:
options(warn=-1)

#
#Parse command line args:
args <- commandArgs(trailingOnly = TRUE)
# 1.bins $scaffoldstobins
scaffolds_to_bins <- as.character(args[1])
# 2.labels $binlabels
bin_set_labels <- as.character(args[2])
# 3.contigs $contigs
assembly <-  as.character(args[3])
# 4.bac_scg $bscg
bac_scg_matrix <- as.character(args[4])
# 5.arc_scg $ascg
arc_scg_matrix <- as.character(args[5])
# 6.out $outputbasename
output_basename <- as.character(args[6])
# 7.score_threshold $score_threshold
score_threshold <- as.numeric(args[7])
# 8.threads $threads
threads <- as.numeric(args[8])
# 9.workingdir $thisdir
workingdir <- as.character(args[9])
# 10.length_table $length_table
length_table <- as.character(args[10])
# 11. write_bin_evals $write_bin_evals
write_bin_evals <- as.logical(args[11])
# 12.create_plots $create_plots
create_plots <- as.logical(args[12])
# 13.debug $debug
debug <- as.logical(args[13])
# 14.b $b
b <- as.numeric(args[14])
# 15.c $c
c <- as.numeric(args[15])

options(show.error.messages=FALSE)

if(debug){
  #enable warnings:
  options(warn=0)
  options(show.error.messages=T)
}

# check DAS Tool package version:
if(!packageVersion("DASTool")>=numeric_version("1.1.1")){
  cat('ERROR: DAS_Tool R-package (version 1.1.1) is not installed\n')
  cat('Please install the current version of DAS_Tool using:\n')
  cat('$ cd DAS_Tool_installation_directory\n')
  cat('$ R CMD INSTALL package/DASTool_1.1.1.tar.gz\n')
  cat('Or read the documentation for more detailed instructions\n')
  quit()
}

if(threads>24) { threads <- 24 } # Limit thread usage to segfaults in data.table
library(data.table) # the bug happened in_pick, which uses data.table. import data.table here so we can use setDTthreads. 
setDTthreads(threads)
registerDoMC(threads)
   
# cat('running DAS Tool with ', threads, ' threads\n', sep = '')
a <- 1
use_Nfifty <- T


if(score_threshold > 1){
  cat('WARNING: score_threshold is set to ', score_threshold, '. Should be between 0 and 1. Setting score_threshold to default (0.5).\n', sep = '')
  score_threshold <- 0.5
}

setwd(workingdir)

bin_evaluations <- cherry_pick(scaffolds_to_bins=scaffolds_to_bins,
                               bin_set_labels=bin_set_labels,
                               bac_scg_matrix=bac_scg_matrix,
                               arc_scg_matrix=arc_scg_matrix,
                               score_threshold=score_threshold,
                               a=a,b=b,c=c,
                               output_basename=output_basename,
                               length_table=length_table,
                               use_N50=use_Nfifty,
                               write_bin_evals=write_bin_evals,
                               debug=debug)

if(create_plots && class(bin_evaluations)=='list'){
  DER_Plot(bin_evaluations,output_basename)
}

