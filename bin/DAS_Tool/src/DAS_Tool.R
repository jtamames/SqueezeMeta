#!/usr/bin/env Rscript

#
# DAS Tool for genome-resolved metagenomics
# by Christian MK Sieber (cmksieber@gmail.com)
#
# Please cite https://doi.org/10.1101/107789
#
#
# DAS Tool Copyright (c) 2017, The Regents of the University of California, through Lawrence
# Berkeley National Laboratory (subject to receipt of any required approvals from the U.S.
# Dept. of Energy).  All rights reserved.
#
# If you have questions about your rights to use or distribute this software, please contact
# Berkeley Lab's Innovation and Partnerships department at IPO@lbl.gov referring to
# "DAS Tool (2017-024)".
#
# NOTICE.  This software was developed under funding from the U.S. Department of Energy.  As
# such, the U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce,
# prepare derivative works, and perform publicly and display publicly. The U.S. Government
# is granted for itself and others acting on its behalf a paid-up, nonexclusive,
# irrevocable, worldwide license in the Software to reproduce, prepare derivative works,
# distribute copies to the public, perform publicly and display publicly, and to permit
# others to do so.
#

doc <- "DAS Tool

Usage:
  DAS_Tool [options] -i <contig2bin> -c <contigs_fasta> -o <outputbasename>
  DAS_Tool -i <contig2bin> -c <contigs_fasta> -o <outputbasename> [--labels=<labels>] [--proteins=<proteins_fasta>] [--threads=<threads>] [--search_engine=<search_engine>] [--score_threshold=<score_threshold>] [--dbDirectory=<dbDirectory> ] [--megabin_penalty=<megabin_penalty>] [--duplicate_penalty=<duplicate_penalty>] [--max_iter_post_threshold=<max_iter>] [--write_bin_evals] [--write_bins] [--write_unbinned] [--resume] [--debug]
  DAS_Tool [--version]
  DAS_Tool [--help]

Options:
   -i --bins=<contig2bin>                   Comma separated list of tab separated contigs to bin tables.
   -c --contigs=<contigs>                   Contigs in fasta format.
   -o --outputbasename=<outputbasename>     Basename of output files.
   -l --labels=<labels>                     Comma separated list of binning prediction names.
   --search_engine=<search_engine>          Engine used for single copy gene identification (diamond/blastp/usearch) [default: diamond].
   -p --proteins=<proteins>                 Predicted proteins (optional) in prodigal fasta format (>contigID_geneNo).
                                            Gene prediction step will be skipped.
   --write_bin_evals                        Write evaluation of input bin sets.
   --write_bins                             Export bins as fasta files.
   --write_unbinned                         Write unbinned contigs.
   -t --threads=<threads>                   Number of threads to use [default: 1].
   --score_threshold=<score_threshold>      Score threshold until selection algorithm will keep selecting bins (0..1) [default: 0.5].
   --duplicate_penalty=<duplicate_penalty>  Penalty for duplicate single copy genes per bin (weight b).
                                            Only change if you know what you are doing (0..3) [default: 0.6].
   --megabin_penalty=<megabin_penalty>      Penalty for megabins (weight c). Only change if you know what you are doing (0..3) [default: 0.5].
   --max_iter_post_threshold=<max_iter>     Maximum number of iterations after reaching score threshold [default: 10].
   --dbDirectory=<dbDirectory>              Directory of single copy gene database [default: db].
   --resume                                 Use existing predicted single copy gene files from a previous run.
   --debug                                  Write debug information to log file.
   -v --version                             Print version number and exit.
   -h --help                                Show this.


Please cite: Sieber et al., 2018, Nature Microbiology (https://doi.org/10.1038/s41564-018-0171-1).
"

version <- 'DAS Tool 1.1.7\n'

if(length(commandArgs(trailingOnly = TRUE)) == 0L) {
   docopt:::help(doc)
   quit()
}

##
## Load packages, define functions
##

suppressMessages(suppressWarnings(library(data.table, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(magrittr, warn.conflicts = F, quietly = T)))
suppressMessages(suppressWarnings(library(docopt, warn.conflicts = F, quietly = T)))

cherry_pick <- function(binTab,scgTab,contigTab,output_basename,score_threshold,duplicate_penalty,megabin_penalty,max_iter_post_threshold,write_unbinned,write_bin_evals,logFile){
   
   # thresholds:
   # score_threshold <- max(score_threshold,-42)
   internal_ratio_threshold <- 0.0
   internal_score_threshold <- min(0.0,score_threshold)
   post_thresh_iter <- 10
   
   # set keys
   setkey(binTab,'contig_id')
   setkey(scgTab,'contig_id')
   setkey(contigTab,'contig_id')
   
   # join tables
   scgTab <- scgTab[ contig_id %in% binTab[, contig_id]]
   binTabScg <- binTab[scgTab, allow.cartesian = T]
   binTabContig <- contigTab[binTab]
   
   # score bins
   binTabEval <- score_bins(bin_tab_scg = binTabScg, bin_tab_contig = binTabContig, b = duplicate_penalty,c = megabin_penalty)
   
   if(write_bin_evals){
      write.table(binTabEval[,.(bin=bin_id,
                                bin_set=binner_name,
                                unique_SCGs=uniqueSCG,
                                redundant_SCGs=duplicatedSCG,
                                SCG_set=protein_set,
                                size=binSize,
                                contigs=nContig,
                                N50=contigN50,
                                bin_score=score,
                                SCG_completeness=round(completeness,2)*100,
                                SCG_redundancy=round(contamination,2)*100)],
                  paste(output_basename,'_allBins.eval',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
   }
   
   if(nrow(binTabEval) == 0 || max(binTabEval[,score]) < score_threshold){
      write.log(message = paste0('No bins with bin-score >', score_threshold, ' found. Adjust score_threshold to report bins with lower quality.'),
                filename = logFile,
                append = T,
                write_to_file = T,
                type = 'stop')
   }
   
   max_score <- max(binTabEval[,score])
   
   append <- F
   sub_bins <- c()
   post_thresh_iter <- 0
   
   while(max_score > internal_score_threshold && post_thresh_iter <= max_iter_post_threshold ){
      
      # identify highest scoring bin
      topBin <- binTabEval[ 1 ]
      
      # identify contigs of highest scoring bin
      topContig2Bin <- binTabContig[bin_id == topBin[,bin_id],.(contig_id,bin_id)]
      
      if(max_score >= score_threshold && max(binTabEval[,completeness]) > internal_ratio_threshold){
         
         # write contig2bins of highest scoring bin
         fwrite(topContig2Bin,paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
         
         # write summary of highest scoring bin
         fwrite(topBin[,.(bin=bin_id,
                          bin_set=binner_name,
                          unique_SCGs=uniqueSCG,
                          redundant_SCGs=duplicatedSCG,
                          SCG_set=protein_set,
                          size=binSize,
                          contigs=nContig,
                          N50=contigN50,
                          bin_score=score,
                          SCG_completeness=round(completeness,2)*100,
                          SCG_redundancy=round(contamination,2)*100)],
                paste0(output_basename,'_DASTool_summary.tsv'),sep='\t',col.names=(!append),row.names = F,quote=F,append = append)
         
         append <- T
      }else{
        post_thresh_iter <- post_thresh_iter+1
        
         # write unbinned contigs if bin-score is below threshold
         if(write_unbinned){
            topContig2Bin[,bin_id:='unbinned'] %>% 
               fwrite(paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
            
         }
      }
      
      # remove contigs of highest scoring bin
      affected_bins <- unique(binTabContig[ contig_id %in% topContig2Bin[,contig_id], bin_id ])
      binTabContig[ bin_id %in% affected_bins, bin_id:= paste(gsub('_sub$','',bin_id),'sub',sep='_')]
      binTabContig <- binTabContig[ !contig_id %in% topContig2Bin[,contig_id] ]
      binTabScg[ bin_id %in% affected_bins, bin_id:= paste(gsub('_sub$','',bin_id),'sub',sep='_')]
      binTabScg <- binTabScg[ !contig_id %in% topContig2Bin[,contig_id] ]
      
      # update scores of remaining bins
      # binTabEval <- score_bins(bin_tab_scg = binTabScg, bin_tab_contig = binTabContig, b = duplicate_penalty,c = megabin_penalty)
      # fwrite(binTabEval,paste0(output_basename,'_DASTool_binTabEval_OLD','_iter_',iter,'.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)

      binTabEval <- rbind(binTabEval[! bin_id %in% affected_bins],
                          score_bins(bin_tab_scg = binTabScg[ bin_id %in% paste(gsub('_sub$','',affected_bins),'sub',sep='_')],
                                   bin_tab_contig = binTabContig[ bin_id %in% paste(gsub('_sub$','',affected_bins),'sub',sep='_')],
                                   b = duplicate_penalty,
                                   c = megabin_penalty)) %>%
      .[ order(score,contigN50,binSize,decreasing = T)]
      
      if(nrow(binTabEval) == 0){
         max_score <- -Inf
      }else{
         max_score <- max(binTabEval[,score])
      }
   }
   
   # write remaining unbinned contigs
   if(write_unbinned){
      binTabContig[,.(contig_id)] %>% 
         .[,bin_id:='unbinned'] %>% 
         fwrite(paste0(output_basename,'_DASTool_contig2bin.tsv'),sep='\t',col.names=F,row.names = F,quote=F,append = append)
      
   }
}


calc_N50 <- function(contig_id,contig_length){
   
   seq_len <- data.table(contig_id,contig_length) %>% 
      .[!duplicated(contig_id)] %>% 
      .[order(contig_length,decreasing = T)] %>% 
      .[,contig_length] %>% 
     as.numeric()
   
   N50 <- seq_len[cumsum(seq_len) > sum(seq_len)/2][1] 
   
   return(N50)
}


calc_bins_size <- function(contig_id,contig_length){
   
   binSize <- data.table(contig_id,contig_length) %>%
      .[!duplicated(contig_id)] %>% 
      .[,contig_length] %>% 
      sum() %>% 
     as.numeric()
   
   return(binSize)
}


calc_duplicates <- function(protein_names){
   tmp <- table(protein_names)
   
   return(max(length(tmp[tmp>1]),0,na.rm = T))
}


calc_bin_score <- function(b,c,protein_set_size,uniqueSCG,duplicatedSCG,sumSCG,additionalSCG){
   
   bin_score <- (uniqueSCG/protein_set_size) - b*(duplicatedSCG / uniqueSCG) - c*(additionalSCG /  protein_set_size)
   
   return(bin_score)
}


score_bins <- function(bin_tab_scg,bin_tab_contig,b=.6,c=.5){
  
   bin_tab_scg_eval <- bin_tab_scg[,.(uniqueSCG=length(unique(protein_name)),
                                      duplicatedSCG=calc_duplicates(protein_name),
                                      sumSCG=.N),by=c('bin_id','protein_set','binner_name','protein_set_size')] %>% 
      .[,additionalSCG:= (sumSCG - uniqueSCG)] %>% 
      setkey(bin_id,binner_name)
    
   bin_tab_contig_eval <- bin_tab_contig[,.(binSize=as.numeric(calc_bins_size(contig_id,contig_length)),
                                            contigN50=as.numeric(calc_N50(contig_id,contig_length)),
                                            nContig=as.numeric(.N)),by=c('bin_id','binner_name')] %>% 
      setkey(bin_id,binner_name)
   
   bin_tab_eval <- bin_tab_contig_eval[bin_tab_scg_eval] %>% 
      .[,score:=calc_bin_score(b=b,c=c,protein_set_size,uniqueSCG,duplicatedSCG,sumSCG,additionalSCG)] %>% 
      .[,completeness:=uniqueSCG/protein_set_size] %>% 
      .[,contamination:=duplicatedSCG/protein_set_size] %>% 
      .[ order(score,contigN50,binSize,decreasing = T)] %>% 
      .[.[, .I[which.max(score)], by=bin_id]$V1]
   
   return(bin_tab_eval)
}

exit <- function() {
   invokeRestart("abort") 
} 

write.log <- function(message,append=T,filename,write_to_file=T,type='none'){
   
   if( 'character' %in% class(message)){
      if(write_to_file){
         cat(message,'\n',file=filename,sep=' ',append = append)
      }
      if(type == 'cat'){
         cat(message,'\n')
      }
      if(type == 'warning'){
         cat('Warning: ',message,'\n')
      }
      if(type == 'stop'){
         cat('Error: ',message,'\n')
         exit()
      }
   }else if("data.frame" %in% class(message)){
      write.table(message,file=filename,sep='\t',col.names=F,row.names = F,quote=F,append = append)
      cat('\n',file=filename,sep=' ',append = append)
   }
}


##
## Run DAS Tool
##

arguments <- docopt(doc, version = version)

# Existence and permissions of output directory
if(! dir.exists(dirname(arguments$outputbasename))){
   write.log(paste0('Output directory does not exist. Attempting to create: ', dirname(arguments$outputbasename),'\n'), 
             filename = '',append = T,write_to_file = F,type = 'warning')
   system(paste0('mkdir -p ',dirname(arguments$outputbasename)))
   
   if(! dir.exists(dirname(arguments$outputbasename))){
      write.log(paste0('Failed to create output directory: ', dirname(arguments$outputbasename)), 
                filename = '',append = T,write_to_file = F,type = 'stop')
   }
}

# Init log file
logFile <- paste0(arguments$outputbasename,'_DASTool.log')
write.log(gsub('\n','',version),filename = logFile,append = F,write_to_file = T,type = 'cat')
write.log(paste(Sys.time()),filename = logFile,append = T,write_to_file = T)
write.log('\nParameters:',filename = logFile,append = T,write_to_file = T)
argsTab <- data.table(opt=names(arguments),args=arguments) %>% 
   .[grepl('^--',opt)]
write.log(argsTab,filename = logFile,append = T,write_to_file = T)

threads <- max(1,as.numeric(arguments$threads))

# Check options
## search engine %in% diamond, usearch, blastp?
searchEngine <- tolower(arguments$search_engine)
if(!searchEngine %in% c('diamond', 'usearch', 'blastp')){
   write.log(paste0('Unknown argument for --search_engine: ',arguments$search_engine,'\n',
                    'Defaulting to diamond'),filename = logFile,append = T,write_to_file = T,type = 'warning')
   searchEngine <- 'diamond'
}

# Check dependencies
dependencies <- Sys.which(c("prodigal", "diamond", "pullseq", "ruby", "usearch", "blastp")) %>% 
   data.table(dependency=names(.),path=.)

min_dependencies <- c("pullseq", "ruby",searchEngine)
if(is.null(arguments$proteins)){
   min_dependencies <- c(min_dependencies,'prodigal')
}

if(any(dependencies[ path == '', dependency ] %in% min_dependencies)){
   
   write.log(paste0('Cannot find dependencies: ', 
                    paste0(dependencies[ path == '' ] %>% .[dependency %in% min_dependencies, dependency] ,collapse = ', ')),filename = logFile,append = T,write_to_file = T,type = 'stop')
}

write.log('Dependencies:',filename = logFile,append = T,write_to_file = T)
write.log(dependencies,filename = logFile,append = T,write_to_file = T,type = 'cat')

searchEngine <- gsub('blastp','blast',searchEngine) # for SCG-finder blast-script

# Get script directory:
cmdArgs <- commandArgs(trailingOnly = FALSE)
scriptDir <- dirname(normalizePath(sub("--file=", "", cmdArgs[grep("--file=", cmdArgs)])))
scriptDir <- gsub('\\/src$','',scriptDir)

# Set default db directory if not defined:
dbDirectory <- ifelse(arguments$dbDirectory == 'db',
                       paste0(scriptDir,'/','db/'),
                       arguments$dbDirectory)

# add trailing slash if needed
if (!endsWith(dbDirectory, '/')) {
    dbDirectory <- paste0(dbDirectory, '/')
}

##
## Check files and directories
##
binSets <- unlist(strsplit(arguments$bins,','))

# Input files existence
## Contig2bin tables:
for(binSet in binSets){
   if(!file.exists(binSet)){
      write.log(paste0('File does not exist: ',binSet),filename = logFile,append = T,write_to_file = T,type = 'stop')
   }
}

## Assembly:
if(!file.exists(arguments$contigs)){
   write.log(paste0('File does not exist: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
}
if(file.size(arguments$contigs) == 0){
   write.log(paste0('File is empty: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
}
tmpfile <- file(arguments$contigs)
if(summary( tmpfile )$class == 'gzfile'){
   write.log(paste0('Unsupported file: Assembly file is compressed. Please extract: ',arguments$contigs),filename = logFile,append = T,write_to_file = T,type = 'stop')
}
close(tmpfile)

## Proteins:
if(!is.null(arguments$proteins) && !file.exists(arguments$proteins)){
   write.log(paste0('File does not exist: ',arguments$proteins),filename = logFile,append = T,write_to_file = T,type = 'stop')
}

# Check bin set labels
binSetLabels <- basename(binSets)
if( !is.null(arguments$labels) ){
   binSetLabels <- unlist(strsplit(arguments$labels,','))

   if(length(binSetLabels) != length(binSets)){
      write.log(paste0('Number of bin sets (', length(binSets),') is different from number of labels (', length(binSetLabels),')'),
                filename = logFile,append = T,write_to_file = T,type = 'warning')
      binSetLabels <- basename(binSets)
   }
}
if(any(duplicated(binSetLabels))){
   if( !is.null(arguments$labels) ){
      write.log('Non-unique bin set labels given. Renaming bin set labels.', filename = logFile,append = T,write_to_file = T,type = 'warning')
   }
   # write.log('Creating bin set labels', filename = logFile,append = T,write_to_file = T,type = 'cat')
   binSetLabels <- paste0(binSetLabels,'.binner.',sprintf("%02d",c(1:length(binSets))))
}

binSetToLabel <- data.table(inputNo=c(1:length(binSets)),binSetLabel=binSetLabels,binSet=binSets)

# Check contig2bin files:
binTab <- data.table()
for(i in 1:nrow(binSetToLabel)){
   if(file.size(binSetToLabel[i,binSet]) == 0){
      ## Warn if contig2bin table is empty, but keep going:
      write.log(paste0('File is empty: ',binSetToLabel[i,binSet]), filename = logFile,append = T,write_to_file = T,type = 'warning')
   }else{
      ## Read contig2bin table and concatenate:
      dt <- fread(binSetToLabel[i,binSet],
                  sep='\t',
                  header=F,
                  colClasses=c(V1='character',V2='character'),
                  col.names=c('contig_id','bin_id')) %>% 
         .[,binner_name:=binSetToLabel[i,binSetLabel]] %>% 
         .[,contig_id:=gsub(' .*','',contig_id)] # trim contig-ids after first whitespace
      binTab <- rbind(binTab,dt)
   }
}

if(arguments$debug){
   write.log('binTab:', filename = logFile,append = T,write_to_file = T,type = 'debug')
   write.log(binTab, filename = logFile,append = T,write_to_file = T,type = 'debug')
}

## Stop if all contig2bin tables are empty:
if(nrow(binTab) == 0){
   write.log('No bins provided. Check contig2bin files.', filename = logFile,append = T,write_to_file = T,type = 'stop')
}

## Check for duplicate bin-IDs
if(nrow(unique(binTab[,.(binner_name,bin_id)])) > length(unique(binTab[,bin_id]))){
   binTab[,bin_id:=paste0(binner_name,'_',bin_id)]
   write.log(paste0('Non-unique bin-ids given. Renaming bin-ids.'), filename = logFile,append = T,write_to_file = T,type = 'warning')
}


# Existence of database directory
if( !file.exists(paste0(dbDirectory,'bac.all.faa')) || 
    !file.exists(paste0(dbDirectory,'arc.all.faa')) || 
    !file.exists(paste0(dbDirectory,'bac.scg.faa')) || 
    !file.exists(paste0(dbDirectory,'arc.scg.faa')) || 
    !file.exists(paste0(dbDirectory,'bac.scg.lookup')) || 
    !file.exists(paste0(dbDirectory,'arc.scg.lookup'))){
   write.log(paste0('Database directory does not exist or is incomplete:', dbDirectory,'\n
                  Attempting to extract ',scriptDir,'/db.zip'), 
             filename = logFile,append = T,write_to_file = T,type = 'warning')
   if(file.exists(paste0(scriptDir,'/db.zip'))){
      system(paste0('unzip -o ',scriptDir,'/db.zip -d ',dbDirectory))
   }else{
      write.log(paste0('File does not exist: ',scriptDir,'/db.zip'), 
                filename = logFile,append = T,write_to_file = T,type = 'stop')
   }
}

##
## Calculate contig lengths
##
write.log('Analyzing assembly',filename = logFile,append = T,write_to_file = T,type = 'cat')
system(paste0('bash ', scriptDir,'/src/contig_length.sh ',
              arguments$contigs,' ',
              arguments$outputbasename,'.seqlength'))

contigTab <- fread(paste0(arguments$outputbasename,'.seqlength'),
                   sep='\t',
                   header=F,
                   colClasses=c(V1='character',V2='integer'),
                   col.names=c('contig_id','contig_length')) %>% 
   .[,contig_id:=gsub(' .*','',contig_id)]

if(arguments$debug){
   write.log('Assembly stats:', filename = logFile,append = T,write_to_file = T,type = 'debug')
   assemblyStatsTab <- data.table(contigs=nrow(contigTab),
                                  size=sum(contigTab[,contig_length]),
                                  median=as.double(median(contigTab[,contig_length])),
                                  min=min(contigTab[,contig_length]),
                                  max=max(contigTab[,contig_length]),
                                  N50=calc_N50(contigTab[,contig_id],contigTab[,contig_length]))
   write.log(colnames(assemblyStatsTab),filename = logFile,append = T,write_to_file = T,type = 'debug')
   write.log(assemblyStatsTab,filename = logFile,append = T,write_to_file = T,type = 'debug')
}

## Check for duplicate contig-IDs
if( any(duplicated(contigTab[,contig_id]))){
   duplicatedContigs <- contigTab[duplicated(contigTab[,contig_id])]
   write.log(paste0('Error: Duplcate contig-IDs found in assembly: ', nrow(duplicatedContigs)), 
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   write.log(paste0(paste0(head(duplicatedContigs[,contig_id],5),collapse = ',\n'), ifelse(nrow(duplicatedContigs) >5,',...','')),
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   exit()
}

## Check if bin-set contigs are in assembly:
missingContigs <- binTab[!contig_id %in% contigTab[,contig_id]]
if( nrow(missingContigs) > 0){
   write.log(paste0('Error: Contigs of contig2bin files not found in assembly: ', nrow(missingContigs)), 
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   write.log(paste0(paste0(head(missingContigs[,contig_id],5),collapse = ',\n'), ifelse(nrow(missingContigs) >5,',...','')),
             filename = logFile,append = T,write_to_file = T,type = 'cat')
   exit()
}

##
## Predict single copy genes
##

# Run prodigal
if(!is.null(arguments$proteins)){
   proteins <- arguments$proteins
   write.log(paste0('Skipping gene prediction\n\t using protein fasta file: ',proteins),filename = logFile,append = T,write_to_file = T,type = 'cat')
}else{
   write.log('Predicting genes',filename = logFile,append = T,write_to_file = T,type = 'cat')
   
   proteins <- paste0(arguments$outputbasename,'_proteins.faa')

   if(threads == 1){
      system(paste('prodigal -i ',arguments$contigs,
                   ' -a ',proteins,
                   ' -p meta -m -q > /dev/null 2>&1'))
   }else{
      system(paste0('bash ', scriptDir,'/src/prodigal_parallel.sh ',
                    arguments$contigs,' ',
                    proteins, ' ',
                    threads))
   }
}
if(!file.exists(proteins) || file.size(proteins) == 0){
   write.log('Gene prediction failed.',filename = logFile,append = T,write_to_file = T,type = 'stop')
}


# Identify single copy genes
scgTab <- data.table()
## Predict bacterial SCGs
bacteria_scg_out <- paste0(proteins,'.bacteria.scg')
if(!arguments$resume || !file.exists(bacteria_scg_out)){
   write.log(paste0('Annotating single copy genes using ',searchEngine),filename = logFile,append = T,write_to_file = T,type = 'cat')
   system(paste0('ruby ', scriptDir,'/src/scg_blank_',searchEngine,'.rb ',
                 searchEngine,' ',
                 proteins,' ',
                 dbDirectory,'/bac.all.faa ',
                 dbDirectory,'/bac.scg.faa ',
                 dbDirectory,'/bac.scg.lookup ',
                 threads,
                 ifelse(arguments$debug,paste0(' 2>&1 | tee -a ',logFile, ' > /dev/null 2>&1'),' > /dev/null 2>&1')))
   if(file.exists(paste0(proteins,'.scg')) && file.size(paste0(proteins,'.scg')) > 0){
      system(paste0('mv ',proteins,'.scg ',bacteria_scg_out))
      scgTab <- rbind(scgTab,
                      fread(paste0(proteins,'.bacteria.scg'),header=F,sep='\t',
                            colClasses=c(V1='character',V2='character'),
                            col.names=c('protein_id','protein_name')) %>% 
                         .[,protein_set:='bacteria'] %>% 
                         .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                         .[,protein_set_size:=51])
   }else{
      write.log(paste0('No SCGs detected for SCG set: bacteria'), 
                filename = logFile,append = T,write_to_file = T,type = 'warning')
   }
}else{
   write.log(paste0('Skipping SCG prediction for SCG set: bacteria','\n\t using ',bacteria_scg_out),
             filename = logFile,type = 'cat')
   if(file.exists(bacteria_scg_out) && file.size(bacteria_scg_out) > 0){
      scgTab <- rbind(scgTab,
                      fread(paste0(proteins,'.bacteria.scg'),header=F,sep='\t',
                            colClasses=c(V1='character',V2='character'),
                            col.names=c('protein_id','protein_name')) %>% 
                         .[,protein_set:='bacteria'] %>% 
                         .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                         .[,protein_set_size:=51])
   }else{
      write.log(paste0('Empty table for SCG set: bacteria: ', bacteria_scg_out), 
                filename = logFile,append = T,write_to_file = T,type = 'warning')
   }
}

## Predict archaeal SCGs
archaea_scg_out <- paste0(proteins,'.archaea.scg')
if(!arguments$resume || !file.exists(archaea_scg_out)){
   system(paste0('ruby ', scriptDir,'/src/scg_blank_',searchEngine,'.rb ',
                 searchEngine,' ',
                 proteins,' ',
                 dbDirectory,'/arc.all.faa ',
                 dbDirectory,'/arc.scg.faa ',
                 dbDirectory,'/arc.scg.lookup ',
                 threads,
                 ifelse(arguments$debug,paste0(' 2>&1 | tee -a ',logFile, ' > /dev/null 2>&1'),' > /dev/null 2>&1')))
   if(file.exists(paste0(proteins,'.scg')) && file.size(paste0(proteins,'.scg')) > 0){
      system(paste0('mv ',proteins,'.scg ',archaea_scg_out))
      scgTab <- rbind(scgTab,
                      fread(paste0(proteins,'.archaea.scg'),header=F,sep='\t',
                            colClasses=c(V1='character',V2='character'),
                            col.names=c('protein_id','protein_name')) %>% 
                        .[,protein_set:='archaea'] %>% 
                        .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                        .[,protein_set_size:=38])
   }else{
      write.log(paste0('No SCGs detected for SCG set: archaea'), 
                filename = logFile,append = T,write_to_file = T,type = 'warning')
   }
}else{
   write.log(paste0('Skipping SCG prediction for SCG set: archaea','\n\t using ',archaea_scg_out),
             filename = logFile,type = 'cat')
   if(file.exists(archaea_scg_out) && file.size(archaea_scg_out) > 0){
      scgTab <- rbind(scgTab,
                      fread(paste0(proteins,'.archaea.scg'),header=F,sep='\t',
                            colClasses=c(V1='character',V2='character'),
                            col.names=c('protein_id','protein_name')) %>% 
                         .[,protein_set:='archaea'] %>% 
                         .[,contig_id:=gsub('_[0-9]+$','',protein_id)] %>% 
                         .[,protein_set_size:=38])
   }else{
      write.log(paste0('Empty table for SCG set: archaea: ', archaea_scg_out), 
                filename = logFile,append = T,write_to_file = T,type = 'warning')
   }
}

## Stop if no single copy genes were predicted:
if(nrow(scgTab) == 0){
   write.log('No single copy genes predicted', 
             filename = logFile,append = T,write_to_file = T,type = 'stop')
}

if(arguments$debug){
   write.log('\nscgTab:', filename = logFile,append = T,write_to_file = T,type = 'debug')
   write.log(scgTab, filename = logFile,append = T,write_to_file = T,type = 'debug')
}

##
## Run bin selection
##
write.log('Dereplicating, aggregating, and scoring bins',filename = logFile,append = T,write_to_file = T,type = 'cat')
cherry_pick(binTab=binTab,
            scgTab=scgTab,
            contigTab=contigTab,
            output_basename=arguments$outputbasename,
            score_threshold=as.numeric(arguments$score_threshold),
            duplicate_penalty=as.numeric(arguments$duplicate_penalty),
            megabin_penalty=as.numeric(arguments$megabin_penalty),
            max_iter_post_threshold=as.numeric(arguments$max_iter_post_threshold),
            write_unbinned=arguments$write_unbinned,
            write_bin_evals=arguments$write_bin_evals,
            logFile=logFile)


##
## Extract bins
##
if(arguments$write_bins){
   write.log('Writing bins',filename = logFile,append = T,write_to_file = T,type = 'cat')
   binDir <- paste0(arguments$outputbasename,'_DASTool_bins')
   if(!dir.exists(binDir)){
      system(paste0('mkdir ',binDir))
   }
   system(paste0('bash ', scriptDir,'/src/extract_bins.sh ',arguments$outputbasename,'_DASTool_contig2bin.tsv ',arguments$contigs,' ',binDir))
}

