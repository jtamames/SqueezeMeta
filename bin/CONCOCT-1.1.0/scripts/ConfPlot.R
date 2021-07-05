#!/usr/bin/Rscript

#load libraries
library(gplots)
library(getopt)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'confile','c',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE)); 
	q(status=1);
}

confFile <- opt$confile

Conf <- read.csv(confFile,header=TRUE,row.names=1)
Conf.t <- t(Conf)
Conf.t <- Conf.t[rowSums(Conf.t) > 0,]
ConfP <- Conf.t/rowSums(Conf.t)

crp <- colorRampPalette(c("blue","red","orange","yellow"))(100)

pdf(opt$ofile)

heatmap.2 (as.matrix(t(ConfP)),col=crp,trace="none",dendrogram="none",Rowv=FALSE,lwid = c(1.25,4.5),lhei = c(1.25,4.5),cexRow=0.5)

dev.off()
