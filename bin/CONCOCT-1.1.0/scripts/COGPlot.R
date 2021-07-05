#!/usr/bin/Rscript

library(ggplot2)
library(reshape)

library(getopt)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'scgfile','s',1,"character",'ofile','o',1,"character"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE)); 
	q(status=1);
}

scgfile <- opt$scgfile


tab <- read.table(scgfile,header=TRUE,row.names=1)
ecogs <- tab[,3:ncol(tab)]

maxecogs <- max(ecogs)

sumCogs <- rowSums(ecogs)

ecogs$sum <- sumCogs

ecogs <- subset(ecogs,ecogs$sum > 0)

ecogs.order <- ecogs[order(ecogs$sum),]

names <- row.names(ecogs.order)

ecogs.order$sum <- NULL

N=nrow(ecogs.order)
alpha = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z')
slist = rep("AA",N)

i=1
j=1
n=1

while (n <= N) {
    slist[n] = (paste(alpha[i],alpha[j],sep=""))
    j = j + 1
    n = n + 1
    if (j > 26) {
        j = 1
        i = i + 1
    }
}

ecogs.order$id <- slist

mecogsorder <- melt(ecogs.order,id=("id"))

color <- c("ivory2", "#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027", "#2166ac",  "#4393c3", "#92c5de", "#d1e5f0")
if (maxecogs > 10) {
  color <- c("ivory2", "#1a9850", "#91cf60", "#d9ef8b", "#fee08b", "#fc8d59", "#d73027", "#2166ac",  "#4393c3", "#92c5de", "#d1e5f0", grey((maxecogs - 10):0/(1.2*(maxecogs-10))))
}
p <- ggplot(mecogsorder, aes(variable,id)) + geom_tile(aes(fill = as.factor(value)), colour = "white") + scale_fill_manual(name="COG counts",values = color)

pdf(opt$ofile)
p + scale_y_discrete(breaks=slist,labels=names) + theme(axis.text.x=element_text(size=7,angle=-90,vjust=0.5)) + ylab("Cluster")
dev.off()
