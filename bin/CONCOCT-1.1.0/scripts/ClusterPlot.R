#!/usr/bin/Rscript

#load libraries
library(ggplot2)
library(ellipse)
library(getopt)
library(grid)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'cfile','c',1,"character",'pcafile','p',1,"character",'mfile','m',1,"character",'proot','r',1,"character",'ofile','o',1,"character",'legend','l',0,"logical"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE)); 
	q(status=1);
}

clusterFile <- opt$cfile
pcaFile <- opt$pcafile
meanFile <- opt$mfile
pcaRoot <- opt$proot


PCA <- read.csv(pcaFile,header=TRUE,row.names=1)
Clusters <- read.csv(clusterFile,header=FALSE,row.names=1)
means <- read.csv(meanFile,header=FALSE)


PCA.df <- data.frame(x=PCA[,1],y=PCA[,2],c=Clusters$V2)
PCA.df$c <- factor(PCA.df$c)

df_ell <- data.frame()
for(c in levels(PCA.df$c))
{
	filename = sprintf("%s%s.csv",pcaRoot,c); 
	
	print(filename)
	
	temp <- read.csv(filename,header=FALSE)
	temp2 <- as.matrix(temp[1:2,1:2])
	
	elt <-  as.data.frame(ellipse(temp2,centre=c(means[strtoi(c) + 1,1],means[strtoi(c) + 1,2])))	
	eltc <- cbind(elt,group=c)
	
	df_ell <- rbind(df_ell,eltc)
	rm(temp)
	rm(temp2)
	
} 

colnames(df_ell)[1] <- "x"
colnames(df_ell)[2] <- "y"

colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00");

shapes <- c(15,16,17,18)

nC <- length(colours);
nS <- length(shapes);

nClust <- length(levels(PCA.df$c))

valuesC <- vector()
valuesS <- rep(16,nClust);

for(i in 1:nClust){
	valuesC[i] <- colours[i %% nC + 1] 
	valuesS[i] <- 15 + i %/% nC
}

print(valuesC);
print(valuesS);

# Order the factor levels
valuesS <- valuesS[as.integer(factor(PCA.df$c, levels = sort(unique(PCA.df$c))))]

pdf(opt$ofile)
theme_set(theme_bw(20))
p <- ggplot(data=PCA.df, aes(x=x, y=y,colour=c)) + geom_point(size=1.0, alpha=.3) + xlab("PCA1") + ylab("PCA2") + scale_colour_manual(values=valuesC) + scale_shape_manual(values=valuesS) + geom_path(data=df_ell, aes(x=x, y=y,colour=as.factor(group)), size=0.5, linetype=2)

if( !is.null(opt$legend)){ p + theme(legend.key.size = unit(0.3, "cm")) + guides(col = guide_legend(ncol = 2,override.aes = list(alpha = 1)))+ theme(legend.text=element_text(size=4));}else{p + theme(legend.position="none");}

dev.off()
