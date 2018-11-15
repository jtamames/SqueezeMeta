library("gplots");
drawmarker <- function(filename, outname)
{
filename;
outname;
	read.table(filename, sep="\t", header=TRUE) -> t
	t2 = t[,seq(4, length(t))]
	t2[is.na(t2)] = 0
	t3 = as.matrix(t2)
	per = paste("Percent = ", format((length(t3[t3==1]) / length(t3) * 100), digits = 2, nsmall=2), "%")
	c = gray((max(t3) + 1) : 0 / (max(t3) + 1))
	pdf(outname);
	heatmap.2(t3, margins=c(5,10), dendrogram="none", trace="none", density.info=c("none"), Rowv = FALSE, Colv = FALSE, col=c)
	text(0.5, 0.9, per)
	dev.off()
}

args <- commandArgs(TRUE);
drawmarker(args[1], args[2])

