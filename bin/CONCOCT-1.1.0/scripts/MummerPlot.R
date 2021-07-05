library(ggplot2)
library(gridExtra)
library(reshape)

contigs<-read.table("newannotations.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
genome<-subset(contigs,OutbreakGenome == 'Y')
ggplot(genome, aes(x=as.numeric(OutbreakGenome_Pos)/1000000, y=Cov, colour=factor(Cluster))) +
    geom_point(alpha=0.5) +
    theme_bw() +
    scale_x_continuous("Genome position (megabases)", breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5)) +
    scale_y_continuous("Total coverage depth")
ggsave("Ecoli_Alignment.pdf", width=10, height=6)
