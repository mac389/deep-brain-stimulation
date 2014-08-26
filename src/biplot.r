
data <- read.csv('pain-v.tsv',sep='\t')
pc <- princomp(data,scale=TRUE)

png(filename='../results/pain-biplot.png')
biplot(pc,xlabs=rep(c('x'),nrow(data)))
dev.off()

png(filename='../results/pain-scree.png')
plot(pc,main=NULL)
dev.off()

