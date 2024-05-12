library(limma)
library(ggplot2)

expFile=" "
clusterFile=" "
setwd(" ") 

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=t(data)

data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="newTab.xls", quote=F, sep="\t")

cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
TFcluster=as.vector(cluster[,1])

bioCol=c("#3969AC", "#E73F74", "#80BA5A", "#A5AA99")
prgCluCol=bioCol[1:length(levels(factor(TFcluster)))]

PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], TFcluster=TFcluster)
PCA.mean=aggregate(PCA[,1:2], list(TFcluster=PCA$TFcluster), mean)
pdf(file="PCA.pdf", width=2.8, height=2.6)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = TFcluster)) +
	scale_colour_manual(name="TFcluster", values =prgCluCol)+
    theme_bw()+
    theme(legend.position=c(0.8,0.8), plot.margin=unit(rep(0.5,4),'lines'))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
