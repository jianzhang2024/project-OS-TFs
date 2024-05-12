library(pheatmap)         
expFile=" " 
clusterFile=" " 
cliFile=" " 

setwd(" ")    

exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>15,">15","<=15"))
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

data=data[order(data$TFcluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
data=t(data[,1:ncol(exp)])

bioCol=c("#3969AC","#E73F74", "#7F3C8D", "#11A579", "#F2B701",  "#80BA5A", "#A5AA99")
ann_colors=list()
prgCluCol=bioCol[1:length(levels(factor(Type$TFcluster)))]
names(prgCluCol)=levels(factor(Type$TFcluster))
ann_colors[["TFcluster"]]=prgCluCol

pdf("heatmap.pdf", width=9, height=9)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("#3969AC",4), "white", rep("#E73F74",4)))(100),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize=18,
         fontsize_row=18,
         fontsize_col=18)
dev.off()
