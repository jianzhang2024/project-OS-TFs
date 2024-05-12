library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
R.utils::setOption("clusterProfiler.download.method",'auto') 
pvalueFilter=0.05 
qvalueFilter=0.05 

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd(" ") 
rt=read.table(" ", header=F, sep="\t", check.names=F) 

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]    

kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=30, color=colorSel)
dev.off()

pdf(file="bubble.pdf", width = 6, height = 2.8)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", 
        label_format=70)
bub=bub + scale_color_continuous(low="#3969AC",high="#E73F74")
print(bub)
dev.off()
