library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library("stringr")

pvalueFilter=0.05 
qvalueFilter=0.05
showNum=5
setwd(" ")
rt=read.table(" ",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
   colorSel="pvalue"
}
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, 
            pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)
if(nrow(GO)<30){
   showNum=nrow(GO)
}

pdf(file="GO_barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,
            split="ONTOLOGY",color = colorSel) + 
  facet_grid(ONTOLOGY~., scale='free')+
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()
pdf(file="GO_bubble.pdf",width = 6.5,height =6.5)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",
            split="ONTOLOGY", color = colorSel) + 
  facet_grid(ONTOLOGY~., scale='free')+
  scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=45))
bub=bub + scale_color_continuous(low="#3969AC",high="#E73F74")
print(bub)
dev.off()

pdf(file="GO_cnet.pdf",width = 10,height = 5)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk, showCategory = 10, categorySize="pvalue",
         circular = TRUE,colorEdge = TRUE,
         cex_label_category=0.65,cex_label_gene=0.6)
dev.off()
pdf(file="GO_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = 20,cex_label_category=0.65,
         color = "pvalue",layout ="nicely" )
dev.off()
