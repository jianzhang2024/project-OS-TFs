library(limma)
library(survival)
library(dplyr)
library(ggplot2)
library(ggthemes)

expFile=" " 
cliFile=" " 
setwd(" ")

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]
data=log2(data+1) 
data1=t(data)
rownames(data1)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data1))

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365

sameSample=intersect(row.names(data1), row.names(cli))
data1=data1[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data1)

outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<1){
		sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

write.table(outTab,file="uniCox volcano.txt",sep="\t",row.names=F,quote=F)

data <- read.table("uniCox volcano.txt", header=TRUE)
data$HR=log(data$HR)

data$type[data$pvalue < 0.05 & 
            data$HR > 0] = "Positively correlated genes"
data$type[data$pvalue < 0.05 &
            data$HR < 0] = "Negatively correlated genes"

data$type[data$pvalue >= 0.05] = "Not significant genes"

p<-ggplot(data,aes(x=HR,y=-1*log10(pvalue),colour=type))+
  xlab("log(HR value)")+
  ylab("-log10(pvalue)")+
  geom_point(size=1.2,alpha=1.2)+
  scale_color_manual(values =c("#3969AC","#A5AA99","#E73F74"))+
  geom_hline(yintercept=1.30103,linetype=3)+
  theme_few()+
  theme(legend.position="right",
        legend.key.size = unit(0.5, "inches"),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) 

top_data <- data %>% group_by(type) %>% 
  top_n(n = 5, wt = -1*log(pvalue)) %>% 
  filter(type !='Not significant genes') 

p<-p+geom_text(data=top_data,
               aes(x=HR,y=-1*log10(pvalue),
                   label=as.character(id)),size=4) 
ggsave(p, file='volcano.pdf',width=5.8,height=2.8)
