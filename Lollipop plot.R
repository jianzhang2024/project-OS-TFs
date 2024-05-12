inputFile=" "  
setwd(" ")  

data = read.table(inputFile, header=T, sep="\t", check.names=F)

p.col = c("#F6A97A","#FA7876","#EA4F88","#C0369D","#872CA2")
fcolor = function(x,p.col){
  color = ifelse(x>0.04,p.col[1],ifelse(x>0.03,p.col[2],ifelse(x>0.02,p.col[3],
                ifelse(x>0.01,p.col[4], p.col[5])
                )))
  return(color)
}

p.cex = seq(3, 8, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<2,p.cex[1],ifelse(x<2.5,p.cex[2],ifelse(x<3,p.cex[3],
              ifelse(x<3.5,p.cex[4],p.cex[5]))))
  return(cex)
}

points.color = fcolor(x=data$pvalue,p.col=p.col)
data$points.color = points.color

points.cex = fcex(x=data$HR)
data$points.cex = points.cex
data=data[order(data$HR),]

xlim = ceiling(max(abs(data$HR))*10)/10   
pdf(file="Lollipop.pdf", width=9, height=9)  
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,8,1,12),cex.axis=2,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="log(HR value)",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)

segments(x0=data$HR,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=2)

points(x=data$HR,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)

text(par('usr')[1],1:nrow(data),data$id,adj=1,xpd=T,cex=2)
axis(1,tick=F)

par(mar=c(0,1,1,0))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(1,2,3,4,5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="|log(HR value)|")

par(mar=c(0,1,1,8),cex.axis=2,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(0.05,0.04,0.03,0.02,0.01,0),tick=F)
dev.off()
