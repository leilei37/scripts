install.packages('hglm')
data_rna@phdata<-data_rna@phdata[-rms,]
ks<-setdiff(1:298,rms)
gtdata1<-gtdata(data_rna)[ks,all1qtl]
gtdata2<-gtdata(data_rna)[ks,all2qtl]
gtdataa<-gtdata(data_rna)
data_rna@gtdata<-gtdataa
data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
#the heritability of eQTL in 1 tissue
data_rna@gtdata<-gtdata1
data1.gkin1 <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
#the heritability of eQTL in >1 tissue
data_rna@gtdata<-gtdata2
data1.gkin2 <- ibs(data_rna[, autosomal(data_rna)], weight="freq")

library(hglm)
hglm_model2 <- hglm(X = as.matrix(X0[na,]),y = y[na],Z = cbind(Z,Z2),RandC = c(ncol(Z),ncol(Z2)),calc.like = T)
c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d))
  return(Z)
}
Z1 <- c.z.hglm(kin = data1.gkin1)
Z1<-Z1[-rms,-rms]
Z2 <- c.z.hglm(kin = data1.gkin2)
Z2<-Z2[-rms,-rms]
X0<-matrix(1,nrow = 263,ncol = 1)
rms<-which(is.na(phdata(data_rna)[,3]))#remove na Samples
outputpve<-list()
for (i in 3:13) {
  y0<-as.numeric(as.character(phdata(data_rna)[,i]))
  hglm_model2 <- hglm(X0,y = y0,Z = cbind(Z1,Z2),RandC = c(ncol(Z1),ncol(Z2)),calc.like = T)
  outputpve[[i]][1:2]<-hglm_model2$varRanef/sum(c(hglm_model2$varRanef ,hglm_model2$varFix))
}
#convert list to dataframe
apve<-do.call(rbind.data.frame, outputpve)
colnames(apve)[2]<-'ConversedQTL'
colnames(apve)[1]<-'SpecificQTL'
apve<-cbind(colnames(phdata(data_rna))[3:13],apve)
#hglm_model2$varFix is residual的variance component
#hglm_model2$varRanef are variance components of 2 random effects
#heribility is equal to sum（c）
write.csv(apve,file = '~/Documents/apve.csv')
apve<-read.csv('~/Documents/apve.csv')
library(ggplot2)
p <- ggplot(apve,aes(x=reorder(Traits,X),y=PVE,fill=QTLtype))+geom_bar(position="dodge",stat="identity")
p+xlab("") + ylab("PVE") + labs(fill="")+theme_classic()+ scale_fill_manual(values = c('#7972A6','#F2B84B'))+
  theme(legend.position = 'bottom',axis.text.x=element_text(vjust=1,size=5,face = "bold"))
