pve<-c()
for (j in 1:2) {
  i<-which(colnames(phdata(data_rna))==ldlist[j])
  data1<-print(as.numeric(as.character(phdata(data_rna)[,i])),digits=6)
  #data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
  hy<- try(polygenic(data1,kinship.matrix =data1.gkin,data =data_rna),silent = T)
  pve<-append(pve,hy$esth2)
}   #calculate gene expression PVE 
library(GenABLE)
library(BGLR)
df.gkin <- ibs(df[, autosomal(df)], weight="freq")
G<-df.gkin       #make G matrix by GenABLE
rr1<-G
rr1[upper.tri(rr1)]<-t(rr1)[upper.tri(rr1)]
G<-rr1
tst<-c(299:648)
G11=G[-tst,-tst] # genomic relationships in the training data
G21=G[tst,-tst]
nowph<-phlb
new350<-matrix(nrow =350 ,ncol = length(nowph))
i<-3
for (i in 3:length(nowph)) {
  fm=BGLR(y=as.numeric(nowph[,i]),ETA=list(list(K=G11,model='RKHS')),nIter=2000,burnIn=1000)
  yHat_4=fm$mu+as.vector(G21%*%solve(G11)%*%fm$ETA[[1]]$u)
  new350[,i]<-yHat_4
}
save(new350,file = '')
