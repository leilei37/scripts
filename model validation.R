#sort by heritability
h1<-sample(which(pvea[,2]>0.01&pvea[,2]<0.1),50)
h2<-sample(which(pvea[,2]>0.1&pvea[,2]<0.2),50)
h3<-sample(which(pvea[,2]>0.2&pvea[,2]<0.3),50)
h4<-sample(which(pvea[,2]>0.3&pvea[,2]<0.4),50)
h5<-sample(which(pvea[,2]>0.4&pvea[,2]<0.5),50)
h6<-sample(which(pvea[,2]>0.5&pvea[,2]<0.6),50)
h7<-sample(which(pvea[,2]>0.6&pvea[,2]<0.7),50)
h8<-sample(which(pvea[,2]>0.7&pvea[,2]<0.8),50)
h9<-sample(which(pvea[,2]>0.8&pvea[,2]<0.9),50)
h10<-sample(which(pvea[,2]>0.9&pvea[,2]<1),50)
##5 fold cross validation
#n=length(Y.rtr)
folds=sample(1:5,size=298,replace=T)
#which(folds==5)
accuracyestimate<-function(h3=h3,accuracymatrixh3=accuracymatrixh3){
accuracymatrixh3<-matrix(ncol = length(h3),nrow = 6)
for (j in h3) {
  y<-phall[,j]
for(i in 1:max(folds)){
  tst=which(folds==i)
  G11=G[-tst,-tst] # genomic relationships in the training data
  G21=G[tst,-tst]
  yTRN=y[setdiff(1:298,tst)]
  fm=BGLR(y=as.numeric(yTRN),ETA=list(list(K=G11,model='RKHS')),nIter=2000,burnIn=1000)
  yHat_4=fm$mu+as.vector(G21%*%solve(G11)%*%fm$ETA[[1]]$u)
  cor1<-cor.test(as.numeric(y[tst]),yHat_4)
  corr<-cor1$estimate
  colnum<-which(h3==j)
  accuracymatrixh3[i+1,colnum]<-corr
  accuracymatrixh3[1,colnum]<-j
}
}
}
accuracyestimate(h3=h4,accuracymatrixh3=accuracymatrixh4)
