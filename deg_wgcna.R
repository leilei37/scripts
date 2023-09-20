i<-0
aHN<-c()
for (i in 0:7) {
  group<-c(i*6+1,i*6+2,i*6+3)
  group<-as.character(lapply(group,function(x){paste0('X',x)}))
  testmatrix<-exp[,group]
  rownames(testmatrix)<-exp[,1]
  u<-testmatrix
  library(edgeR)
  keep <- rowSums(cpm(u) > 0.5) >= 2
  u<-u[keep,]
  library(psych)
  png(filename = paste0('~/Downloads/dd/',i*6+1,'-',i*6+3,'cor.png'))
  pairs.panels(u,
               gap = 0,
               bg = c("grey", "blue","black"),
               pch=21,
               digits=3,
               cex.cor = 1)
  dev.off()
  aHN<-append(aHN,group)
}
aHNm<-exp[,aHN]
rownames(aHNm)<-exp[,1]
u<-aHNm
library(edgeR)
keep <- rowSums(cpm(aHNm) > 0.5) >= 2
aHNm<-aHNm[keep,]
#pca
dat_scale=scale(aHNm,scale=F)
dat_cor=cor(dat_scale)
#options(digits=4, scipen=4)
#kable(dat_cor)
dat_eigen=eigen(dat_cor)
dat_var=dat_eigen$values ## 相关系数矩阵的特征值
dat_var
pca_var=dat_var/sum(dat_var)
pca_var
pca_cvar=cumsum(dat_var)/sum(dat_var)
pca_cvar
library(ggplot2)
p=ggplot(,aes(x=1:12,y=pca_var))
p1=ggplot(,aes(x=1:12,y=pca_cvar))
#p+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
#p1+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
pca_vect= dat_eigen$vector  ## 相关系数矩阵的特征向量
loadings=data.frame(sweep(pca_vect,2,sqrt(pca_var),"*"))
rownames(loadings)=colnames(aHNm)

qplot(loadings[,1], loadings[,2], data = loadings, 
      color=aHN,
      xlab = 'PC1 77.67%',ylab = 'PC2 10.70%')+
  geom_text(aes(loadings[,1], loadings[,2],label = aHN2))+
  theme(legend.title=element_blank(),legend.key.size = unit(10, "pt"),axis.title.x = element_text(size = 10))
dev.off()
aHN1<-c()
for (i in 1:8) {
#  for (j in 1:3) {
  aHN1<-append(aHN1,lapply(i,function(x){rep(paste0('HN',x,'_R'),3)}))
}
aHN1<-as.character(unlist(aHN1))
aHN2<-c()
for (i in 1:24) {
  aHN2[i]<-paste0(aHN1[i],rep(1:3,8)[i])
}

i<-0
aHB<-c()
for (i in 0:7) {
  group<-c(i*6+4,i*6+5,i*6+6)
  group<-as.character(lapply(group,function(x){paste0('X',x)}))
  testmatrix<-exp[,group]
  rownames(testmatrix)<-exp[,1]
  u<-testmatrix
  library(edgeR)
  keep <- rowSums(cpm(u) > 0.5) >= 2
  u<-u[keep,]
  library(psych)
  png(filename = paste0('~/Downloads/dd/',i*6+4,'-',i*6+6,'cor.png'))
  pairs.panels(u,
               gap = 0,
               bg = c("grey", "blue","black"),
               pch=21,
               digits=3,
               cex.cor = 1)
  dev.off()
  aHB<-append(aHB,group)
}
aHBm<-exp[,aHB]
rownames(aHBm)<-exp[,1]
u<-aHBm
library(edgeR)
keep <- rowSums(cpm(aHBm) > 0.5) >= 2
aHBm<-aHBm[keep,]
#pca
dat_scale=scale(aHBm,scale=F)
dat_cor=cor(dat_scale)
#options(digits=4, scipen=4)
#kable(dat_cor)
dat_eigen=eigen(dat_cor)
dat_var=dat_eigen$values ## 相关系数矩阵的特征值
dat_var
pca_var=dat_var/sum(dat_var)
pca_var
pca_cvar=cumsum(dat_var)/sum(dat_var)
pca_cvar
library(ggplot2)
p=ggplot(,aes(x=1:12,y=pca_var))
p1=ggplot(,aes(x=1:12,y=pca_cvar))
#p+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
#p1+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
pca_vect= dat_eigen$vector  ## 相关系数矩阵的特征向量
loadings=data.frame(sweep(pca_vect,2,sqrt(pca_var),"*"))
rownames(loadings)=colnames(aHBm)

qplot(loadings[,1], loadings[,2], data = loadings, 
      color=aHB,
      xlab = 'PC1 79.90%',ylab = 'PC2 8.29%')+
  geom_text(aes(loadings[,1], loadings[,2],label = aHB2))+
  theme(legend.title=element_blank(),legend.key.size = unit(10, "pt"),axis.title.x = element_text(size = 10))
dev.off()
aHB1<-c()
for (i in 1:8) {
  #  for (j in 1:3) {
  aHB1<-append(aHB1,lapply(i,function(x){rep(paste0('H',x,'B_R'),3)}))
}
aHB1<-as.character(unlist(aHB1))
aHB2<-c()
for (i in 1:24) {
  aHB2[i]<-paste0(aHB1[i],rep(1:3,8)[i])
}
colnames(aHBm)<-aHB2
colnames(aHNm)<-aHN2

write.table(aHBm,'~/Downloads/dd/iseqqc/HBmatrix.txt',sep = '\t',quote = F)
write.table(aHNm,'~/Downloads/dd/iseqqc/HNmatrix.txt',sep = '\t',quote = F)


i<-0
aLN<-c()
for (i in 8:17) {
  group<-c(i*6+1,i*6+2,i*6+3)
  group<-as.character(lapply(group,function(x){paste0('X',x)}))
  testmatrix<-exp[,group]
  rownames(testmatrix)<-exp[,1]
  u<-testmatrix
  library(edgeR)
  keep <- rowSums(cpm(u) > 0.5) >= 2
  u<-u[keep,]
  library(psych)
  png(filename = paste0('~/Downloads/dd/',i*6+1,'-',i*6+3,'cor.png'))
  pairs.panels(u,
               gap = 0,
               bg = c("grey", "blue","black"),
               pch=21,
               digits=3,
               cex.cor = 1)
  dev.off()
  aLN<-append(aLN,group)
}
aLNm<-exp[,aLN]
rownames(aLNm)<-exp[,1]
u<-aLNm
library(edgeR)
keep <- rowSums(cpm(aLNm) > 0.5) >= 2
aLNm<-aLNm[keep,]
#pca
dat_scale=scale(aLNm,scale=F)#################
dat_cor=cor(dat_scale)
#options(digits=4, scipen=4)
#kable(dat_cor)
dat_eigen=eigen(dat_cor)
dat_var=dat_eigen$values ## 相关系数矩阵的特征值
dat_var
pca_var=dat_var/sum(dat_var)
pca_var
pca_cvar=cumsum(dat_var)/sum(dat_var)
pca_cvar
library(ggplot2)
p=ggplot(,aes(x=1:12,y=pca_var))
p1=ggplot(,aes(x=1:12,y=pca_cvar))
#p+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
#p1+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
pca_vect= dat_eigen$vector  ## 相关系数矩阵的特征向量
loadings=data.frame(sweep(pca_vect,2,sqrt(pca_var),"*"))
rownames(loadings)=colnames(aLNm)

qplot(loadings[,1], loadings[,2], data = loadings, 
      color=aLN,
      xlab = 'PC1 83.63%',ylab = 'PC2 5.83%')+
  geom_text(aes(loadings[,1], loadings[,2],label = aLN2))+
  theme(legend.title=element_blank(),legend.key.size = unit(10, "pt"),axis.title.x = element_text(size = 10))
dev.off()
aLN1<-c()
for (i in 1:10) {
  #  for (j in 1:3) {
  aLN1<-append(aLN1,lapply(i,function(x){rep(paste0('LN',x,'_R'),3)}))
}
aLN1<-as.character(unlist(aLN1))
aLN2<-c()
for (i in 1:30) {
  aLN2[i]<-paste0(aLN1[i],rep(1:3,8)[i])
}

aLB<-c()
for (i in 8:17) {
  group<-c(i*6+4,i*6+5,i*6+6)
  group<-as.character(lapply(group,function(x){paste0('X',x)}))
  testmatrix<-exp[,group]
  rownames(testmatrix)<-exp[,1]
  u<-testmatrix
  library(edgeR)
  keep <- rowSums(cpm(u) > 0.5) >= 2
  u<-u[keep,]
  library(psych)
  png(filename = paste0('~/Downloads/dd/',i*6+4,'-',i*6+6,'cor.png'))
  pairs.panels(u,
               gap = 0,
               bg = c("grey", "blue","black"),
               pch=21,
               digits=3,
               cex.cor = 1)
  dev.off()
  aLB<-append(aLB,group)
}
aLBm<-exp[,aLB]
rownames(aLBm)<-exp[,1]
u<-aLBm
library(edgeR)
keep <- rowSums(cpm(aLBm) > 0.5) >= 2
aLBm<-aLBm[keep,]
#pca
dat_scale=scale(aLBm,scale=F)
dat_cor=cor(dat_scale)
#options(digits=4, scipen=4)
#kable(dat_cor)
dat_eigen=eigen(dat_cor)
dat_var=dat_eigen$values ## 相关系数矩阵的特征值
dat_var
pca_var=dat_var/sum(dat_var)
pca_var
pca_cvar=cumsum(dat_var)/sum(dat_var)
pca_cvar
library(ggplot2)
p=ggplot(,aes(x=1:12,y=pca_var))
p1=ggplot(,aes(x=1:12,y=pca_cvar))
#p+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
#p1+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
pca_vect= dat_eigen$vector  ## 相关系数矩阵的特征向量
loadings=data.frame(sweep(pca_vect,2,sqrt(pca_var),"*"))
rownames(loadings)=colnames(aLBm)

qplot(loadings[,1], loadings[,2], data = loadings, 
      color=aLB,
      xlab = 'PC1 80.71%',ylab = 'PC2 7.49%')+
  geom_text(aes(loadings[,1], loadings[,2],label = aLB2))+
  theme(legend.title=element_blank(),legend.key.size = unit(10, "pt"),axis.title.x = element_text(size = 10))
dev.off()
aLB1<-c()
for (i in 1:10) {
  #  for (j in 1:3) {
  aLB1<-append(aLB1,lapply(i,function(x){rep(paste0('L',x,'B_R'),3)}))
}
aLB1<-as.character(unlist(aLB1))
aLB2<-c()
for (i in 1:30) {
  aLB2[i]<-paste0(aLB1[i],rep(1:3,8)[i])
}


colnames(aLBm)<-aLB2
colnames(aLNm)<-aLN2

write.table(aLBm,'~/Downloads/dd/iseqqc/LBmatrix.txt',sep = '\t',quote = F)
write.table(aLNm,'~/Downloads/dd/iseqqc/LNmatrix.txt',sep = '\t',quote = F)


dev.off()
n.sample=ncol(aHBm)
par(cex = 0.9)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(1,2))
setwd('~')
install.packages('rstatix')
boxplot(aHBm, col = cols,main="expression value",las=2)
boxplot(expr_new, col = cols,main="expression value",las=2)
?boxplot



exp<-read.table('~/Downloads/dd/ddmatrix.txt',header = T)
deglist<-read.csv('~/Downloads/dd/deglist.csv',header = F)
up<-list(c())
down<-list(c())
up2<-list(c())
down2<-list(c())
#i<-1
for (i in 1:length(deglist[,1])) {
  col<-as.character(unlist(deglist[i,1:6]))
testmatrix<-exp[,col]
rownames(testmatrix)<-exp[,1]
u<-testmatrix
library(edgeR)
keep <- rowSums(cpm(u) > 0.5) >= 2
u<-u[keep,]
type <- factor(c(rep(as.character(deglist[i,7]),3), rep(as.character(deglist[i,8]),3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u), type)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ type)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', as.character(deglist[i,7]), as.character(deglist[i,8])))
#7/8
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 2 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -2 & res1$padj < 0.05),'sig'] <- 'down'
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.csv(res1_select, file = paste0('~/Downloads/dd/degresults4fold/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = ',', col.names = NA, quote = FALSE)
up[[i]]<-rownames(res1)[which(res1$log2FoldChange >= 2 & res1$padj < 0.05)]
down[[i]]<-rownames(res1)[which(res1$log2FoldChange <= -2 & res1$padj < 0.05)]
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.05),'sig'] <- 'down'
#length(which(res1$log2FoldChange <= -1 & res1$padj < 0.05))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.csv(res1_select, file = paste0('~/Downloads/dd/degresults2fold/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = ',', col.names = NA, quote = FALSE)
up2[[i]]<-rownames(res1)[which(res1$log2FoldChange >= 1 & res1$padj < 0.05)]
down2[[i]]<-rownames(res1)[which(res1$log2FoldChange <= -1 & res1$padj < 0.05)]
}
save(up2,file = '~/Downloads/dd/uplog2fold1padj0.05.RData')
save(down2,file = '~/Downloads/dd/downlog2fold1padj0.05.RData')
unlist(up2[1:160])
union160<-unique(unlist(up2[1:160]))
ints160<-Reduce(intersect,up2[1:160])
union160f4<-unique(unlist(up[1:160]))
########set for 3 groups
Nxianup<-Reduce(intersect,up2[c(1,2,5,8,31,32,35,38)])
Nxiandown<-Reduce(intersect,down2[c(1,2,5,8,31,32,35,38)])
Bxianup<-Reduce(intersect,up2[c(81,82,85,88,111,112,115,118)])
Bxiandown<-Reduce(intersect,down2[c(81,82,85,88,111,112,115,118)])
Btmpup<-Reduce(intersect,up2[c(93,123,133,143,153)])
Btmpdown<-Reduce(intersect,down2[c(93,123,133,143,153)])
Ntmpup<-Reduce(intersect,up2[c(13,43,53,63,73)])
Ntmpdown<-Reduce(intersect,down2[c(13,43,53,63,73)])
Btmpup<-Reduce(intersect,up2[c(93,123,133,143,153)])
Btmpdown<-Reduce(intersect,down2[c(93,123,133,143,153)])
Ntropup<-Reduce(intersect,up2[c(30,24,26,27,29)])
Ntropdown<-Reduce(intersect,down2[c(30,24,26,27,29)])
Btropup<-Reduce(intersect,up2[c(110,104,106,107,109)])
Btropdown<-Reduce(intersect,down2[c(110,104,106,107,109)])
save(Nxianup,Nxiandown,Ntmpdown,Ntmpup,Ntropdown,Ntropup,Bxiandown,Bxianup,Btmpdown,Btmpup,Btropdown,Btropup,file = '~/Downloads/dd/3groupDEG.RData')

Reduce(union,list(up2[[1]]))
#H1N L1-10N
length(Reduce(union,list(up2[[1]],up2[[2]],up2[[3]],up2[[4]],up2[[5]],up2[[6]],up2[[7]],up2[[8]],up2[[9]],up2[[10]])))
listInput <- list(L1 = up2[[1]], L2 = up2[[2]], L3 = up2[[3]],L4 = up2[[4]],L5 = up2[[5]],L6 = up2[[6]],L7 = up2[[7]],L8 = up2[[8]],L9 = up2[[9]],L10 = up2[[10]])
#L1N H1-8N
length(Reduce(union,list(up2[[1]],up2[[11]],up2[[21]],up2[[31]],up2[[41]],up2[[51]],up2[[61]],up2[[71]])))
listInput <- list(H1 = up2[[1]], H2 = up2[[11]], H3 = up2[[21]],H4 = up2[[31]],H5 = up2[[41]],H6 = up2[[51]],H7 = up2[[61]],H8 = up2[[71]])
library(UpSetR)
upset(fromList(listInput), nsets = 8,order.by = "freq",
      mb.ratio=c(0.7, 0.3),
      #      set_size.scale_max = 500,
      #      set_size.show = TRUE,
      nintersects = 100,color.pal = 1,
      keep.order = TRUE,
      #      matrix.color ="#b35806", main.bar.color = colorPalette,
      #      sets.bar.color = c("#e41a1c","#377eb8","#4daf4a"), 
)



up2l<-c()
down2l<-c()
for (i in 1:178) {
  up2l<-append(up2l,length(up2[[i]]))
  down2l<-append(down2l,length(down2[[i]]))
}
write.csv(up2l,file = '~/Downloads/dd/up2foldlength.csv')
write.csv(down2l,file = '~/Downloads/dd/down2foldlength.csv')

testmatrix<-exp[,c('X52','X53','X54','X4','X5','X6')]
rownames(testmatrix)<-exp[,1]
u<-testmatrix
library(edgeR)
keep <- rowSums(cpm(u) > 0.5) >= 2
u<-u[keep,]
type <- factor(c(rep("LB1",3), rep("HB1",3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u), type)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ type)
plotDispEsts(dds1)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', 'LB1', 'HB1'))
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.1),'sig'] <- 'up'
length(which(res1$log2FoldChange >= 1 & res1$padj < 0.1))
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.1),'sig'] <- 'down'
length(which(res1$log2FoldChange <= -1 & res1$padj < 0.1))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = paste0('~/Downloads/test/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = '\t', col.names = NA, quote = FALSE)
up1<-rownames(res1)[which(res1$log2FoldChange >= 2 & res1$padj < 0.01)]
down1<-rownames(res1)[which(res1$log2FoldChange <= -2 & res1$padj < 0.01)]

testmatrix<-exp[,c('X49','X50','X51','X1','X2','X3')]
rownames(testmatrix)<-exp[,1]
u<-testmatrix
keep <- rowSums(cpm(u) > 0.5) >= 2
u<-u[keep,]
type <- factor(c(rep("LN1",3), rep("HN1",3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u), type)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ type)
plotDispEsts(dds1)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', 'LN1', 'HN1'))
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 2 & res1$padj < 0.01),'sig'] <- 'up'
length(which(res1$log2FoldChange >= 2 & res1$padj < 0.01))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = paste0('~/Downloads/test/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = '\t', col.names = NA, quote = FALSE)
up2<-rownames(res1)[which(res1$log2FoldChange >= 2 & res1$padj < 0.01)]
down2<-rownames(res1)[which(res1$log2FoldChange <= -2 & res1$padj < 0.01)]

length(intersect(up1,up2))
length(intersect(down1,down2))

testmatrix<-exp[,c('X10','X11','X12','X7','X8','X9')]
testmatrix<-exp[,c('X4','X5','X6','X1','X2','X3')]
rownames(testmatrix)<-exp[,1]
u<-testmatrix
library(edgeR)
keep <- rowSums(cpm(u) > 0.5) >= 2
u<-u[keep,]
type <- factor(c(rep("HB1",3), rep("HN1",3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u),type)
meta<-cbind(meta,batch)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ batch+type)
plotDispEsts(dds1)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', 'HB1', 'HN1'))
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.05),'sig'] <- 'up'
length(which(res1$log2FoldChange >= 1 & res1$padj < 0.05))
length(which(res1$log2FoldChange <= -1 & res1$padj < 0.05))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = paste0('~/Downloads/test/2',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = '\t', col.names = NA, quote = FALSE)

testmatrix<-exp[,c('X52','X53','X54','X49','X50','X51')]
rownames(testmatrix)<-exp[,1]
u<-testmatrix
keep <- rowSums(cpm(u1) > 0.5) >= 2
u1<-u1[keep,]
type <- factor(c(rep("LB1",3), rep("LN1",3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u), type)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ type)
plotDispEsts(dds1)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', 'LB1', 'LN1'))
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 2 & res1$padj < 0.01),'sig'] <- 'up'
length(which(res1$log2FoldChange >= 2 & res1$padj < 0.01))
length(which(res1$log2FoldChange <= -2 & res1$padj < 0.01))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = paste0('~/Downloads/test/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = '\t', col.names = NA, quote = FALSE)
###################################################3
rld <- rlog(dds, blind = FALSE)
dds <- estimateSizeFactors(dds)
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                          colData=colData(dds))
dds <- assay(dds)
dds<-cor(dds)
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
#pca
colData(rld)
plotPCA(rld,intgroup=c('batch','type'))
#########correct batch effect
assay(rld) <- limma::removeBatchEffect(assay(rld), rld$batch)
pcaData <- plotPCA(rld, intgroup=c("batch", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=type, shape = batch)) +
  geom_point(size=3) +
#  xlim(-12, 12) +
#  ylim(-10, 10) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text(aes(label=name),vjust=2)
ggsave("myPCABatchEffectRemoved.png")

pcaData <- plotPCA(rld, intgroup = c( "type"), returnData = TRUE)
pcaData <- plotPCA(DESeqTransform( se ), intgroup = c( "type"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar
library(ggplot2)
pdf(file = '~/Downloads/PCAr_by_ggplot.pdf',width = 8,height = 4)
ggplot(pcaData, aes(x = PC1, y = PC2, color = type, shape = name)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with rld data")
dev.off()
#distribution
  pdf("Counts_distribution.pdf", height = 5, width = 6)
  boxplot(rld@assays@data@listData[[1]], pch=".", horizontal=TRUE, cex.axis=0.9, main= "Distribution of counts per sample", las=1, xlab="normalized counts by rlog", col="red")
  dev.off()
  colnames(u1)[49:54]<-c('L9B_R1','L9B_R2','L9B_R3','L10B_R1','L10B_R2','L10B_R3')
  pdf("Downloads/LHCounts_distribution.pdf", height = 10, width = 12)
  boxplot(log2(u1 + 1), pch=".", horizontal=TRUE, cex.axis=0.9, main= "Distribution of counts per sample", las=1, xlab="log2(counts+1)", col="red")
  dev.off()
# Plot heatmap
  library(pheatmap)
  pheatmap(rld_cor, annotation = meta[, c("type"), drop=F])
#corelation
   u<-rld_mat
   u<-dds
   u<-exp
   library(psych)
pairs.panels(u,
               gap = 0,
               bg = c("grey", "blue","black"),
               pch=21,
               digits=3,
               cex.cor = 1)
u<-apply(exp[,2:109], 2, as.numeric)
u<-cbind(exp[,1],u)
u[,1]<-exp[,1]
u<-u[,-1]


#pca
dat_scale=scale(lh18,scale=F)
dat_cor=cor(dat_scale)
#options(digits=4, scipen=4)
#kable(dat_cor)
dat_eigen=eigen(dat_cor)
dat_var=dat_eigen$values ## 相关系数矩阵的特征值
dat_var
pca_var=dat_var/sum(dat_var)
pca_var
pca_cvar=cumsum(dat_var)/sum(dat_var)
pca_cvar
library(ggplot2)
p=ggplot(,aes(x=1:12,y=pca_var))
p1=ggplot(,aes(x=1:12,y=pca_cvar))
#p+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
#p1+geom_point(pch=2,lwd=3,col=2)+geom_line(col=2,lwd=1.2)
pca_vect= dat_eigen$vector  ## 相关系数矩阵的特征向量
loadings=data.frame(sweep(pca_vect,2,sqrt(pca_var),"*"))
rownames(loadings)=colnames(lh18)

ghl<-read.table('~/Downloads/dd/iseqqc/LBHBphe1.txt',header = T)
ghl<-read.csv('~/Downloads/dd/iseqqc/group18.csv')
pdf('~/Downloads/lh18pca.pdf')
qplot(loadings[,1], loadings[,2], data = loadings, 
      color=ghl$groups,
      xlab = 'PC1 74.50%',ylab = 'PC2 14.17%')+
  geom_text(aes(loadings[,1], loadings[,2],label = colnames(lh18)))+
  theme(legend.title=element_blank(),legend.key.size = unit(10, "pt"),axis.title.x = element_text(size = 10))
dev.off()

write.csv(u1,'~/Downloads/lh.csv')

lh18<-data.frame()
for (i in 1:length(u1[,1])) {
  lh18[i,2]<-mean(u1[i,4:6])
  lh18[i,3]<-mean(u1[i,7:9])
  lh18[i,4]<-mean(u1[i,10:12])
  lh18[i,5]<-mean(u1[i,13:15])
  lh18[i,6]<-mean(u1[i,16:18])
  lh18[i,7]<-mean(u1[i,19:21])
  lh18[i,8]<-mean(u1[i,22:24])
  lh18[i,9]<-mean(u1[i,25:27])
  lh18[i,10]<-mean(u1[i,28:30])
  lh18[i,11]<-mean(u1[i,31:33])
  lh18[i,12]<-mean(u1[i,34:36])
  lh18[i,13]<-mean(u1[i,37:39])
  lh18[i,14]<-mean(u1[i,40:42])
  lh18[i,15]<-mean(u1[i,43:45])
  lh18[i,16]<-mean(u1[i,46:48])
  lh18[i,17]<-mean(u1[i,49:51])
  lh18[i,18]<-mean(u1[i,52:54])
}
rownames(lh18)<-rownames(u1)
colnames(lh18)<-c('H1B','H2B','H3B','H4B','H5B','H6B','H7B','H8B','L1B','L2B','L3B','L4B','L5B','L6B','L7B','L8B','L9B','L10B')

VennDiagram::venn.diagram(
  x= list("Gj-tmp" = Btmpup, "Gj-trop" =Btropup,"xian"=Bxianup),
  filename = "~/Downloads/up.png", 
  height = 800, 
  width = 800,
  resolution =300, 
  imagetype="png", 
  col="transparent",
  fill=c("blue","green",'yellow'),
  alpha = 0.50, 
  cex=1,
  cat.cex=0.5
)
VennDiagram::venn.diagram(
  x= list("Gj-tmp" = Btmpdown, "Gj-trop" =Btropdown,"xian"=Bxiandown),
  filename = "~/Downloads/down.png", 
  height = 800, 
  width = 800,
  resolution =300, 
  imagetype="png", 
  col="transparent",
  fill=c("blue","green",'yellow'),
  alpha = 0.50, 
  cex=1,
  cat.cex=0.5
)
VennDiagram::venn.diagram(
  x= list("Gj-tmp up" = Btmpup, "Gj-trop up" =Btropup,"xian up"=Bxianup,"Gj-tmp down" = Btmpdown, "Gj-trop down" =Btropdown),
  filename = "~/Downloads/all.png", 
  height = 800, 
  width = 800,
  resolution =300, 
  imagetype="png", 
  col="transparent",
  fill=c("blue","green",'yellow','red','brown'),
  alpha = 0.50, 
  cex=1,
  cat.cex=0.5
)
############go kegg
load('~/Downloads/dd/3groupDEG.RData')
LOC <- read.csv("~/Downloads/研究生其他/mengyun/mengyun-小明/annotnow.csv")
library(clusterProfiler)
# install.packages("https://github.com/xuzhougeng/org.Osativa.eg.db/releases/download/v0.01/org.Osativa.eg.db.tar.gz", 
#                  repos = NULL, 
#                  type="source")
library(org.Osativa.eg.db)
library(ggplot2)
library(stringr)
org <- org.Osativa.eg.db
# msu_id: MSU ID
rap_id<-Btmpdown
now<-'Btmpdown'
rap_id<-Btmpup
now<-'Btmpup'
rap_id<-Btropdown
now<-'Btropdown'
rap_id<-Btropup
now<-'Btropup'

rap_id<-Bxiandown
now<-'Bxiandown'
rap_id<-Bxianup
now<-'Bxianup'

map_id <- AnnotationDbi::select(org, keys = rap_id, 
                                columns=c("RAP"), keytype = "GID") 
rapid<-map_id$RAP[-which(is.na(map_id$RAP))]
#rap<-read.table('~/Downloads/btmpd.txt_out')
ego <- enrichGO(rapid,
                OrgDb = org,
                keyType = "RAP",
                ont="ALL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                qvalueCutoff =0.1
                )
go_df <- as.data.frame(ego)
write.table(go_df, file=paste0("~/Downloads/dd/output/go_df",now,".txt"),  sep="\t", row.names=FALSE, quote=FALSE)
png(filename =paste0("~/Downloads/dd/output/go_bar",now,".png"),width = 600,height = 600,units = "px",res = 80)
barplot(ego,drop = TRUE, showCategory =10,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
# ego <- enrichGO(rap$V2,
#                 OrgDb = org,
#                 keyType = "RAP",
#                 ont="ALL",
#                 pvalueCutoff = 1,
#                 pAdjustMethod = "BH",
#                 qvalueCutoff = 1
# )
#plotGOgraph(ego)
#heatplot(ego)
#goplot(ego)       #igraph布局的DAG
#library(enrichplot)
#library(ggnewscale)
png(filename =paste0("~/Downloads/dd/output/go_emap",now,".png"),,width = 1000,height = 1000)
ego2 <- pairwise_termsim(ego)
emapplot(ego2)
dev.off()
#emapplot(ego, showCategory = 30)#GO terms关系网络图（通过差异基因关联）
png(filename =paste0("~/Downloads/dd/output/go_cnet",now,".png"),width = 1000,height = 1000)
#ego2 <- pairwise_termsim(ego)
cnetplot(ego, showCategory = 5)#GO term与差异基因关系网络图
dev.off()
png(filename =paste0("~/Downloads/dd/output/go_bar",now,".png"),width = 600,height = 600,units = "px",res = 80)
barplot(ego,drop = TRUE, showCategory =10,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
#library(DO.db)

load('~/Downloads/研究生其他/mengyun/mengyun-小明/AnnotationMerge.RData')
keggid<-c()
for (i in rapid) {
  keggid<-append(keggid,AnnotationMerge[which(AnnotationMerge[,1]==i),2])
}
# dosaall<-read.table('~/Downloads/dosa.txt')
# keggid1<-intersect(dosaall[,2],keggid)

engid<-c()
for (i in rapid) {
  engid<-append(engid,AnnotationMerge[which(AnnotationMerge[,1]==i),7])
}
# osaall<-read.table('~/Downloads/osaall.txt')
# engid1<-intersect(as.character(osaall[,2]),engid)

kk <- enrichKEGG(
  organism = "dosa",
  gene = keggid,
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff =1
)
keg<-read.csv("~/Downloads/dosa00001.keg",sep = ';')
geneList = kk@result$pvalue
## feature 2: named vector
names(geneList) = kk@result$geneID
## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)
kk2 <- gseKEGG(
  geneList  = geneList,
  keyType  = 'kegg',
  organism = 'dosa',
  nPerm  = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod     = "BH"
)
kegg_df <- as.data.frame(kk)
write.table(kegg_df, file=paste0("~/Downloads/dd/output/kkdosa_df",now,".txt"),  sep="\t", row.names=FALSE, quote=FALSE)

png(filename =paste0("~/Downloads/dd/output/kkdosa_bar",now,".png"),width = 600,height = 600,units = "px",res = 80)
barplot(kk,drop = TRUE, showCategory =20)
dev.off()

kk <- enrichKEGG(
  gene = engid,
  organism = "osa",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff =1
)
kegg_df <- as.data.frame(kk)
write.table(kegg_df, file=paste0("~/Downloads/dd/output/kkosa_df",now,".txt"),  sep="\t", row.names=FALSE, quote=FALSE)

png(filename =paste0("~/Downloads/dd/output/kkosa_bar",now,".png"),width = 600,height = 600,units = "px",res = 80)
barplot(kk,drop = TRUE, showCategory =20)
dev.off()
#fold change配套表格？


require(AnnotationHub)
hub <- AnnotationHub()
query(hub, "oryza sativa")
rice <- hub[['AH94060']]
columns(rice)
keys(rice)
rice$export
require(clusterProfiler)
bitr(map_id$RAP, 'GID', c("REFSEQ", "GO", "ONTOLOGY"), rice)
ego <- enrichGO(map_id$RAP,
                OrgDb = rice,
                keyType = "GID",
                ont="ALL",
                pvalueCutoff = 1,
                pAdjustMethod = "bonferroni",
                qvalueCutoff = 1
)
alldeg<-Reduce(union,list(Btmpdown,Btmpup,Bxiandown,Bxianup,Btropdown,Btropup))
write.table(alldeg,'~/Downloads/alldeg.txt',col.names = F,row.names = F,quote = F)
alldeg<-as.character(read.table('~/Downloads/alldeg.txt')[,1])
#wcgna
exp<-read.table('~/Downloads/dd/ddmatrix.txt',header = T)
name18<-read.csv('~/Downloads/18x3name.csv',header = T)
dataExpr<-cbind(mxian,mtmp,mtrop)
colnames(dataExpr)<-name18[,1]
rownames(dataExpr)<-exp$Geneid
save(dataExpr,file='~/Downloads/54sampleforwcgna.RData')
load('~/Downloads/54sampleforwcgna.RData')
library(DESeq2)
u<-dataExpr[alldeg,]
u1<-u
u<-u1
u<-dataExpr
# library(edgeR)
# keep <- rowSums(cpm(u) > 0.5) >= 2
# u<-u[keep,]
#type <- factor(c(rep("xianL",6), rep("xianH",12),rep("tmpH",15),rep("tmpL",3),rep("tropH",3),rep("tropL",15)))
#u<-u1[,1:6]
#type <- factor(c(rep("xianH",2), rep("xianL",4)))
type <- factor(c(rep("xianH",2), rep("xianL",4),rep("tmpH",5),rep("tmpL",1),rep("tropH",1),rep("tropL",5)))
type <- factor(c(rep("HighDormancy",2), rep("LowDormancy",4),rep("HighDormancy",5),rep("LowDormancy",1),rep("HighDormancy",1),rep("LowDormancy",5)))
meta <- data.frame(row.names=colnames(u), type)
#u<-apply(u, 2, round)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~type)
dds<- estimateSizeFactors(dds)  
dataExpr<- counts(dds, normalized=TRUE)
 
load('~/Downloads/dd/nor18.RData')
sizeFactors(dds.sizefactor )
ddss<-est
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- DESeq2::results(dds1, contrast = c('type', 'xianL', 'xianH'))
res1<-res@listData

dds2<-dds1@assays@data@listData[[1]]
dataExpr<-vst(dds,blind = T)
dataExpr<-dataExpr@assays@data@listData[[1]]
which(rownames(dataExpr)==Bxianup[1])
#dataExpr<-dataExpr[[1]]
save(dataExpr,file='~/Downloads/nor18.RData')
#
load('~/Downloads/vst54.RData')
dataExpr<-t(dataExpr)
## 检测缺失值
library(WGCNA)
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
#dataExpr<-dataExpr[-which(row.names(dataExpr)=='shoot282set_CML311'),]
powers = c(c(1:10), seq(from = 12, to=30, by=2))
type = "unsigned"
#signed方法来削减负相关的影响
#因为一群基因正负调控得到的结果，这时候unsigned的比较合适。
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5,dataIsExpr = T,
                        corFnc = 'bicor')
png(file='~/Downloads/dd/WGCNA/sft.png',width = 1000,height = 1000)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
corType = "bicor"
corFnc = 'bicor'
#对二元变量,如样本性状信息计算相关性时,# 或基因表达严重依赖于疾病状态时
#maxPOutliers = ifelse(corType=="pearson",1,0)
maxPOutliers<-0
exprMat <- "~/Downloads/dd/WGCNA/dd"
cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,deepSplit = 2,
                       #                       detectCutHeight = 0.995,
                       TOMType = 'unsigned', minModuleSize = 30,
                       #                       maxCoreScatter = 1,maxAbsCoreScatter=1,
                       pamStage=T,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.1,
                       #                       trapErrors = T,
                       #                      stabilityCriterion='Common fraction',
                       numericLabels = F, 
                       #                      pamRespectsDendro = T,
                       saveTOMs=TRUE, 
                       corType = 'bicor', 
                       #corFnc = cor,corOptions = list(use='p',method='spearman')
                       #corFnc = cor,corOptions = list(use='p',method='kendall'),
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3
                       #                       useBranchEigennodeDissim = T,minBranchEigennodeDissim = 0.4,
                       #                       reassignThreshold= 1e-6,
                       #                       minCoreKME = 0.6
)
table(net$colors)
write.table(table(net$colors),file = paste0('~/Downloads/WGCNA/',i,'table.txt'),quote = F)

## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting

#moduleLabels = net$colors
moduleColors =net$colors
#moduleColors = labels2colors(moduleLabels)

# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
png('~/Downloads/dd/WGCNA/tree.png')
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
#MEs_col = MEs
#library(tidyverse)
#colnames(MEs_col) = paste0("ME", labels2colors(
#  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
#MEs_col = orderMEs(MEs_col)
write.csv(MEs,file = paste0('~/Downloads/WGCNA/',i,'me.csv'))
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
png('~/Downloads/dd/WGCNA/eigengene_adjacency_heatmap_.png')
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
name18<-read.csv('~/Downloads/dd/18name.csv',header = T,row.names = 1)
traitData<-read.csv('~/Downloads/dd/trait18.csv',header = T,row.names = 1)
rownames(name18)<-name18[,1]
name18<-name18[,-1]
traitData<-name18
moduleTraitCor = cor(MEs, traitData, use = "p")
###########
### 模块与表型数据关联
if (corType=="pearson") {
  modTraitCor = cor(MEs, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(oma=c(1,1,1,1),mar=c(1,1,1,1))
par(mar=c(6,15,3,3))
pdf('~/Downloads/novstmoduletrait.pdf',paper = 'a4r')
pdf('~/Downloads/dd/WGCNA/moduletrait.pdf',paper = 'a4r')
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
               yLabels = colnames(MEs),
               cex.lab = 0.8,
               cex.lab.y=0.5,
               ySymbols = colnames(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
#输出每个模块基因
for (mod in 1:nrow(table(moduleColors))){
  modules = names(table(moduleColors))[mod] 
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules) 
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0('~/Downloads/dd/WGCNA/vst',modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)}
#表达图谱
"magenta"
"black"
"greenyellow"
"pink"
"salmon"
"turquoise"
"green"

"yellow"
'red'
"turquoise"
"green"
"magenta"
"pink"

"black"
"green"
'red'
"yellow"
"pink"

pdf(paste0('~/Downloads/dd/WGCNA/vst',which.module,'expressionplot.pdf'),height = 500,width = 800)
sizeGrWindow(8,7);
which.module="turquoise"
which.module="blue"
which.module="brown"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(dataExpr[,moduleColors==which.module])),
        nrgcols=30,rlabels=F,clabels=rownames(dataExpr),rcols=which.module,
        main="", cex.main=1)
#        ,cex.axes=0.05,cex.lab=0.05)
dev.off()

datKME=signedKME(dataExpr, MEs, outputColumnName="MM.")
write.csv(datKME,'~/Downloads/dd/12KME.csv')

hierTOM = hclust(as.dist(dissTOM),method="average");
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                    deepSplit=2, pamRespectsDendro = FALSE))
ADJ1=abs(cor(dataExpr,use="p"))^6
dissTOM=TOMdist(ADJ1)
colorh1<-colorDynamicHybridTOM 
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)
colorlevels=unique(colorh1)
GS1=as.numeric(cor(,dataExpr, use="p"))
class(dataExpr[1,1])
GeneSignificance=abs(GS1)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (colorh1==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=colorh1[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
table(FilterGenes)
dimnames(data.frame(datExpr))[[2]][FilterGenes]

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(dataExpr, power = sft$powerEstimate); 
# Select module
module = "red";
# Select module probes
probes = colnames(dataExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("~/Downloads/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("~/Downloads/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);

###############################################################3








data(geneList, package="DOSE")
head(geneList)
geneList
gene <- names(geneList)[abs(geneList) > 2]

exp<-read.table('~/Downloads/dd/ddmatrix.txt',header = T)
exp1<-apply(exp, 2, as.numeric)
#前6高 
mxian<-exp1[,c('X4','X5','X6','X22','X23','X24','X52','X53','X54','X58','X59','X60','X76','X77','X78','X94','X95','X96')]
mxian1<-data.frame()
for (i in 1:length(mxian[,1])) {
#  i<-1
  mxian1[i,1]<-mean(mxian[i,1:6])
}
for (i in 1:length(mxian[,2])) {
  #  i<-1
  mxian1[i,2]<-mean(mxian[i,7:18])
}
#前15高
mtmp<-exp1[,c('X10','X11','X12','X28','X29','X30','X34','X35','X36','X40','X41','X42','X46','X47','X48','X64','X65','X66')]
mtmp1<-data.frame()
for (i in 1:length(mtmp[,1])) {
  #  i<-1
  mtmp1[i,1]<-mean(mtmp[i,1:15])
}
for (i in 1:length(mtmp[,2])) {
  #  i<-1
  mtmp1[i,2]<-mean(mtmp[i,16:18])
}
#前3高
mtrop<-exp1[,c('X16','X17','X18','X70','X71','X72','X82','X83','X84','X88','X89','X90','X100','X101','X102','X106','X107','X108')]
mtrop1<-data.frame()
for (i in 1:length(mtrop[,1])) {
  #  i<-1
  mtrop1[i,1]<-mean(mtrop[i,1:3])
}
for (i in 1:length(mtrop[,2])) {
  #  i<-1
  mtrop1[i,2]<-mean(mtrop[i,4:18])
}
m3<-cbind(mxian,mtmp,mtrop)
rownames(m3)<-exp$Geneid
g3<-Reduce(union,list(Bxiandown,Bxianup,Btmpdown,Btmpup,Btropdown,Btropup))
m3<-m3[g3,]
makemeanmatrix<-function(input=m3){
m3o<-data.frame()
for (i in 1:length(input[,1])) {
 for (j in 1:18) {
   m3o[i,j]<-mean(input[i,(j*3-2):(j*3)])
 }
}
return(m3o)
}
m3o<-makemeanmatrix(input = m3)
rownames(m3o)<-rownames(m3)
name18<-read.csv('~/Downloads/dd/18name.csv')
colnames(m3o)<-name18$id
write.csv(m3o,file = '~/Downloads/dd/s18g4703.csv')









rownames(testmatrix)<-exp[,1]
u<-m3o
library(edgeR)
keep <- rowSums(cpm(u) > 0.5) >= 2
u<-u[keep,]
type <- factor(c(rep("LN1",3), rep("HN1",3)))
#database <- round(as.matrix(u))
meta <- data.frame(row.names=colnames(u), type)
dds <- DESeqDataSetFromMatrix(u, colData = meta, design = ~ type)
plotDispEsts(dds1)
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)
res <- results(dds1, contrast = c('type', 'LN1', 'HN1'))
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res1[which(res1$log2FoldChange >= 2 & res1$padj < 0.01),'sig'] <- 'up'
length(which(res1$log2FoldChange >= 2 & res1$padj < 0.01))
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = paste0('~/Downloads/test/',as.character(type[1]),as.character(type[4]),'deg.csv'), sep = '\t', col.names = NA, quote = FALSE)
up2<-rownames(res1)[which(res1$log2FoldChange >= 2 & res1$padj < 0.01)]
down2<-rownames(res1)[which(res1$log2FoldChange <= -2 & res1$padj < 0.01)]


write.table(Btmpdown,'~/Downloads/btmpd.txt',col.names = F,row.names = F,quote = F)
egoALL.5055 <- enrichGO(ALL.5055.vs.mt$EntrezgeneID,
                        OrgDb = org,
                        keyType = "GID",
                        pAdjustMethod = "none",
                        ont="BP")

