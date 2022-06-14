load('~/Downloads/7tiphe.RData')
write.csv(pheroot,file = '~/Downloads/pheroot.csv')
write.csv(pheshoot,file = '~/Downloads/pheshoot.csv')
write.csv(phel3b,file = '~/Downloads/phel3b.csv')
write.csv(phel3t,file = '~/Downloads/phel3t.csv')
write.csv(phekern,file = '~/Downloads/phekern.csv')
write.csv(pheld,file = '~/Downloads/pheld.csv')
write.csv(pheln,file = '~/Downloads/pheln.csv')

samp2000<-sample(row.names(exp_72),2000,replace = F)
ti<-c('pheroot','pheshoot','phel3b','phel3t','phekern','pheld','pheln')
#ti<-c('phel3b','phel3t','phekern','pheld','pheln')
i<-'pheroot'
i<-'pheshoot'
i<-'phel3t'
i<-'pheld'
for (i in ti) {
datae<-read.csv(paste0('~/Downloads/WGCNA/data/',i,'.csv'))
datae<-datae[,-1]
colnames(datae)<-unlist(datae[1,])
row.names(datae)<-unlist(datae[,1])
dataExpr<-datae[-1,-1]
datae<-apply(dataExpr, 2, as.numeric)
row.names(datae)<-row.names(dataExpr)
dataExpr<-datae
dataExpr <- as.data.frame(t(dataExpr))
#col2000<-sample(length(dataExpr[1,]),2000)
dataExpr1<-dataExpr[,samp2000]
dataExpr<-dataExpr1
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
png(paste0('~/Downloads/WGCNA/',i,'sft.png'))
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
exprMat <- paste0("~/Downloads/WGCNA/",i)
cor <- WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,deepSplit = 2,
#                       detectCutHeight = 0.995,
                       TOMType = 'unsigned', minModuleSize = 30,
#                       maxCoreScatter = 1,maxAbsCoreScatter=1,
                       pamStage=F,
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
png(paste0('~/Downloads/WGCNA/tree',i,'.png'))
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
png(paste0('~/Downloads/WGCNA/eigengene_adjacency_heatmap_',i,'.png'))
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
#计算第二主成分特征值#
MEs0 = moduleEigengenes(expr = dataExpr, colors = moduleColors,nPC = 2,verbose = 6)
MEs2
#输出每个模块的基因
for (mod in 1:nrow(table(moduleColors))){
  modules = names(table(moduleColors))[mod] 
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules) 
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0('~/Downloads/WGCNA/',i,modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)}
}
mert<-MEs_col
#eigengwas#练习
phroot<-MEs
load('~/Downloads/LEI/data/maf5pergroot50.RData')
all<-as.character(phdata(data_rna)[,1])
now<-rownames(MEs)
nowc<-c()
for (i in 1:length(now)) {
  nowc[i]<-substr(now[i],6,nchar(now[i]))
}
needm<-matrix(nrow = length(need),ncol = length(MEs[1,]))
phdata<-rbind(MEs,needm)
needm<-matrix(nrow = 298,ncol = 2)
colnames(needm)<-colnames(MEs)
need<-setdiff(all,nowc)
needm[,1]<-c(nowc,need)
needm[,2]<-1
ph<-cbind(needm,phdata)
colnames(ph)[1]<-'id'
colnames(ph)[2]<-'sex'
ph<-ph[match(phdata(data_rna)$id,ph$id),]
phdata(data_rna)<-ph
library(GenABEL)
lambdas<-c()
data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
for (i in 4:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/sh",colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/sh",colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/shlam",i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}

#组织间模块对比
allmodule<-as.character(read.table('~/Downloads/WGCNA/allmodule.txt')[,1])
allmodule<-as.character(read.table('~/Downloads/WGCNA/allmodule1.txt')[,1])
i<-allmodule[1]
for (i in 1:length(allmodule)) {
  for(j in i:length(allmodule)){
    a<-as.character(read.csv(paste0('~/Downloads/WGCNA/',allmodule[i]))[,1])
    b<-as.character(read.csv(paste0('~/Downloads/WGCNA/',allmodule[j]))[,1])
    if(length(intersect(a,b))/length(a)>0.7&i!=j) print(paste0(allmodule[i],'&',allmodule[j],' are similar'))
    if(length(intersect(a,b))/length(b)>0.7&i!=j) print(paste0(allmodule[i],'&',allmodule[j],' are similar'))
    }
}
olm<-as.data.frame(matrix(nrow = length(allmodule),ncol = length(allmodule)))
allmodule1<-as.character(read.table('~/Downloads/WGCNA/allmodule1.txt')[,1])
allmodule2<-as.character(read.table('~/Downloads/WGCNA/allmodule2.txt')[,1])

colnames(olm)<-allmodule2
rownames(olm)<-allmodule2
outma<-matrix(nrow = 1,ncol = 7)
for (i in 1:length(allmodule)) {
  for(j in (i+1):length(allmodule)){
    a<-as.character(read.csv(paste0('~/Downloads/WGCNA/',allmodule[i]),header = F)[,1])
    b<-as.character(read.csv(paste0('~/Downloads/WGCNA/',allmodule[j]),header = F)[,1])
#    olm[i,j]<-paste0('(',length(a),',',length(b),')',length(intersect(a,b)))
    if(length(intersect(a,b))/length(a)>0.7|length(intersect(a,b))/length(b)>0.7) {print(paste0(allmodule[i],'(',length(a),')&',allmodule[j],'(',length(b),')',length(intersect(a,b)),' are similar'))
#    olm[i,j]<-paste0('(',length(a),',',length(b),')',length(intersect(a,b)))}
#    olm[i,j]<-paste0(round(length(intersect(a,b))/length(a)*100,2),'%,',round(length(intersect(a,b))/length(b)*100,2),'%')}
    outl<-matrix(nrow = 1,ncol = 7)
    outl[1,1]<-allmodule1[i]
    outl[1,2]<-length(a)
    outl[1,3]<-allmodule1[j]
    outl[1,4]<-length(b)
    outl[1,5]<-length(intersect(a,b))
    outl[1,6]<-round(length(intersect(a,b))/length(a)*100,2)
    outl[1,7]<-round(length(intersect(a,b))/length(b)*100,2)
    outma<-rbind(outma,outl)
    #    if(length(intersect(a,b))/length(b)>0.7&i!=j) print(paste0(allmodule[i],'&',allmodule[j],' are similar'))
    }
    }
}
outma<-outma[-1,]
write.csv(outma,'~/Downloads/WGCNA/23relation_of_modules_7col.csv')

write.csv(olm,'~/Downloads/WGCNA/relation_of_modules_percentage.csv')
write.csv(olm,'~/Downloads/WGCNA/relation_of_modules_peronly.csv')
#eigengwas
i<-7
MEs<-read.csv(paste0('~/Downloads/WGCNA/',ti[i],'me.csv'))
#phroot<-MEs
#load('~/Downloads/LEI/data/maf5pergroot50.RData')
#all<-as.character(phdata(data_rna)[,1])
now<-as.character(MEs[,1])
#now<-rownames(MEs)
nowc<-c()
for (i in 1:length(now)) {
  nowc[i]<-substr(now[i],3,nchar(now[i]))
}
needm<-matrix(nrow = (298-(length(MEs[,1]))),ncol = length(MEs[1,]))
colnames(needm)<-colnames(MEs)
phdata<-rbind(MEs,needm)
needm<-matrix(nrow = 298,ncol = 2)
#colnames(needm)<-colnames(MEs)
need<-setdiff(all,nowc)
needm[,1]<-c(nowc,need)
needm[,2]<-1
ph<-cbind(needm,phdata)
colnames(ph)[1]<-'id'
colnames(ph)[2]<-'sex'
ph<-ph[match(phdata(data_rna)$id,ph$id),]
ph<-ph[,-3]
phsh<-ph
phl3b<-ph
phl3t<-ph
phke<-ph
phld<-ph
phln<-ph
phrt<-ph
phdata(data_rna)<-ph
library(GenABEL)
lambdas<-c()



data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
j<-1
lambdas<-c()
phdata(data_rna)<-phrt
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}

j<-2
lambdas<-c()
phdata(data_rna)<-phsh
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
j<-3
lambdas<-c()
phdata(data_rna)<-phl3b
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
j<-4
lambdas<-c()
phdata(data_rna)<-phl3t
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
j<-5
lambdas<-c()
phdata(data_rna)<-phke
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
j<-6
lambdas<-c()
phdata(data_rna)<-phld
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
j<-7
lambdas<-c()
phdata(data_rna)<-phln
for (i in 3:length(phdata(data_rna)[1,])) {
  data1<-phdata(data_rna)[,i]
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm <- mmscore(hy, data=data_rna)
  lambdas[i]<-lambda(mm)$estimate
  png(filename=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna))[i],".PNG"))
  par(mfrow=c(3,1))
  plot(mm)
  thres<-2.04E-8
  abline(h = -log10(thres),lty="dashed",col="red")
  estlambda(mm[, "P1df"], plot=TRUE)
  hist(phdata(data_rna)[1:length(phdata(data_rna)[,1]),i],main=paste("Histogram of",colnames(phdata(data_rna)[i]),xlab = "expression"))
  dev.off()
  save(mm,file=paste0("~/Downloads/WGCNA/",ti[j],colnames(phdata(data_rna)[i]),".RData"))
  write.table(lambdas[i],file=paste0("~/Downloads/WGCNA/lam",ti[j],i,colnames(phdata(data_rna)[i])),row.names = F,col.names = F)
}
