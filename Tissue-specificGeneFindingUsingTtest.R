#this script including
#t statstics to find tissue-specific gene
#go &kegg for genes
#t-test for tissue-specific gene and remaining genes in each tissue using dnds

load('~/Downloads/maize_GTEx/data/15505_7ti.RData')
Gene_Expression_T<-exp_71
###center and scale all the columens, i.e., normalize across smaples
Gene_Expression_T <- t(Gene_Expression)

head(Gene_Expression_T[,1:5])

###loged and scaled RPKM
Gene_Expression_T_scaled <- log(Gene_Expression_T+0.25)
Gene_Expression_T_scaled <- apply(Gene_Expression_T_scaled,2,scale)

dim(Gene_Expression_T_scaled);class(Gene_Expression_T_scaled);head(Gene_Expression_T_scaled[,1:5])
rownames(Gene_Expression_T_scaled) <- row.names(Gene_Expression_T)
head(Gene_Expression_T_scaled[,1:10])
X <- Gene_Expression_T_scaled
rt<-c()
rt[1:273]<-1
rt[274:1771]<--1
Sample_class<-rt
#sh
Sample_class<-replicate(1771,-1)
Sample_class[274:550]<-1
#lb
Sample_class<-replicate(1771,-1)
Sample_class[551:813]<-1
#lt
Sample_class<-replicate(1771,-1)
Sample_class[814:1078]<-1
#kn
Sample_class<-replicate(1771,-1)
Sample_class[1079:1307]<-1
#ld
Sample_class<-replicate(1771,-1)
Sample_class[1308:1511]<-1
#ln
Sample_class<-replicate(1771,-1)
Sample_class[1512:1771]<-1

##for loop to compute the t-statistics for each  gene in a tissue
Gene_number<-15505
Myres<-matrix(nrow = 15505,ncol = 2)
for(i in 1:Gene_number){
  Y <- as.numeric(X[i,])
  myfit <- lm(Y~Sample_class)
  a <- summary(myfit)
  t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
  Myres[i,c(1,2)] <- t
}
row.names(Myres) <- row.names(X)
write.table(Myres,file = "~/Downloads/maize_GTEx/diff_expre_genes/ltldln_t.txt",sep = " ",quote = F,row.names = T,col.names = F)

rankrt<-order(Myres[,1],decreasing = TRUE)
Myres[rankrt[1]]
outlier_rt<-rownames(Myres)[rankrt[1:775]]

ranksh<-order(Myres[,1],decreasing = TRUE)
Myres[ranksh[1]]
outlier_sh<-rownames(Myres)[ranksh[1:775]]

ranklb<-order(Myres[,1],decreasing = TRUE)
Myres[ranklb[1]]
outlier_lb<-rownames(Myres)[ranklb[1:775]]

ranklt<-order(Myres[,1],decreasing = TRUE)
Myres[ranklt[1]]
outlier_lt<-rownames(Myres)[ranklt[1:775]]

rankkn<-order(Myres[,1],decreasing = TRUE)
Myres[rankkn[1]]
outlier_kn<-rownames(Myres)[rankkn[1:775]]

rankld<-order(Myres[,1],decreasing = TRUE)
Myres[rankld[1]]
outlier_ld<-rownames(Myres)[rankld[1:775]]
png(file="~/Downloads/maize_GTEx/diff_expre_genes/ldt775.png")
hist(ldt[rankld[1:775],2])
dev.off()

rankln<-order(Myres[,1],decreasing = TRUE)
Myres[rankln[1]]
outlier_ln<-rownames(Myres)[rankln[1:775]]
lnt<-read.table("~/Downloads/maize_GTEx/diff_expre_genes/ln_t.txt")
png(file="~/Downloads/maize_GTEx/diff_expre_genes/lnt775.png")
hist(ldt[rankln[1:775],2])
dev.off()

rankl3<-order(Myres[,1],decreasing = TRUE)
Myres[rankl3[1]]
outlier_l3<-rownames(Myres)[rankl3[1:775]]
# 设置文件目录
setwd("~/Downloads/maize_GTEx/图数据/gene_specific_expression/t/mr_l32")

# 每个outlier_leaf创建一个新画布
#775 5%
#465 3%
paird_out_fpkm<-exp_71
for (i in 1:640) {
  # 创建新画布 
#  title = rownames(Myres)[which(rankrt==1)]
  title = outlier_l3mr[i]
#  outlier_rt[i]<-title
  png(paste(title, ".png", sep = ""), width = 800, height = 480)
  
  # 当前图片的title和数据
  rt_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 1:273])
  sh_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 274:550])
  lb_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 551:813])
  lt_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 814:1078])
  kn_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 1079:1307])
  ld_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 1308:1511])
  ln_gene <- as.numeric(paird_out_fpkm[outlier_l3mr[i], 1512:1771])
  
  # 绘制箱线图
  boxplot(rt_gene,sh_gene, lb_gene, lt_gene, kn_gene, ld_gene,ln_gene,
          at = c(1,2,3,4,5,6,7),
          col = c("gray","gray","gray","red","gray","red","red"),
          names = c("RT","SH","LB","LT","KN","LD","LN"),
          ylab = "", main = title)
  
  # 关闭画布并保存到文件
  dev.off()
}
mr_l3<-c()
for (i in 1:15505) {
  mr_l3[i]<-min(tissue_median[i, c(4,6,7)]) / max(tissue_median[i, c(1, 2, 3, 5)])
}
topl3m<-mr_l3[rankl3[1:775]]

on_l3<-rankl3[which(topl3m>2)]
outlier_l3mr<-rownames(Myres)[on_l3]

save(outlier_rtmr,outlier_shmr,outlier_lbmr,outlier_ltmr,outlier_knmr,outlier_ldmr,outlier_lnmr,outlier_l3mr,file = '~/Downloads/maize_GTEx/diff_expre_genes/ttop5median2.RData')

png("median_rank2.png", width = 800, height = 480)
barplot(toprtm,)
abline(h=3)
abline(h=4)
dev.off()

length(which(toprtm>0&toprtm<2))
length(which(toprtm>2&toprtm<3))
length(which(toprtm>3&toprtm<4))
length(which(toprtm>4&toprtm<5))
length(which(toprtm>5&toprtm<6))
length(which(toprtm>6))

#ltldln 640gene go富集
#Problematic cache解决
library(AnnotationHub)
package = "AnnotationHub"
BiocManager::install("AnnotationHub")
oldcache = path.expand(rappdirs::user_cache_dir(appname=package))
setAnnotationHubOption("CACHE", oldcache)
ah = AnnotationHub(localHub=TRUE)
## removes old location and all resources
removeCache(ah, ask=FALSE)

## create the new default caching location
newcache = tools::R_user_dir(package, which="cache")
setAnnotationHubOption("CACHE", newcache)
#
ah = AnnotationHub()
query(ah, "zea mays")
maize <- ah[['AH111691']]
length(keys(maize))
columns(maize)
require(clusterProfiler)
idcon<-bitr(outlier_knmr, 'ALIAS', c("ENTREZID", "GO", "ONTOLOGY"), maize)
res = enrichGO(idcon[,2], OrgDb=maize, keyType = 'ENTREZID',ont='ALL',pvalueCutoff=0.05, qvalueCutoff=0.05)
save(res,file = "~/Downloads/maize_GTEx/diff_expre_genes/kngo.RData")
#
#sh<-res@result[1:35,]
library(ggplot2)
pdf(file="~/Downloads/maize_GTEx/diff_expre_genes/kngo05.pdf")
#, width = 480, height = 900
ggplot(data = res@result, aes(x = reorder(Description,-p.adjust), y = Count, fill= p.adjust))+geom_bar(stat = 'identity') +coord_flip() +
  scale_fill_gradient(low="black", high="white",
                        limits = c(0, 0.05),
                        breaks = c(0.01, 0.02, 0.03, 0.04))+
  xlab('KN')
#barplot(res, showCategory =10,color = 'p.adjust')
#  scale_fill_manual(breaks = levels(c(0,0.01,0.02,0.03,0.04,0.05)))
#,drop = TRUE, split="ONTOLOGY"
#  facet_grid(ONTOLOGY~., scale='free')
dev.off()

##kegg
library(clusterProfiler)
BiocManager::install(version = "3.16")
BiocManager::install("biocLite")
install.packages('biocLite')
biocLite("KEGG.db")
install.packages('Biostrings')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")

load('~/Downloads/maize_GTEx/diff_expre_genes/ttop5median2.RData')
columns(maize)
idcon<-bitr(outlier_rtmr, 'ALIAS', c("UNIGENE", "GID", "ENTREZID",'REFSEQ','PMID'), maize)
kk <- enrichKEGG(
  organism = "zma",
  gene = idcon[,3],
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff =0.05,
  use_internal_data=F
)
#rm(kk)
pdf(file="~/Downloads/maize_GTEx/diff_expre_genes/l3kegg05.pdf")
#, width = 480, height = 900
ggplot(data = kk@result[which(kk@result$p.adjust<0.05),], aes(x = reorder(Description,-p.adjust), y = Count, fill= p.adjust))+geom_bar(stat = 'identity') +coord_flip() +
  scale_fill_gradient(low="black", high="white",
                      limits = c(0, 0.05),
                      breaks = c(0.01, 0.02, 0.03, 0.04))+
  xlab('L3')
dev.off()

sessionInfo("clusterProfiler")
#add foldenrichment = GeneRatio/bgRatio
ego3 <- mutate(res, foldenrichment = (as.numeric(sub("/\\d+", "", res@result$GeneRatio))/as.numeric(sub("\\d+/", "", res@result$GeneRatio))) / (as.numeric(sub("/\\d+", "", res@result$BgRatio))/as.numeric(sub("\\d+/", "", res@result$BgRatio))))
go_rt5<-ego3@result[1:5,]
go_sh5<-ego3@result[1:5,]
go_ld5<-ego3@result[1:5,]
##t检验 tissue-specific gene vs remaining gene in each tissue
g15<-rownames(exp_71)
length(g15)
a<- read.csv('~/Downloads/maize_GTEx/diff_expre_genes/dndsYN.csv')
#get remaining gene
rg<-setdiff(g15,outlier_l3mr)
length(rg)
#get dnds of remaining gene
dndsa<-matrix(ncol = 2)
colnames(dndsa)<-c('query_id','dNdS')
for (i in 1:length(rg)) {
  if(length(which(a[,1]==rg[i]))!=0){
    if(length(which(a[,1]==rg[i]))==1){
      dndsa<-rbind(dndsa,a[which(a[,1]==rg[i]),c(1,6)])}
    else { 
      posmatch<-a[which(a[,1]==rg[i]),14]
      names(posmatch)<-which(a[,1]==rg[i])
      linen<-names(which(posmatch==max(posmatch)))
      dndsa<-rbind(dndsa,a[linen,c(1,6)])
    }
  }
}
dndsb<-dndsa[-which(is.na(dndsa[,2])),]
dndsd<-dndsb[-(which(dndsb[,2]==0)),]
dndse<-dndsd[-(which(dndsd[,2]>1)),]

sg<-read.csv('~/Downloads/maize_GTEx/diff_expre_genes/dndsee.csv')
t.test(dndse[,2],sg[which(sg$tissue=='RT'),3])
t.test(dndse[,2],sg[which(sg$tissue=='SH'),3])
t.test(dndse[,2],sg[which(sg$tissue=='LB'),3])
t.test(dndse[,2],sg[which(sg$tissue=='LT'),3])
t.test(dndse[,2],sg[which(sg$tissue=='LD'),3])
t.test(dndse[,2],sg[which(sg$tissue=='LN'),3])
t.test(dndse[,2],sg[which(sg$tissue=='KN'),3])
t.test(dndse[,2],sg[which(sg$tissue=='L3'),3])
