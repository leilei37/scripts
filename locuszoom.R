source('~/Downloads/int.R')
install.packages(c("devtools","remotes"))

#install depended packages, including ggplot2, SNPRelate, ggrepel, gdsfmt and reshape2 #ggplot2, ggrepel, and reshape2 are installed from CRAN

install.packages(c("ggplot2","ggrepel","reshape2"))

#SNPRelate and gdsfmt are installed from Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c("SNPRelate","gdsfmt"))
##### 下载IntAssoPlot from Github
library(remotes) # version 2.1.0

#download, build, and install IntAssoPlot without creating vignette

install_github("whweve/IntAssoPlot",force = T)
或者
#download, build, and install IntAssoPlot with creating vignette

install_github("whweve/IntAssoPlot",build=TRUE,build_vignettes = TRUE,force = T)


export.plink(df,filebasename = '~/Downloads/why',transpose = F,phenotypes = NULL,export012na=T)
export.merlin(df, pedfile = "~/Downloads/merlin.ped", datafile = "~/Downloads/merlin.dat", 
              mapfile = "~/Downloads/merlin.map", format = "merlin", fixstrand = "no", 
              extendedmap = TRUE, traits = 1, order = TRUE, stepids = 100, 
              dpieceFun = "new")
##
phe<-read.csv('~/Downloads/whyphe.csv')
phe2<-read.csv('~/Downloads/whyphe2.csv')

dy<-read.csv('~/Downloads/why对应.csv')
phe2<-phe2[match(dy[,2],phe2[,1]),]
phdata(df)[299:648,3]<-as.numeric(phe2[,2])
phdata(df)[299:648,3:66]<-phe2[,-1]
phdata(df)<-phdata(df)[,-67]
colnames(phdata(df))[3:66]<-colnames(phe2)[-1]
phdata(df)[1:298,3]<-NA
save(df,file = '~/Downloads/why.RData')

phdata(df)[,3:66]<-apply(phdata(df)[,3:66], 2, as.numeric)
class(phdata(df)[,3])
for (i in 1:648) {
  phdata(df)[i,67]<-mean(unlist(phdata(df)[i,c(12,21,51,63)]),na.rm = T)
}
phdata(df)[1:298,67]<-NA
which(snpnames(df)=='S1_207091226')
which(snpnames(df)=='S7_19634030')

gtdata(dfn)[,169869:169870]
dfn@gtdata<-gtdata(df)[,-(which(duplicated(snpnames(df))))]
phdata(df)<-phdata(dfn)
i<-9
data1<-as.numeric(phdata(df)[,i])
data1.gkin <- ibs(df[, autosomal(df)], weight="freq")
hy<- try(polygenic(data1,kinship.matrix =data1.gkin,data =df),silent = T)
#write.table(hy$esth2,file=paste0('/Users/leilei/Downloads/GNdata_maf001_geno01_pairld01/output/',colnames(phdata(data_rna))[i],'heritability.txt'))
if (!inherits(hy, "try-error")) {
  mm <- try(mmscore(hy, data=df), silent = T)
mmep<-mm
mm63<-mm
mm@results[which(rownames(mm@results)=='S2_397909'),] 
mm@results[which(rownames(mm@results)=='S2_49451722'),] 
mml<-mm
plot(mm)
mm
mm1<-cbind(mm@annotation[,1:2],mm@results$P1df)
write.csv(mm1,file = '~/Downloads/llp.csv',row.names = T)

library(IntAssoPlot)
association<-read.csv('~/Downloads/locuszoom/epp.csv')
association<-read.csv('~/Downloads/llp.csv')
colnames(association)<-c('Marker','Locus','Site','p')
#gtf<-read_csv('~/Downloads/LEI/data/usedanno/GCA_000005005.5_B73_RefGen_v3_genomic.gtf',header = F)
gtf<-as.matrix(GCA_000005005_5_B73_RefGen_v3_genomic)
write.table(hapmap_am368,file = '~/Downloads/locuszoom/why.hmp.txt')
hapmap_am368<-read.table('~/Downloads/locuszoom/why.hmp.txt')
#colnames(hapmap_am368)[1]<-'rs#'
colnames(hapmap_am368)<-hapmap_am368[1,]
hapmap_am368<-hapmap_am368[-1,]
asso<-association
gtf[which(gtf[,1]=='GK000031.3'),1]<-1
gtf[which(gtf[,1]=='GK000032.3'),1]<-2
gtf[which(gtf[,1]=='GK000033.3'),1]<-3
gtf[which(gtf[,1]=='CM000780.3'),1]<-4
gtf[which(gtf[,1]=='CM000781.3'),1]<-5
gtf[which(gtf[,1]=='CM000782.3'),1]<-6
gtf[which(gtf[,1]=='GK000034.3'),1]<-7
gtf[which(gtf[,1]=='CM000784.3'),1]<-8
gtf[which(gtf[,1]=='CM000785.3'),1]<-9
gtf[which(gtf[,1]=='CM000786.3'),1]<-10
rm(GCA_000005005_5_B73_RefGen_v3_genomic)
colnames(gtf)<-c('V1','V2','V3','V4','V5','V6','V7','V8','V9')
gtf<-as.data.frame(gtf)
class(hapmap$V4)
hapmap<-hapmap_am368
association1<-association[which(association$Site>left&association$Site<right),]
class(hapmap$pos)
hapmap$pos<-as.numeric(hapmap$pos)
hapmap[which(hapmap[,1]=='S2_49452364'),]

#intersect(which(hapmap$pos>left),which(hapmap$pos<right))
hapmap1<-hapmap[which(hapmap[,4]>left&hapmap[,4]<right),]
highlight<-association[,1:3]
highlight[,4]<-17
highlight[,5]<-'black'
highlight[,6]<-1.5
highlight[,7]<-'black'

colnames(highlight)<-c('rs','chrom','pos','shape','colour','size','fill')
highlight[match(c('S2_49452364','S2_49451722'),highlight$rs),4]<-16
highlight[match(c('S2_49452364','S2_49451722'),highlight$rs),5]<-'red'
highlight[match(c('S2_70066','S2_186514','S2_186514','S2_397909'),highlight$rs),4]<-16
highlight[match(c('S2_70066','S2_186514','S2_186514','S2_397909'),highlight$rs),5]<-'red'
class(maker2link$pos)
link<-highlight[match(c('S2_70066','S2_186514','S2_186514','S2_397909'),marker2highlight$rs),1:3]
link<-highlight[match(c('S2_49452364','S2_49451722'),highlight$rs),1:3]

IntRegionalPlot(chr=2,left=70000,right=400000,gtf=gtf,association=association1,hapmap=hapmap1,hapmap_ld=hapmap1,threshold=4,marker2highlight=highlight,link2gene=link,link2LD=link,leadsnpLD = FALSE,marker2label=link,label_gene_name = TRUE)

IntRegionalPlot(chr=2,left=49450000,right=49460000,gtf=gtf,association=association1,hapmap=hapmap1,hapmap_ld=hapmap1,threshold=4,marker2highlight=highlight,link2gene=link,link2LD=link,leadsnpLD = FALSE,marker2label=link,label_gene_name = TRUE)

IntRegionalPlot(chr=7,left=10695000,right=10697000,gtf=gtf,association=association1,hapmap=hapmap1,hapmap_ld=hapmap1,threshold=4,marker2highlight=highlight,link2gene=link,link2LD=link,leadsnpLD = FALSE,marker2label=link,label_gene_name = TRUE)

IntRegionalPlot(chr=3,left=229990000,right=230000000,gtf=gtf,association=association1,hapmap=hapmap1,hapmap_ld=hapmap1,threshold=4,marker2highlight=highlight,link2gene=link,link2LD=link,leadsnpLD = FALSE,marker2label=link,label_gene_name = TRUE)
