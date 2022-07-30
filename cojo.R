# all_gbl 是 所有要做GCTA cojo 的文件的全连接 "/Volumes/4TZS/yeast1/results/cojo/input/AT1Gxx.RData"
# path_cojo 是 输出文件的文件夹 例如 "/Volumes/4TZS/yeast1/results/cojo/output/"
#p0<-read.table('~/Downloads/testgrootgene.txt')
p0<-read.table('/vol3/agis/yeguoyou_group/leimengyu/results/mm/groot/testgrootgene.txt')
p0<-read.table("~/Downloads/WGCNA/rdatalist.txt")
p0<-as.character(p0[,1])
p0<-p0[-(1:84)]
p0<-p0[-(1:11)]

path_cojo<-'/vol3/agis/yeguoyou_group/leimengyu/results/cojo_input/'
path_cojo<-'~/Downloads/WGCNA/cojo_input/input/'
setwd('~/Downloads/WGCNA/')
pre_cojo <- function(all_gbl,path_cojo){
  for( i in 1:length(all_gbl)){
    #cat("preparing " , i," out of ", length(all_gbl)," ",all_gbl[i],"\n ")
    inFile1 <- all_gbl[i]
    load(inFile1)
    #emmax <- data.frame(fread(inFile1))
    q2<-summary(data_rna)
    cojo <- data.frame("SNP"=rownames(mm@results),"A1"=mm@annotation$A1,"A2"=mm@annotation$A2,"freq"=q2[rownames(mm@results),"Q.2"],
                       "b"=mm@results$effB,"se"=mm@results$se_effB,"p"=mm@results$P1df,"N"=mm@results$N)
    path_out <- paste0(path_cojo,sub(".Rdata","",basename(all_gbl)[i]),".txt")
    fwrite(cojo,file=path_out,sep="\t",quote=F,row.names = F,col.names = T)
    cat(i,"\n")
    cat(inFile1,"\n")
  }
}
require(data.table)
pre_cojo(all_gbl = list36[1],path_cojo = path_cojo)
#for check
ll<-read.table('~/Downloads/WGCNA/cojo_input/input/phekernMEgreenpc2.RData.txt')
#delete
write.table(list36,'~/Downloads/list36.txt',quote = F,row.names = F,col.names = F)
for line in `cat /Users/leilei/Downloads/list36.txt`;do
awk 'NF>7{print $0}' $line > ${line}n
done

ll<-read.table('~/Downloads/WGCNA/cojo_input/phekernMEgreenpc2.RData.txtn')
ll2<-read.table('~/Downloads/WGCNA/cojo_input/phel3bMEyellow.RData.txtn',header = T)
sn<-as.character(ll2[,1])

e<-gregexpr(pattern ='_', sn)
library(stringr)
a<-c()
i<-1
for (i in 1:4260519) {
  a[i]<-str_sub(sn[i],2,(e[[i]][length(e[[i]])])-1)
}
b<-c()
for (i in 1:4260519) {
  b[i]<-str_sub(sn[i],(e[[i]][length(e[[i]])])+1,nchar(sn[i]))
}
b<-as.numeric(b)
a<-as.numeric(a)
trans2<-cbind(sn,a,b,as.numeric(ll2[,7]))
trans2<-data.frame(trans2)
trans2[,2]<-as.numeric(trans2[,2])
trans2[,3]<-as.numeric(b)
trans2[,4]<-as.numeric(unlist(ll2[,7]))
class(trans2[,4])
colnames(trans2)[1]<-'SNP'
colnames(trans2)[2]<-'CHR'
colnames(trans2)[3]<-'BP'
colnames(trans2)[4]<-'P'
library(qqman)
png(file="~/Downloads/manhattantrans.png", width=12, height=8)
manhattan(trans2, col = c("royalblue4", "darksalmon"), suggestiveline = -log10(2.04E-8), annotateTop = T)
dev.off()

list36<-read.table('/Users/leilei/Downloads/list36.txt')
list36<-as.character(list36[,1])
cojo_gcta <- function(i){
  gcta <- "/vol3/agis/yeguoyou_group/leimengyu/src/gcta_1.93.2beta/gcta64" # gcta path
  bfile <- "/vol3/agis/yeguoyou_group/leimengyu/data/maf05" # using GenABEL to convert to plink， export.plink() function
  cojofile <- paste0("/vol3/agis/yeguoyou_group/leimengyu/results/cojo_input/",cojo_files[i])# input
  cojoout <- paste0("/vol3/agis/yeguoyou_group/leimengyu/results/cojo_output/",sub(".txt","",cojo_files[i]))# output input
  #snpfile <- "results/FT_include_snp_list.txt"
  #paste(gcta,"--bfile",bfile,"--cojo-file",cojofile,"--cojo-p 1.53e-8","--extract",snpfile,"--cojo-slct","--out",cojoout)
  cmd <- paste(gcta,"--bfile",bfile,"--cojo-file",cojofile,"--cojo-p 1.04E-8","--cojo-collinear 0.1","--cojo-wind 200","--cojo-slct","--out",cojoout)
  print(cmd)
  system(cmd)
}
cojo_gcta <- function(i){
  gcta <- "/Users/leilei/Downloads/gcta_v1.94.0Beta_macOS/gcta_v1.94.0Beta_macOS" # gcta 的卢静
  bfile <- "~/Downloads/plink_mac_20220305/maf05t" # 要把GenABEL 转换成plink， 用export.plink() 函数
  cojofile <- paste0("~/Downloads/WGCNA/cojo_input/",cojo_files[i])# 刚导出的input
  cojoout <- paste0("~/Downloads/WGCNA/cojo_output/",sub(".txt","",cojo_files[i]))# output input
  #snpfile <- "results/FT_include_snp_list.txt"
  #paste(gcta,"--bfile",bfile,"--cojo-file",cojofile,"--cojo-p 1.53e-8","--extract",snpfile,"--cojo-slct","--out",cojoout)
  cmd <- paste(gcta,"--bfile",bfile,"--cojo-file",cojofile,"--cojo-p 2E-7","--cojo-collinear 0.2","--cojo-wind 200","--cojo-slct","--out",cojoout)
  print(cmd)
  system(cmd)
}
cojo_files_all <- read.table('/vol3/agis/yeguoyou_group/leimengyu/results/cojo_input/inputlist')
cojo_files_all <- read.table('~/Downloads/list36.txt')

cojo_files_all <- read.table('~/Downloads/WGCNA/mahsnplist.txt')
cojo_files_all<-as.character(cojo_files_all[,1])
setwd('~/Downloads/WGCNA/cojo_input/')
setwd('~')
#cojo_files_all <- list.files("/vol3/agis/yeguoyou_group/leimengyu/results/cojo_input/",pattern = "dot")
#cojo_files_all <- list.files("~/Downloads/gcta_1.93.2beta",pattern = "dot")

#cojo_files <- cojo_files_all[!(cojo_files_all %in% cojo_files_dtt)]
cojo_files <- cojo_files_all
cojo_files<-list36
#cojo_files_dtt <- list.files("/Volumes/4TZS/yeast1/results/cojo/input/",pattern = "DTT")
library(parallel)
cojo_gcta(6)
r <- mclapply(X = 1:length(cojo_files),FUN = cojo_gcta, mc.cores = 2)
trans2[which(trans2[,4]<1.02E-7),]
for (i in start:end) {
  
}
###Results will be saved in a *.cma

a1<-summary(data_rna)
