install.packages("DSS")
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("R.methodsS3")
BiocManager::install("R.oo")
BiocManager::install("bsseq")
BiocManager::install("DSS")
library("DSS")
T1<-read.table('~/Downloads/epi/meth/snp10l3')
T2 = data.frame(chr = T1[,1], pos = T1[,3], N = T1[,5] + T1[,6], X = T1[,5]) # 列名必须为这四个
snp1l3<-T2
BSobj <- makeBSseqData(list(snp1rt,snp1sh,snp1l3),c('RT',"SH",'L3') )
BSobj <- makeBSseqData(list(snp1rt,snp1sh),c('RT',"SH") )

dmlTest <- DMLtest(BSobj, group1 = c("RT",'L3'), group2 = "SH",smoothing = TRUE)
dmrs    <- callDMR(dmlTest, delta = 0.1, p.threshold = 0.05)
showOneDMR(dmrs[1,], BSobj)
