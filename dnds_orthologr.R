#this script including
#dnds calculate by orthologr package
#t-test across tissues using dnds
#visualization typical expression across tissue, go&kegg and t test results
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# Install package dependencies
BiocManager::install(c("Biostrings", "GenomicRanges", "GenomicFeatures", "Rsamtools", "rtracklayer"))

# install orthologr from GitHub
devtools::install_github("HajkD/metablastr")

# install orthologr from GitHub
devtools::install_github("HajkD/orthologr")
#
library(orthologr)
# get a dNdS table using:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single core processing 'comp_cores = 1'
read.cds('/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/zvtest.fasta',format = 'fasta')
read.cds('/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/b73test.txt',format = 'fasta')
test <- Biostrings::readDNAStringSet("/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/zvtest.fasta")
test
test <- Biostrings::readDNAStringSet("/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/b73test.fasta")

a<-dNdS(query_file      = '/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/maizeb73v3_cds1.fasta',
     subject_file    = '/Users/leilei/Downloads/maize_GTEx/diff_expre_genes/Zv-TIL01-REFERENCE-PanAnd-1.0_Zv00001aa.1.cds.fa',
     ortho_detection = "RBH", 
     blast_path      = "/Users/leilei/Downloads/packages/ncbi-blast-2.14.0+/bin/",
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     aa_aln_path     = "/Users/leilei/Downloads/packages/clustalw-2.1-macosx/",
     codon_aln_tool  = "pal2nal", 
     kaks_calc_path  = "/Users/leilei/Downloads/packages/KaKs_Calculator2.0/src/KaKs_Calculator",
     dnds_est.method = "Comeron", 
     comp_cores      = 6)
#各组织dnds计算
rt<-cbind(outlier_rtmr,'RT')
sh<-cbind(outlier_shmr,'SH')
lb<-cbind(outlier_lbmr,'LB')
lt<-cbind(outlier_ltmr,'LT')
kn<-cbind(outlier_knmr,'KN')
ld<-cbind(outlier_ldmr,'LD')
ln<-cbind(outlier_lnmr,'LN')
l3<-cbind(outlier_l3mr,'L3')
colnames(l3)<-c('gene','tissue')
dndsgene<-rbind(rt,sh,lb,lt,kn,ld,ln,l3)
write.csv(a,'~/Downloads/maize_GTEx/diff_expre_genes/dndsYN.csv')
write.csv(dndsgene,'~/Downloads/maize_GTEx/diff_expre_genes/dndsgene.csv')
a<- read.csv('~/Downloads/maize_GTEx/diff_expre_genes/dnds.csv')
a<- read.csv('~/Downloads/maize_GTEx/diff_expre_genes/dndsYN.csv')
dndsa<-matrix(ncol = 3)
colnames(dndsa)<-c('query_id','dNdS','as.character(dndsgene[i, 2])')
i<-70
for (i in 1:length(dndsgene[,1])) {
  if(length(which(a[,1]==dndsgene[i,1]))!=0){
    if(length(which(a[,1]==dndsgene[i,1]))==1){
      dndsa<-rbind(dndsa,cbind(a[which(a[,1]==dndsgene[i,1]),c(1,6)],as.character(dndsgene[i,2])))}
    else { 
      posmatch<-a[which(a[,1]==dndsgene[i,1]),14]
      names(posmatch)<-which(a[,1]==dndsgene[i,1])
      linen<-names(which(posmatch==max(posmatch)))
      dndsa<-rbind(dndsa,cbind(a[linen,c(1,6)],as.character(dndsgene[i,2])))
      }
    }
  }
#boxplot
class(dndsa[,2])
colnames(dndsa)<-c('gene','dnds','tissue')
#根据比对百分比 删除1241 GRMZM2G103900
dndsa<-dndsa[-1241,]
dndsb<-dndsa[-which(is.na(dndsa[,2])),]
#dndsb<-dndsb[-which(dndsb[,1]=='GRMZM2G050982'),]

#dndsc<-dndsb[-(which(dndsb[,2]>1)),]

dndsd<-dndsb[-(which(dndsb[,2]==0)),]
dndse<-dndsd[-(which(dndsd[,2]>1)),]
dndse<-dndse[-(which(duplicated(dndse[,c(1,3)]))),]
write.csv(dndse,file = '~/Downloads/maize_GTEx/diff_expre_genes/dndsee.csv')

#boxplot(dndsb[which(dndsb$tissue=='RT'),2])
install.packages('ggsignif')
library(ggsignif)
png(file="~/Downloads/maize_GTEx/diff_expre_genes/dndsYNee.png", width = 900, height = 900)
pdf(file="~/Downloads/maize_GTEx/diff_expre_genes/dndsYNee.pdf")
pdf(file="~/Downloads/maize_GTEx/图数据/gene_specific_expression/dndsYNee.pdf")
ggplot(dndse,aes(x=factor(dndse$tissue,levels = c('LD','LN','LT','L3','RT','SH','LB','KN')),y=dnds))+
  geom_boxplot(width=0.6,alpha=0.8,outlier.shape = NA)+theme_classic()+
  stat_boxplot(geom = "errorbar",
               aes(ymin = ..ymax..), 
               width = 0.2, size = .3) +
  stat_boxplot(geom = "errorbar",
               aes(ymax = ..ymin..), width = 0.2, size = .3)+
  geom_signif(comparisons = list(c('RT','L3')),
#                                 c('LD','RT'),
#                                 c('LD','SH'),
#                                 c('LD','LB'),
#                                 c('LD','LT'),
#                                 c('LD','L3'),
#                                 c('LD','KN'),
#                                 c('LD','LN')),
              test = 't.test')
#            map_signif_level = T)
#  +geom_hline(aes(yintercept=0.15))
dev.off()

ldt<-read.table("~/Downloads/maize_GTEx/diff_expre_genes/ltldln_t.txt")
png(file="~/Downloads/maize_GTEx/diff_expre_genes/ldt.png", width = 900, height = 900)
hist(ldt[,2])
dev.off()

dndse<-sg
t.test(dndse$dnds[which(dndse$tissue=='L3')],dndse$dnds[which(dndse$tissue=='RT')])
#LD对所有都显著
#RT和L3显著<0.01

##fig2 
columns(maize)
tsg<-read.csv('~/Downloads/maize_GTEx/diff_expre_genes/dndsgene.csv',row.names = 1)
intersect(tsg$gene,idcon$ALIAS)
setdiff(idcon$ALIAS,tsg$gene)
unique(idcon$ALIAS)
unique(tsg$gene)

idcon<-bitr(tsg$gene, 'ALIAS', c("ENTREZID", "EVIDENCE", "SYMBOL",'GENENAME'), maize,drop = F)
sst<-merge(tsg,idcon,by.x='gene',by.y='ALIAS',all.x=T)
write.csv(sst,file = '~/Downloads/maize_GTEx/diff_expre_genes/supplementaltable_tissuespecificgene.csv')
dndsganno<-idcon

title = 'GRMZM2G045318'
i<-which(rownames(exp_71)==title)
pdf(file="~/Downloads/maize_GTEx/图数据/gene_specific_expression/knexpression.pdf",height = 3,width = 6)
paird_out_fpkm<-exp_71
#for (i in 1:640) {
  # 创建新画布 
  #  title = rownames(Myres)[which(rankrt==1)]
  #  outlier_rt[i]<-title
#  png(paste(title, ".png", sep = ""), width = 800, height = 480)
  
  # 当前图片的title和数据
  rt_gene <- as.numeric(paird_out_fpkm[i, 1:273])
  sh_gene <- as.numeric(paird_out_fpkm[i, 274:550])
  lb_gene <- as.numeric(paird_out_fpkm[i, 551:813])
  lt_gene <- as.numeric(paird_out_fpkm[i, 814:1078])
  kn_gene <- as.numeric(paird_out_fpkm[i, 1079:1307])
  ld_gene <- as.numeric(paird_out_fpkm[i, 1308:1511])
  ln_gene <- as.numeric(paird_out_fpkm[i, 1512:1771])
  
  # 绘制箱线图
  boxplot(rt_gene,sh_gene, lb_gene, lt_gene, kn_gene, ld_gene,ln_gene,
          at = c(1,2,3,4,5,6,7),width=c(0.5,0.5,0.5,0.5,0.5,0.5,1),
          col = c("gray","gray","gray",'gray',"red","gray","gray"),
          names = c("RT","SH","LB","LT","KN","LD","LN"),
          ylab = "", main = title)
  
  # 关闭画布并保存到文件
  dev.off()
#}
#fig2 c
rtgo<-res@result[1:5,]
shgo<-res@result[1:5,]
kngo<-res@result[1:5,]
l3go<-res@result[1:5,]
go4<-rbind(rtgo,shgo,kngo,l3go)
rtkk<-kk@result[1:3,]
shkk<-kk@result[1:3,]
knkk<-kk@result[1:3,]
l3kk<-kk@result[1:3,]
go4<-rbind(rtgo[1:3,-1],rtkk,shgo[1:3,-1],shkk,kngo[1:3,-1],knkk,l3go[1:3,-1],l3kk)

ego3 <- mutate(go4, foldenrichment = (as.numeric(sub("/\\d+", "", go4$GeneRatio))/as.numeric(sub("\\d+/", "", go4$GeneRatio))) / (as.numeric(sub("/\\d+", "", go4$BgRatio))/as.numeric(sub("\\d+/", "", go4$BgRatio))))
ego3<-cbind(ego3,1:24)
ego3<-cbind(ego3,rep(c(1,1,1,2,2,2),4))
colnames(ego3)[11]<-'number'
ego3[23,2]<-'Carbon metabolism2'
write.csv(ego3,file = '~/Downloads/maize_GTEx/diff_expre_genes/gokegg6x4.csv')
pdf(file="~/Downloads/maize_GTEx/diff_expre_genes/fig2ckk2.pdf",width = 7,height = 5)
#, width = 480, height = 900
#749BCE blue
#FBB3B0 pink
ggplot(data = ego3, aes(x = reorder(Description,-number), y = foldenrichment,fill= p.adjust))+geom_bar(stat = 'identity') +coord_flip() +
  scale_fill_gradient(low="#FBB3B0", high="white",
                      limits = c(0, 0.05),
                      breaks = c(0.01, 0.02, 0.03, 0.04))+theme(axis.text = element_text(size = 13))
#  xlab('RT')
#barplot(res, showCategory =10,color = 'p.adjust')
#  scale_fill_manual(breaks = levels(c(0,0.01,0.02,0.03,0.04,0.05)))
#,drop = TRUE, split="ONTOLOGY"
#  facet_grid(ONTOLOGY~., scale='free')
dev.off()
