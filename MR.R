gdddta<-read.csv('~/Downloads/MR/trait1.csv',header = F)
nowtrait<-gdddta
bxx<-list()
byy<-list()
bb<-list()
i<-1
for (i in 1:length(nowtrait[,2])) {
  library(GenABEL)
  load(paste0('~/Downloads/LEI/data/no0maf05',as.character(nowtrait[i,4]),'.RData'))
  j<-which(colnames(phdata(data_rna))==nowtrait[i,3])
  data1<-as.numeric(phdata(data_rna)[,j])
  data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mmx<- mmscore(hy,data=data_rna)
  nowsnp<-which(row.names(mmx)==nowtrait[i,2])
  bx<-mmx@results[(nowsnp-30):(nowsnp+30),]
  library(genetics)
  now60<-gtdata(data_rna)[,(nowsnp-30):(nowsnp+30)]
  b <- LD(as.genotype(now60))$"R^2"
  bt<-t(b)
  for (h in 1:length(b)) {
    if(is.na(b[h])){b[h]<-bt[h]}
  }
  for (h in 1:length(b)) {
    if(is.na(b[h])){b[h]<-1}
  }
  load(paste0('~/Downloads/MR/',as.character(nowtrait[i,1]),'.RData'))
  by<-mm@results[(nowsnp-30):(nowsnp+30),]
  bxx[i]<-list(bx)
  byy[i]<-list(by)
  bb[i]<-list(b)
}
IVWresults<-list()
for (i in 1:length(bxx)) {
  library(MendelianRandomization)
  MRInputObject <- mr_input(bx = as.numeric(bxx[[i]][,2]),
                            bxse = as.numeric(bxx[[i]][,3]),
                            by = as.numeric(byy[[i]][,2]),
                            byse = as.numeric(byy[[i]][,3]),
                            corr = bb[[i]])
  IVWObject <- mr_ivw(MRInputObject,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = TRUE,
                      weights = "simple",
                      psi = 0,
                      distribution = "normal",
                      alpha = 0.05/13)
  IVWresults[i]<-IVWObject
}
  library(MendelianRandomization)
  MRInputObject <- mr_input(bx = as.numeric(bx[,2]),
                            bxse = as.numeric(bx[,3]),
                            by = as.numeric(by[,2]),
                            byse = as.numeric(by[,3]),
                            corr = b)
  IVWObject <- mr_ivw(MRInputObject,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = TRUE,
                      weights = "simple",
                      psi = 0,
                      distribution = "normal",
                      alpha = 0.05/13)

#by性状
by<-mm@results[(which(row.names(mm)=='S1_294014160')-30):(which(row.names(mm)=='S1_294014160')+30),]

#require(GenABEL.data)
#data(ge03d2)
# r2s using r2fast
a <- r2fast(ge03d2,snps=c(1:10))
#install.packages('MendelianRandomization')
library(MendelianRandomization)
MRInputObject <- mr_input(bx = as.numeric(bx[16:45,2]),
                          bxse = as.numeric(bx[16:45,3]),
                          by = as.numeric(by[16:45,2]),
                          byse = as.numeric(by[16:45,3]),
                          corr = b[16:45,16:45])
IVWObject <- mr_ivw(MRInputObject,
                    model = "default",
                    robust = FALSE,
                    penalized = FALSE,
                    correl = TRUE,
                    weights = "simple",
                    psi = 0,
                    distribution = "normal",
                    alpha = 0.05/13)
#IVWObject <- mr_ivw(mr_input(bx = ldlc, bxse = ldlcse,
#                             by = chdlodds, byse = chdloddsse))
rslt1<-IVWObject
calc.rho
bt<-t(b)
length(b)
for (i in 1:length(b)) {
  if(is.na(b[i])){b[i]<-bt[i]}
}
for (i in 1:length(b)) {
  if(is.na(b[i])){b[i]<-1}
}
EggerObject <- mr_egger(MRInputObject,
                        robust = FALSE,
                        penalized = FALSE,
                        correl = FALSE,
                        distribution = "normal",
                        alpha = 0.05)

EggerObject <- mr_egger(mr_input(bx = ldlc, bxse = ldlcse,
                                 by = chdlodds, byse = chdloddsse))

#gwas结果
mm@results[which(row.names(mm)=='S5_200490842'),]
mm@results[which(row.names(mm)=='S5_133167347'),]

mm@results[which(row.names(mm)=='S1_86775549'),]

mm@results[which(row.names(mm)=='S8_138416368'),]
j<-8
j<-18
j<-15
colnames(phdata(data_rna))[8]
for (j in c(15,16,17,22,25)) {
data1<-as.numeric(as.character(phdata(data_rna)[,j]))
data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
mm<- mmscore(hy,data=data_rna)
save(mm,paste0('~/Documents/gwas/',colnames(phdata(data_rna))[j],'.RData'))
}
j<-18
for (j in c(16,17,22,25)) {
  data1<-as.numeric(as.character(phdata(data_rna)[,j]))
  data1.gkin <- ibs(data_rna[, autosomal(data_rna)], weight="freq")
  hy<- polygenic(data1,kinship.matrix =data1.gkin,data =data_rna)
  mm<- mmscore(hy,data=data_rna)
  save(mm,file=paste0('~/Documents/gwas/',colnames(phdata(data_rna))[j],'.RData'))
}      
load()
      
      