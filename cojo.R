path_cojo <- '/vol3/agis/yeguoyou_group/leimengyu/results/cojo_input/'
all_gbl <- '' #mm file name(which have peak)
load('/vol3/agis/yeguoyou_group/leimengyu/data/maf5pergroot50.RData') #load gwaa data to calculate q2
pre_cojo <- function(all_gbl, path_cojo) {
  for( i in 1:length(all_gbl)){
    #cat("preparing " , i," out of ", length(all_gbl)," ",all_gbl[i],"\n ")
    inFile1 <- all_gbl[i]
    load(inFile1)
    q2 <- summary(data_rna)
    #emmax <- data.frame(fread(inFile1))
    cojo <- data.frame("SNP"=rownames(mm@results), "A1"=mm@annotation$A1, "A2"=mm@annotation$A2, "freq"=q2[rownames(mm@results), "Q.2"],
                       "b"=mm@results$effB, "se"=mm@results$se_effB, "p"=mm@results$P1df, "N"=mm@results$N)
    path_out <- paste0(path_cojo, sub(".Rdata", "", basename(all_gbl)[i]), ".txt")
    fwrite(cojo, file=path_out, sep="\t", quote=F, row.names=F, col.names=T)
    cat(i, "\n")
    cat(inFile1, "\n")
  }
}
require(data.table)
pre_cojo(all_gbl=all_gbl, path_cojo=path_cojo)
cojo_gcta <- function(i) {
  gcta <- "/vol3/agis/leimengyu/src/gcta_1.93.2beta/gcta64" # gcta path
  bfile <- "/vol3/agis/leimengyu/data/maf05" # using GenABEL to convert to plink， export.plink() function
  cojofile <- paste0("/vol3/agis/leimengyu/results/cojo_input/", cojo_files[i])# input path
  cojoout <- paste0("/vol3/agis/leimengyu/results/cojo_output/", sub(".txt", "", cojo_files[i])) # output path
  #snpfile <- "results/FT_include_snp_list.txt"
  #paste(gcta,"--bfile",bfile,"--cojo-file",cojofile,"--cojo-p 1.53e-8","--extract",snpfile,"--cojo-slct","--out",cojoout)
  cmd <- paste(gcta, "--bfile", bfile, "--cojo-file", cojofile, "--cojo-p 1.04E-8", "--cojo-collinear 0.1", "--cojo-wind 200", "--cojo-slct", "--out", cojoout)
  print(cmd)
  system(cmd)
}
cojo_files_all <- read.table('~/Downloads/WGCNA/mahsnplist.txt') #files list after pre_cojo
cojo_files <- cojo_files_all
r <- mclapply(X = 1:length(cojo_files), FUN=cojo_gcta, mc.cores=2) #parallel running
###Results will be saved in a *.cma

a1<-summary(data_rna)
