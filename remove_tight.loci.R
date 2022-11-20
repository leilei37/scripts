
remove_tight.loci <- function(SNPs=full_name_g,genable=data,cut_dis=10e3,cut_r2=0.9){
  ############cut_dis=10e6
  if(!require(igraph))
    require(igraph)
  if(!require(GenABEL))
    require(GenABEL)
  # rule1 physical distance
  pos.all <- map(genable[,unique(SNPs)])
  names(pos.all) <- unique(SNPs)
  pos <- pos.all[SNPs]
  diffs <- abs(outer(pos, pos, FUN = "-")) #all pairwise differences in physical position
  diffs[lower.tri(diffs)] <- 1000000
  diag(diffs) <- 1000000
  #identical(rownames(diffs),rownames(r2))
  ## rule 2 linkage
  r2 <- r2fast(data = genable,snpsubset = SNPs)
  diag(r2) <- 0
  r2[lower.tri(r2)] <- 0
  ## rule3 chromsome
#  chrs <- gsub(pattern = ".*_.*_.*_.*",replacement = "\\1",x = SNPs)
  e<-gregexpr(pattern ='_', SNPs)
  library(stringr)
  a<-c()
  for (i in 1:length(SNPs)) {
    a[i]<-str_sub(SNPs[i],2,(e[[i]][length(e[[i]])])-1)
  }
  chrs<-a
#  chrs <- gsub(pattern = ".*_(chr.*)_.*_.*_.*",replacement = "\\1",x = SNPs)
  chrs.order <- match(chrs,unique(chrs))
  diffs.chr <- abs(outer(chrs.order, chrs.order, FUN = "-")) #all pairwise differences in index
  diffs.chr[lower.tri(diffs.chr)] <- 1
  diag(diffs.chr) <- 1
  
  clump <- which(diffs < cut_dis & r2 > cut_r2 & diffs.chr==0) #differences smaller than n
  if(length(clump) <=1)
    stop("nothing")
  ## diffs.col here is wrong should be row and should be clump %% nrow(diffs)
  ## diffs.row should be ceiling(clump/ncol(diffs))
  diffs.col <- ceiling(clump/ncol(diffs)) #differences smaller than n, columns in the distance matrix
  diffs.row <- clump - (diffs.col - 1)*nrow(diffs) #differences smaller than n, rows in the distance matrix
  #mycol <- clump %% nrow(diffs)
  #myrow <- ceiling(clump/ncol(diffs))
  
  #require(igraph)
  diffs.graph <- graph_from_data_frame(data.frame(diffs.col, diffs.row))
  diffs.clust <- clusters(diffs.graph)
  
  snps <- SNPs
  relation <- list()
  name <- c()
  for(i in 1:diffs.clust$no){
    nodesInCluster <- as.numeric(names(diffs.clust$membership[diffs.clust$membership == i]))
    #chrs[nodesInCluster] <- chrs[nodesInCluster[1]]
    relation[[i]] <- snps[nodesInCluster]
    snps[nodesInCluster] <- snps[nodesInCluster[1]]
    name <- c(name, snps[nodesInCluster[1]])
    #pos[nodesInCluster] <- pos[nodesInCluster[1]]
  }
  names(relation) <- name
  if(identical(snps,SNPs)){
    warning("no change ")
    return(list("snp"=snps,"re"=relation))
  }else{
    old <- SNPs#c()
    new <- snps#c()
    # for(i in 1:length(relation)){
    #   list_now <- relation[[i]]
    #   list_now_name <- names(relation)[i]
    #   new <- c(new,rep(list_now_name,length(list_now)))
    #   old <- c(old,list_now)
    # }
    return(list("snp"=snps,"re"=relation,"cor"=data.frame("new"=new,"old"=old)))
  }
}
