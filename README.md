This repo include analysis for RNA-seq, GWAS, GS, and epigenomics(below).
Some theorical interpretations for code can be found in my [blog](https://leilei37.github.io/2022/11/17/leilei-blog/).

### GP (Genotype preparation)

### PP (Phenotype preparation)

### pre-analysis
1. pca

### GWAS
1. calculate heritability.
2. calculate independent snp and bonferroni-p value threshold
3. GWAS
4. Fine mapping
   conditional analysis
6. colocalization(for post-GWAS to uniform peak snp)
7. LD block (to find candidate gene)
   locus zoom	  [[code]](https://github.com/leilei37/scripts/blob/main/locuszoom.R)
   
### PA (post-analysis)
1. SNP annotation	   [[code]](https://github.com/leilei37/scripts/blob/main/snpEff)
2. calculate PVE for specific cluster of snps.
3. MR
4. WGCNA
5. **evolutionary analysis** dn/ds

### GS
1. Genomic prediction(G matrix)		[[code]](https://github.com/leilei37/meta/blob/main/GS.R)
2. Prediction acuracy(cross validation)	[[code]](https://github.com/leilei37/meta/blob/main/model%20validation.R)

### epi (Epigenomics)
1. ATAC-seq
2. Bisfule-seq
   find DMR

### RNA-seq
1. Find DEG

