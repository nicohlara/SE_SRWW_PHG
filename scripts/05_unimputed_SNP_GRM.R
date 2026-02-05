.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")

# options(java.parameters = "-Xmx100g")
# library(rJava)
# library(rPHG2)
library(glue)
# library(GenomicRanges)
library(dplyr)
library(gaston)
# library(parallel)

# dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"
geno_dir="/project/guedira_seq_map/nico/SunRILs_population_description"
# geno_dir="c:/Users/nalara/Documents/GitHub/SunRILs_population"

# blues <- read.delim(glue("{geno_dir}/data/blues.csv"), sep=",")
blues <- read.delim("/project/guedira_seq_map/nico/pangenome/data/blues.csv", sep=",")
ugeno <- read.vcf(glue("{geno_dir}/data/processed_vcf_allsequencedata_20250817/SunRILs_filtered.vcf.gz"), convert.chr=F)
pedigree <- read.csv(glue("{geno_dir}/data/cross_info.csv"), header=F, col.names = c("Cross_ID", "Parent_1", "Parent_2"))
ugeno@ped$famid <- ifelse(grepl("UX", ugeno@ped$id), str_sub(ugeno@ped$id, 1, 6), 'Parent')

##Basic filtering for unusable data based on previous parameters tested
dim(ugeno)
hist(ugeno@snps$callrate)
ugeno <- select.snps(ugeno, callrate >= 0.7)
dim(ugeno)
ugeno <- select.inds(ugeno, id %in% blues$Entry)
dim(ugeno)
hist(ugeno@snps$N1/nrow(ugeno))
ugeno <- select.snps(ugeno, N1/nrow(ugeno) < .1) 
dim(ugeno)
hist(ugeno@ped$N1/ncol(ugeno))
ugeno <- select.inds(ugeno, N1/ncol(ugeno) < 0.1)
dim(ugeno)
##LD thin to save processing power
ugeno <- LD.thin(ugeno, threshold = 0.95, max.dist=50e6)
dim(ugeno)

##Massage SNP into format best suited for fast processing
# H[is.na(H)] <- "NA"
G <- as.matrix(ugeno)
G <- G[, apply(G, 2, function(x) length(unique(x))) > 1]
Gr <- rownames(G); Gc <- colnames(G)
# H <- sapply(seq_len(ncol(H)), function(j) {
#   as.integer(factor(H[, j]))
# })
# rownames(H) <- Hr; colnames(H) <- Hc

##function for calculating haplotype GRM
hap_GRM_calculator <- function(H) {
  samples <- rownames(H)
  n <- nrow(H)
  m <- ncol(H)
  G <- matrix(0, n, n)
  
  for (k in seq_len(m)) {
    col <- H[, k]
    ok <- !is.na(col)
    if (sum(ok) < 2) next
    
    Z <- model.matrix(~ factor(col[ok]) - 1)
    G[ok, ok] <- G[ok, ok] + tcrossprod(Z)
  }
  
  G <- G / m 
  rownames(G) <- samples; colnames(G) <- samples
  return(G)
}

print(glue("Calculating GRM for {nrow(G)} samples over {ncol(G)} haplotype regions"))
start.time <- Sys.time()
print(start.time)
GRM <- hap_GRM_calculator(G)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

write.csv(GRM, glue("{dir}/output/rPHG/uSNP_GRM.csv"), quote = F)
