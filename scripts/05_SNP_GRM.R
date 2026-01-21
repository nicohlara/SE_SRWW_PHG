
library(gaston)
library(popkin)
library(glue)

dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
geno_dir="c:/Users/nalara/Documents/GitHub/SunRILs_population"

##read in SNP data
genotype <- read.bed.matrix(glue("{geno_dir}/data/SunRILs_imp_filtmerge"))

X <- t(as.matrix(genotype))
SNP_GRM <- popkin(X)
write.csv(SNP_GRM, glue("{dir}/output/GRM.csv"), quote = F)
