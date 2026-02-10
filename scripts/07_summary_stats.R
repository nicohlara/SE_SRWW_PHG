.libPaths( "C:/Users/nalara/AppData/Local/R/win-library/4.3")
library(glue)
library(dplyr)
library(sommer)
library(vegan)


###test similarity of Haplotype and SNP-based GRM
GRM <- read.delim(glue("{PHG_dir}/GRM.csv"), sep=",", row.names = 1)
colnames(GRM) <- rownames(GRM)
HRM <- read.delim(glue("{PHG_dir}/HRM.csv"), sep=",", row.names = 1)
colnames(HRM) <- rownames(HRM)
print(mantel(GRM, HRM))


###Get heritabilities
dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
geno_dir = "C:/Users/nalara/Documents/GitHub/SunRILs_population"
##summarize identity table
blues <- read.delim(glue("{dir}/data/blues.csv"), sep=",")
genotype <- read.bed.matrix(glue("{geno_dir}/data/SunRILs_imp_filtmerge"))
geno <- as.matrix(genotype)
geno_map <- genotype@snps[,c("id", "pos", "chr")]
A <- A.mat(geno) # additive relationship matrix
D <- D.mat(geno) # dominance relationship matrix
E <- E.mat(geno) # epistatic relationship matrix
bl <- blues |> filter(Entry %in% colnames(A))
bl$EntryD <- bl$Entry
heritability <- data.frame()
for (trait in c("HD", "WDR", "Height", "PM")) {
  ans.ADE <- mmes(as.formula(glue("{trait} ~ 1")), 
                  # random=~vsm(ism(Entry),Gu=A) + vsm(ism(EntryD),Gu=D), 
                  random=~vsm(ism(Entry),Gu=A), 
                  
                  rcov=~units, nIters=10,
                  data=bl,verbose = FALSE)
  (summary(ans.ADE)$varcomp)
  h2 <- vpredict(ans.ADE, h2 ~ (V1) / ( V1+V2) )$Estimate[1]
  # h2 <- vpredict(ans.ADE, h2 ~ (V1) / ( V1+V3) )$Estimate[1] # narrow sense
  # H2 <- vpredict(ans.ADE, h2 ~ (V1+V2) / ( V1+V2+V3) )$Estimate[1] # broad sense
  H2 <- NA
  heritability <- rbind(heritability, data.frame(trait = trait, h2 = h2, H2 = H2, prediction_limit=NA))
}
heritability$prediction_limit <- sqrt(heritability$h2)
write.csv(heritability, glue("{dir}/output/blue_heritability.csv"),quote = F, row.names = F)
