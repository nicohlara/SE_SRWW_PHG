#.libPaths( "C:/Users/nalara/AppData/Local/R/win-library/4.3")
.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")
library(glue)
library(dplyr)
library(sommer)
library(vegan)

PHG_dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"

##read in Relationship Matrices
GRM <- read.delim(glue("{PHG_dir}/output/rPHG/GRM.csv"), sep=",", row.names = 1)
colnames(GRM) <- rownames(GRM)
print(dim(GRM))
HRM <- read.delim(glue("{PHG_dir}/output/rPHG/HRM.csv"), sep=",", row.names = 1)
HRM <- HRM[grep("SRGBS|UX", rownames(HRM), value=T),grep("SRGBS|UX", colnames(HRM), value=T)]
rownames(HRM) <- gsub("_SRGBS", "", rownames(HRM))
colnames(HRM) <- rownames(HRM)
print(dim(HRM))
uGRM <- read.delim(glue("{PHG_dir}/output/rPHG/uSNP_GRM.csv"), sep=",", row.names = 1)
colnames(uGRM) <- rownames(uGRM)
print(dim(uGRM))

overlap_samples <- intersect(intersect(rownames(GRM),rownames(HRM)),rownames(uGRM))
GRM <- GRM[overlap_samples, overlap_samples]
HRM <- HRM[overlap_samples, overlap_samples]
uGRM <- uGRM[overlap_samples, overlap_samples]
print(dim(GRM)); print(dim(HRM)); print(dim(uGRM))

print("Comparing GRM and HRM")
print(mantel(HRM, GRM))
print("Comparing HRM and uGRM")
print(mantel(uGRM, HRM))
print("Comparing GRM and uGRM")
print(mantel(uGRM, GRM))


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
