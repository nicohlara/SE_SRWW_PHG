.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")
options(java.parameters = "-Xmx100g")
library(rJava)
library(rPHG2)
library(sommer)
library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)
library(rrBLUP)
library(glue)
library(data.table)
library(GenomicRanges)

##basic implementation of an input-blind GWAS for haplotype to SNP comparisons
##inputs: BLUEs file with column 1 'ID' designating line name and X phenotype columns
##        phgv2 database haplotype calls (rPHG: phg@hapIds, write to csv)
##        SNP dataset from GASTON

##set up to be run on-server
#project_dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
project_dir="/project/guedira_seq_map/nico/pangenome"
PHG_dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"

## read in BLUEs
blues <- read.delim(glue("{project_dir}/data/blues.csv"), sep=",") |>
  dplyr::rename(ID = Entry)
print(dim(blues))
##read in genotype files. Rows=ID, columns = sequence region/marker
#haplotypes <- data.frame(fread(glue("{PHG_dir}/haplotypes.tsv")), row.names=1)
initPhg(glue("{PHG_dir}/../phgv2_v2.4/lib"))

##create PHG database
locCon <- PHGLocalCon(as.character(glue("{PHG_dir}/vcf_dbs/hvcf_files")))
graph <- locCon |> buildHaplotypeGraph()
phgDs <- graph |> readPhgDataSet()

##characterize haploblocks
haploblocks <- phgDs |> numberOfHaplotypes(byRefRange=T)
haploblocks <- haploblocks |>
  mutate(group = sub("chr([0-9]+)[A-Z]$", "\\1", seqnames),
         subgenome = sub("chr[0-9]+([A-Z])$", "\\1", seqnames))
min(haploblocks$width); mean(haploblocks$width); max(haploblocks$width)

##get genic region names
ref_ranges <- read.delim(glue("{project_dir}/data/ref_ranges.bed"), header = F) |>  
  dplyr::rename(start = V2, end = V3) |>
  dplyr::mutate(start = start + 1)
print(dim(ref_ranges))
haplo_regionames <- merge(haploblocks, ref_ranges, by = c("start", "end"))
genic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4),] |>
  dplyr::filter(seqnames != 'chrUnknown')
genic_r <- GRanges(
  seqnames = genic_ranges$seqnames,
  ranges = IRanges(genic_ranges$start, genic_ranges$end)
)

##filter PHG haplotypes to genic regions
phg_genic <-  phgDs |> filterRefRanges(genic_r)
print(length(phg_genic |> readSamples()))

##Extract haplotypes
haplotypes <- phg_genic |> readHapIds()
rownames(haplotypes) <- gsub("_G1", "", rownames(haplotypes))
haplotypes <- haplotypes[grep("SRGBS|UX", rownames(haplotypes), value=T),]
rownames(haplotypes) <- gsub("_SRGBS", "", rownames(haplotypes))
print(rownames(haplotypes)[1:25])
print(dim(haplotypes))

SNPs <- data.frame(fread(glue("{project_dir}/data/SNP.csv")), row.names=1)
print(dim(SNPs))

##read in Relationship Matrices, rather than recalculate (especially the haplotype)
GRM <- read.delim(glue("{PHG_dir}/output/rPHG/GRM.csv"), sep=",", row.names = 1)
colnames(GRM) <- rownames(GRM)
HRM <- read.delim(glue("{PHG_dir}/output/rPHG/HRM.csv"), sep=",", row.names = 1)
HRM <- HRM[grep("SRGBS|UX", rownames(HRM), value=T),grep("SRGBS|UX", colnames(HRM), value=T)]
rownames(HRM) <- gsub("_SRGBS", "", rownames(HRM))
colnames(HRM) <- rownames(HRM)
print(dim(HRM))

## GWAS function
##phenotype is 2-column dataframe, ID and y
##K = relationship matrix
##geno = dataframe, ID rows x genetic columns

input_agnostic_GWAS <- function(phenotype, K, geno, na_informative = TRUE) {
  phenotype <- phenotype |> 
    dplyr::filter(ID %in% row.names(geno))
  geno <- geno[phenotype$ID,]
  geno <- as.matrix(geno)
  geno <- geno[, apply(geno, 2, function(x) length(unique(x))) > 1]
  if (na_informative == TRUE) {geno[is.na(geno)] <- "NA"}
  K <- K[phenotype$ID, phenotype$ID]
  K <- as.matrix(K)
  print(dim(phenotype))
  print(dim(geno))
  print(dim(K))
  null <- mmer(as.formula(glue("{colnames(phenotype)[2]} ~ 1")), random = ~ vs(ID, Gu = K),
               data = phenotype)
  V  <- as.numeric(null$sigma$`u:ID`) * K + diag(as.numeric(null$sigma$units), nrow(K))
  Vi <- solve(V)
  y_t <- Vi %*% phenotype[,2]
  print("Size transformed phenotype vector")
  print(dim(y_t))
  fast_assoc <- function(y_t, G, Vi) {
    pvals <- numeric(ncol(G))
    
    for (i in seq_len(ncol(G))) {
      g <- G[,i]
      if (length(unique(g)) < 2) {
        pvals[i] <- NA
        next
      }
      if (!is.numeric(g)) {g <- as.numeric(as.factor(g))}
      g_t <- Vi %*% g
      
      fit <- lm(y_t ~ g_t)
      pvals[i] <- summary(fit)$coefficients[2,4]
    }
    pvals
  }
  p <- fast_assoc(y_t, geno, Vi)
  plot_df <- data.frame(id = colnames(geno), p = p)
  return(plot_df)
}

for (trait in colnames(blues)[-c(1:2)]) {
  print(trait)
  ph <- blues[,c('ID', trait)]
  #snp_gwas <- input_agnostic_GWAS(ph, GRM, SNPs, na_informative = FALSE)
  #print(dim(snp_gwas))
  #write.csv(snp_gwas, glue("{PHG_dir}/output/rPHG/GWAS_{trait}_SNP.csv"), row.names=F) 
  
  hap_gwas <- input_agnostic_GWAS(ph, HRM, haplotypes, na_informative = TRUE)
  print(dim(hap_gwas))
  write.csv(hap_gwas, glue("{PHG_dir}/output/rPHG/GWAS_{trait}_haplotype.csv"), row.names=F)  
}


