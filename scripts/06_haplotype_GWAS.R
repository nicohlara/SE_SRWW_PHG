.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")
library(sommer)
library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)
library(rrBLUP)

##basic implementation of an input-blind GWAS for haplotype to SNP comparisons
##inputs: BLUEs file with column 1 'ID' designating line name and X phenotype columns
##        phgv2 database haplotype calls (rPHG: phg@hapIds, write to csv)
##        SNP dataset from GASTON

##set up to be run on-server
project_dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
PHG_dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom/rPHG/"

## read in BLUEs
blues <- read.delim(glue("{project_dir}/data/blues.csv"), sep=",") |>
  dplyr::rename(ID = Entry)

##read in genotype files. Rows=ID, columns = sequence region/marker
haplotypes <- read.delim(glue("{PHG_dir}/haplotypes.tsv"))
SNPs <- read.delim(glue("{project_dir}/data/SNPs.tsv"))

##read in Relationship Matrices, rather than recalculate (especially the haplotype)
GRM <- read.delim(glue("{PHG_dir}/GRM.csv"), sep=",", row.names = 1)
colnames(GRM) <- rownames(GRM)
HRM <- read.delim(glue("{PHG_dir}/HRM.csv"), sep=",", row.names = 1)
colnames(HRM) <- rownames(HRM)

## GWAS function
##phenotype is 2-column dataframe, ID and y
##K = relationship matrix
##geno = dataframe, ID rows x genetic columns

input_agnostic_GWAS <- function(phenotype, K, geno) {
  phenotype <- phenotype |> 
    dplyr::filter(ID %in% row.names(geno))
  geno <- geno[phenotype$ID,]
  geno <- as.matrix(geno)
  geno <- geno[, apply(geno, 2, function(x) length(unique(x))) > 1]
  geno[is.na(geno)] <- "NA"
  K <- K[phenotype$ID, phenotype$ID]
  K <- as.matrix(K)
  y <- phenotype[,2]
  null <- mmer(as.formula(glue("{colnames(phenotype)[2]} ~ 1")), random = ~ vs(ID, Gu = K),
                data = phenotype)
  y_resid <- y - c(null$Beta$Estimate) - null$U$`u:ID`[[1]]
  fast_assoc <- function(y, G) {
    pvals <- numeric(ncol(G))
    for (i in seq_len(ncol(G))) {
      g <- G[, i]
      if (length(unique(g)) < 2) {
        pvals[i] <- NA
        next
      }
      fit <- lm(y ~ g)
      pvals[i] <- summary(fit)$coefficients[2, 4]
    }
    pvals
  }
  p <- fast_assoc(y_resid, geno)
  plot_df <- data.frame(id = colnames(geno), p = p)
  # plot_df$chr <- str_replace(str_replace(plot_df$id, "^S", ""), "_\\d*$", "")
  # plot_df$pos <- as.numeric(str_replace(plot_df$id, "^S\\d[ABD]_", ""))
  return(plot_df)
}


gwas_results <- data.frame()
for (trait in colnames(blues)[-c(1:2)]) {
  ph <- blues[,c('ID', trait)]
  snp_gwas <- input_agnostic_GWAS(ph, GRM, SNPs)
  snp_gwas$chr <- str_replace(str_replace(snp_gwas$id, "^S", ""), "_\\d*$", "")
  snp_gwas$pos <- as.numeric(str_replace(snp_gwas$id, "^S\\d[ABD]_", ""))
  snp_gwas$end = NA; snp_gwas$K <- "SNP"; snp_gwas$trait <- trait; 
  
  
  hap_gwas <- input_agnostic_GWAS(ph, HRM, haplotypes)
  # hap_gwas$chr <- str_replace(str_replace(hap_gwas$id, "^chr", ""), "\\:\\d*$", "")
  hap_gwas <- hap_gwas %>%
    mutate(haploregion = id) |>
    separate(haploregion, into = c("chr", "range"), sep = ":") %>%
    separate(range, into = c("pos", "end"), sep = "-") %>%
    mutate(chr = sub("^chr", "", chr))
  
  hap_gwas$K <- "Haplotype"; hap_gwas$trait <- trait 
  
  

  gwas_results <- rbind(gwas_results, snp_gwas, hap_gwas)
}

ggplot(plot_df, aes(x = pos, y = -log10(p))) +
  geom_point() +
  facet_grid(cols = vars(chr))

