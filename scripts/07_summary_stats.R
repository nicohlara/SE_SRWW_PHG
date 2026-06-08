#.libPaths( "C:/Users/nalara/AppData/Local/R/win-library/4.3")
.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")
library(glue)
library(dplyr)
library(sommer)
library(vegan)
library(stringr)

PHG_dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"

##read in Relationship Matrices
#GRM <- read.delim(glue("{PHG_dir}/output/rPHG/GRM.csv"), sep=",", row.names = 1)
group_order <- c("UX1989", "UX1991", "UX1994", "UX2029", 
           "UX1992", "UX1993", "UX2010", "UX2012", "UX2013", "UX2023", "UX2026", 
           "UX1995", "UX1997", "UX2000", "UX2031")

###test similarity of Haplotype and SNP-based GRM
GRM <- read.delim(glue("{PHG_dir}/GRM.csv"), sep=",", row.names = 1)
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

sharecol <- grep("UX", intersect(rownames(GRM), rownames(HRM)), value=T)
# sharecol <- sharecol[order(sharecol)]

sharecol <- sharecol[order(factor(sub("-.*", "", sharecol), levels = group_order),
        as.numeric(sub(".*-", "", sharecol)))]


GRM <- GRM[sharecol, sharecol]
HRM <- HRM[sharecol, sharecol]

GRM <- as.matrix(GRM)
HRM <- as.matrix(HRM)

print(mantel(GRM, HRM))
cor(diag(GRM), diag(HRM))

min(as.vector(GRM)); max(as.vector(GRM)); 
min(as.vector(HRM)); max(as.vector(HRM)); 

res_val=3
png(filename = paste0(glue("{PHG_dir}/../figures/kinship_map.png")),  
    width=750*res_val, height=300*res_val, res=72*res_val,
    bg='transparent')
nf <- layout(matrix(c(1,2,3), 1,3 ,byrow=T), widths=c(rep(ncol(GRM)+2, 2), ncol(GRM)/10), heights = c(nrow(GRM)+2), respect=T)
par(mar=c(1.5,1,2,0))
image(GRM , col=monochrome_palette, axes=F)
title("SNP GRM")
image(HRM , col=monochrome_palette, axes=F)
title("Haplotype GRM")
image(t(matrix(1:100)), col=monochrome_palette, axes=F)
axis(side=4, at=seq(0, 1, length.out=10), labels=seq(.1,1, length.out=10))#, col.axis=dark_hue)
mtext('Degree of relatedness', side=4, line=2.5, at=.5,  srt=270)#, col=dark_hue)
dev.off()

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


###evaluate significant QTL
sig_hap <- read.delim("output/significant_haplotype_calls.csv", sep=",", row.names=1)

hapassign_percent <- function(founder, vector) {vec <- grepl(founder, vector); return(length(vec[vec == T]))}

crosses <- c("UX1989", "UX1991", "UX1992", "UX1993", "UX1994",
             "UX1995", "UX1997", "UX2000", "UX2010", "UX2012",
             "UX2013", "UX2023", "UX2026", "UX2029", "UX2031")
crosses <- c("UX1992", "UX1993", "UX2010", "UX2013", "UX2023", "UX2026")
for (fam in crosses) {
  print(fam)
  h <- sig_hap[grep(fam, row.names(sig_hap)),]
  print(dim(h))
  hv <- as.vector(as.matrix(h))
  if (fam == "UX2023") {
    print("AGS2000 haplotype percentage:")
    print(hapassign_percent("AGS2000", hv)/length(hv))
  }
  print("HILLIARD haplotype percentage:")
  print(round(hapassign_percent("HILLIARD", hv)/length(hv), 2))
  print("")
}

subpop <- "UX2029"
chrom <- "chr4D"
gene_exam <- sig_hap[grep(subpop, row.names(sig_hap)), grep(chrom, colnames(sig_hap))]
##remove fixed haplotypes
gene_exam <- gene_exam[,as.vector(apply(gene_exam, 2, function(x) length(unique(x))) > 1)]
gene_exam[,"chr4D.18875554.18886103"]


##zoom in on the Ppd-D1 region
haplotype_collapse <- function(haplotypes) {
  tb <- table(haplotypes)
  tb <- round(tb/sum(tb), 2)
  paste0(names(tb), ":", tb, collapse = "/")
}

##select reference ranges
zoom_region <- sig_hap[,c("chr2D.31411462.31415806",
                          "chr2D.31698093.31705264",
                          "chr2D.31485532.31490516",
                          "chr2D.36627456.36629723",
                          "chr2D.36203982.36210203")]
imputed_lines <- gsub("_SRGBS", "", rownames(zoom_region))

##shorten haplotype names
# zoom_region[] <- lapply(zoom_region, function(x) gsub("AGS2000", "A", x))
replacements <- c(
  "AGS2000" = "A",
  "CLARK" = "C",
  "HILLIARD" = "H",
  "IL0218228" = "I",
  "LA03136E71" = "L",
  "MO080104" = "M",
  "_" = ""
)

zoom_region[] <- lapply(zoom_region, function(x) {
  str_replace_all(x, replacements)
})
zoom_region$Entry <- imputed_lines
zoom_region$Cross_ID <- ifelse(grepl("UX", imputed_lines), 
                               str_sub(imputed_lines, 1, 6), 
                               imputed_lines)

zoom_region <- zoom_region |>
  filter(!grepl("UX", Cross_ID) | Entry %in% blues$Entry)

# zoom_region$Entry <- rownames(zoom_region)
# zoom_region <- merge(zoom_region, blues[c("Cross_ID", "Entry")], by="Entry")
zoom_table <- zoom_region |> 
  select(-Entry) |>
  group_by(Cross_ID) |>
  summarize_all((haplotype_collapse)) |>
  arrange(Cross_ID) 
write.csv(zoom_table, "output/haplotype_regional_summary_PpdD1.csv",
          row.names=F)

