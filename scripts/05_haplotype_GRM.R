.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")

options(java.parameters = "-Xmx100g")
library(rJava)
library(rPHG2)
library(glue)
library(GenomicRanges)
library(dplyr)
library(parallel)

dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"

initPhg(glue("{dir}/../phgv2_v2.4/lib"))

##create PHG database
locCon <- PHGLocalCon(as.character(glue("{dir}/vcf_dbs/hvcf_files")))
graph <- locCon |> buildHaplotypeGraph()
phgDs <- graph |> readPhgDataSet()

##characterize haploblocks
haploblocks <- phgDs |> numberOfHaplotypes(byRefRange=T)
haploblocks <- haploblocks |>
  mutate(group = sub("chr([0-9]+)[A-Z]$", "\\1", seqnames),
         subgenome = sub("chr[0-9]+([A-Z])$", "\\1", seqnames))
min(haploblocks$width); mean(haploblocks$width); max(haploblocks$width)

##get genic region names
ref_ranges <- read.delim("/project/guedira_seq_map/nico/pangenome/data/ref_ranges.bed", header = F) |>  
  dplyr::rename(start = V2, end = V3) |>
  mutate(start = start + 1)
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

##Extract haplotypes, save copy
H <- phg_genic |> readHapIds()
#write.csv(H, glue("{dir}/output/rPHG/haplotypes.tsv"), quote = F)
	
##Massage haplotypes into format best suited for fast processing
H[is.na(H)] <- "NA"
H <- H[, apply(H, 2, function(x) length(unique(x))) > 1]
Hr <- rownames(H); Hc <- colnames(H)
H <- sapply(seq_len(ncol(H)), function(j) {
  as.integer(factor(H[, j]))
})
rownames(H) <- Hr; colnames(H) <- Hc

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

print(glue("Calculating GRM for {nrow(H)} samples over {ncol(H)} haplotype regions"))
start.time <- Sys.time()
print(start.time)
HRM <- hap_GRM_calculator(H)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

write.csv(HRM, glue("{dir}/output/rPHG/HRM.csv"), quote = F)
