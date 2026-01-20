.libPaths( "C:/Users/nalara/AppData/Local/R/win-library/4.3")

options(java.parameters = "-Xmx12g")
library(rJava)

library(rPHG2)
library(glue)
library(GenomicRanges)
library(dplyr)

dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
initPhg(glue("{dir}/phg/lib"))


##create PHG database
locCon <- PHGLocalCon(as.character(glue("{dir}/hvcf_files")))
graph <- locCon |> buildHaplotypeGraph()
phgDs <- graph |> readPhgDataSet()

##characterize haploblocks
haploblocks <- phgDs |> numberOfHaplotypes(byRefRange=T)
haploblocks <- haploblocks |>
  mutate(group = sub("chr([0-9]+)[A-Z]$", "\\1", seqnames),
         subgenome = sub("chr[0-9]+([A-Z])$", "\\1", seqnames))
min(haploblocks$width); mean(haploblocks$width); max(haploblocks$width)

##get genic region names
ref_ranges <- read.delim(glue("{dir}/ref_ranges.bed"), header = F) |>
  dplyr::rename(start = V2, end = V3) |>
  mutate(start = start + 1)
haplo_regionames <- merge(haploblocks, ref_ranges, by = c("start", "end"))
genic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4),] |>
  dplyr::filter(seqnames != 'chrUnknown')
genic_r <- GRanges(
  seqnames = genic_ranges$seqnames,
  ranges = IRanges(genic_ranges$start, genic_ranges$end)
)

##filter PHG haplotypes
phg_genic <-  phgDs |> filterRefRanges(genic_r)

hap_GRM_calculator <- function(phg) {
  samples <- phg |> readSamples()
  n <- length(samples)
  GRM <- matrix(NA_real_, n, n)
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      if(samples[i] == samples[j]) {
        val <- 1
      } else {
        samp_haps <- phg |> filterSamples(c(samples[i], samples[j]))
        m <- samp_haps |> readHapIds()
        m2 <- m[,!is.na(m[1,]) & !is.na(m[2,]) ]
        same <- m2[,m2[1,] == m2[2,]]
        val <- ncol(same)/ncol(m2)
      }
      
      val <- compare_overlap(phg, samples[[i]], samples[[j]])
      GRM[i, j] <- val
      GRM[j, i] <- val
    }
  }
  dimnames(GRM) <- list(samples, samples)
  return(GRM)
}

HRM <- hap_GRM_calculator(phg_genic)
write.csv(HRM, glue("{dir}/output/HRM.csv"), quote = F)
