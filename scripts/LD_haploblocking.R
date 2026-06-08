##adapted from LD-based technique presented in Weber et al. 2023

library(gaston)
library(snpStats)
library(stringr)
library(purrr)

g <- genotype
chromosomes <- unique(g@snps$chr)
chrom_size <- read.delim("data/chromosome_lengths.tsv")

# LD_matrix <- matrix(nrow=9, ncol=9)
# LD_matrix[lower.tri(LD_matrix)] <- c(9,8,6,1,1,1,1,1,9,6,6,1,1,1,1,2,9,3,1,1,1,3,5,2,1,1,7,5,1,1,6,5,1,7,8,9)
# diag(LD_matrix) <- 10
# LD_matrix <- t(LD_matrix)
# threshold <- 8
# rownames(LD_matrix) <- 1:9
# colnames(LD_matrix) <- 1:9

haplo_collapse <- function(LD_matrix, threshold) {
  r <- 1
  block <- TRUE
  while (block == TRUE) {
    if (r == nrow(LD_matrix)) {
      block <- FALSE
    } else if (LD_matrix[r,r+1] >= threshold) {
      r <- r+1
    } else if (r + 1 == nrow(LD_matrix)) {
      block <- FALSE
    } else if (LD_matrix[r,r+2] >= threshold) {
      r <- r+2
    } else {
      block <- FALSE
    }
  }
  if (r > 1) { 
    return(c(rownames(LD_matrix)[1], rownames(LD_matrix)[r])) 
  } else {
      return(rownames(LD_matrix)[1])
  }
}

get_reference_range <- function(LD_region, chromosome = NA, start=FALSE, end=FALSE) {
  if (start == TRUE) { start = 1 }
  if (end == TRUE) { end = chrom_size[chrom_size$Chromosome == chromosome, 'End'] }
  
  
}


haploblocks <- list()
reference_ranges <- data.frame()
thresh <- 0.5

for (chrom in chromosomes) {
  print(chrom)
  gs <- select.snps(g, chr == chrom)
  LD <- LD(gs, lim=c(1, ncol(gs)), measure = "r2", trim=TRUE)
  infomarker <- rowSums(is.na(LD)) != ncol(LD)
  LD <- LD[infomarker,infomarker]
  LD[is.nan(LD)] <- 0
  LDs <- LD
  i <- 1
  while (nrow(LDs) >= 1) {
    block <- haplo_collapse(LDs, thresh)
    if (length(block) == 1) {
      ins <- which(rownames(LDs) == block)
    } else {
      ins <- (which(rownames(LDs) == block[1]):which(rownames(LDs) == block[2]) )
    }
    # haploblocks[[glue("{chrom}_{i}")]] <- rownames(LDs)[ins]
    haploblocks[[glue("{chrom}_{i}")]] <- c(block)
    i <- i+1
    rr <- 
    LDs <- LDs[-ins, -ins, drop=FALSE ]
  }
}

length(haploblocks)

block_info <- purrr::imap_dfr(haploblocks, function(markers, block_name) {
  chrom <- sub("_\\d+$", "", block_name)
  pos <- as.numeric(str_extract(markers, "\\d+$"))
  tibble(block = block_name,
         chromosome = chrom,
         marker_start = min(pos),
         marker_end = max(pos))
})

block_ranges <- block_info %>%
  left_join(chrom_size |> dplyr::select(Chromosome, chrom_end = End),
            by=c("chromosome" = "Chromosome")) |>
  group_by(chromosome) %>%
  arrange(marker_start, .by_group = TRUE) %>%
  mutate(start = c(1, floor((lag(marker_end)[-1] + marker_start[-1]) / 2) + 1),
         end = c(floor((marker_end[-n()] + lead(marker_start)[-n()]) / 2), chrom_end[n()])
  ) %>%
  ungroup |>
  dplyr::select(chromosome, start, end, block) |>
  mutate(chromosome = paste0("chr", chromosome),
         block = paste0("chr", block),
         c4 = 0, c5 = ".")
write.table(block_ranges, "../SE_SRWW_PHG/output/LD_block_algorithm_haploblocks.bed",
            quote=F, row.names=F, col.names=F)


