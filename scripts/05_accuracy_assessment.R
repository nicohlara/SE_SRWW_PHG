.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")

options(java.parameters = "-Xmx100g")
library(rJava)

library(rPHG2)
library(glue)
library(GenomicRanges)
library(dplyr)

dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"
initPhg(glue("{dir}/../phgv2_v2.4/lib"))


##create PHG database
locCon <- PHGLocalCon(as.character(glue("{dir}/vcf_dbs/hvcf_files")))
graph <- locCon |> buildHaplotypeGraph()
phgDs <- graph |> readPhgDataSet()

##subset down to parent lines
samples <- phgDs |> readSamples()
samples <- samples[grep("UX", samples, invert=T)]
phgDs <- phgDs |> filterSamples(c(samples))
print(phgDs |> readSamples())
print(dim(phgDs |> readHapIds()))

##filter out fixed regions
multimorphic <- phgDs |> numberOfHaplotypes(byRefRange=T) |> dplyr::filter(n_haplo > 1)
ranges <- multimorphic |>
  dplyr::filter(seqnames != 'chrUnknown')
gr <- GRanges(
  seqnames = ranges$seqnames,
  ranges = IRanges(ranges$start, ranges$end)
  )
phgDs <-  phgDs |> filterRefRanges(gr)
print(dim(phgDs |> readHapIds()))

##filter to genic regions
haploblocks <- phgDs |> numberOfHaplotypes(byRefRange=T)
haploblocks <- haploblocks |>
  mutate(group = sub("chr([0-9]+)[A-Z]$", "\\1", seqnames),
         subgenome = sub("chr[0-9]+([A-Z])$", "\\1", seqnames))
ref_ranges <- read.delim("/project/guedira_seq_map/nico/pangenome/data/ref_ranges.bed", header = F) |>
  dplyr::rename(start = V2, end = V3) |>
  mutate(start = start + 1)
print(dim(ref_ranges))

haplo_regionames <- merge(haploblocks, ref_ranges, by = c("start", "end"))
genic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4),] |>
  dplyr::filter(seqnames != 'chrUnknown')
genic_r <- GRanges(
  seqnames = genic_ranges$seqnames,
  ranges = IRanges(genic_ranges$start, genic_ranges$end)
)
phg_genic <-  phgDs |> filterRefRanges(genic_r)
print("Genic PHG size:")
print(dim(phg_genic |> readHapIds()))

##create haplotype overlap table
compare_overlap <- function(phg, sample1, sample2) {
  if(sample1 == sample2) {return(1)}
  samp_haps <- phg |> filterSamples(c(sample1, sample2))
  m <- samp_haps |> readHapIds()
  m2 <- m[,!is.na(m[1,]) & !is.na(m[2,]) ]
  same <- m2[,m2[1,] == m2[2,]]
  return(ncol(same)/ncol(m2))
}

exome <- grep("exm", phgDs |> readSamples(), value=T)
gbs <- grep("GBS", phgDs |> readSamples(), value=T, ignore.case=T) 
founders <- grep("exm|GBS", phgDs |> readSamples(), value=T, invert=T, ignore.case=T)

identity_table <- data.frame(seq_sample = c(exome, gbs))
for (sample in founders) {
  print(sample)
  comp_vec_genic <- c()
  for (shortread in c(exome, gbs)) {
    print(shortread)
    comp_vec_genic <- c(comp_vec_genic, compare_overlap(phg_genic, sample, shortread))
    print(tail(comp_vec_genic, 1))
  }
  identity_table[[glue("{sample}_genic")]] = round(as.numeric(comp_vec_genic), 2)
}
print("Genic identity table calculated")
print(dim(identity_table))
write.csv(identity_table, glue("{dir}/output/rPHG/haplotype_identity_table.csv"), quote = F, row.names=F)


##calculate intergenic ranges
intergenic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4, invert=T),] |>
  dplyr::filter(seqnames != 'chrUnknown')
intergenic_r <- GRanges(
  seqnames = intergenic_ranges$seqnames,
  ranges = IRanges(intergenic_ranges$start, intergenic_ranges$end)
)
phg_intergenic <-  phgDs |> filterRefRanges(intergenic_r)
print("Intergenic PHG size:")
print(dim(phg_intergenic |> readHapIds()))

identity_table <- data.frame(seq_sample = c(exome, gbs))
for (sample in founders) {
  print(sample)
  comp_vec_intergenic <- c()
  for (shortread in c(exome, gbs)) {
    print(shortread)
    comp_vec_intergenic <- c(comp_vec_intergenic, compare_overlap(phg_intergenic, sample, shortread))
    print(tail(comp_vec_intergenic, 1))
  }
  identity_table[[glue("{sample}_intergenic")]] = round(as.numeric(comp_vec_intergenic), 2)
}
print("Intergenic identity table calculated")
print(dim(identity_table))
write.csv(identity_table, glue("{dir}/output/rPHG/haplotype_identity_table_intergenic.csv"), quote = F, row.names=F)

phg_subset <- phg_genic

phg_subset@hapIds[,1:9]

remap_haplotype <- function(hap_vec, founders) {
  founder_calls <- hap_vec[founders]

  # unique hex codes among founders
  uniq_hex <- unique(founder_calls)

  # build mapping: hex -> label made of founder names
  mapping <- vapply(uniq_hex, function(h) {
    fns <- names(founder_calls[founder_calls == h])
    paste(fns[!is.na(fns)], collapse = "_")
  }, FUN.VALUE = character(1))

  # replace across all samples
  mapped <- mapping[hap_vec]
  return(mapped)
}

hapmat_remapped <- apply(phg_subset@hapIds, 2, remap_haplotype, founders = paste0(founders, "_G1"))
rownames(hapmat_remapped) <- rownames(phg_subset@hapIds)

phg_stats <- list()
for (line in c("AGS2000", "HILLIARD", "CLARK", "IL02", "LA03136", "MO080104" )) {

  exome_line <-  hapmat_remapped[grep(glue("exm{line}"), rownames(hapmat_remapped), value=T),]
  ex_table <- table(exome_line)
  ontype <- round(sum(ex_table[grep(line, names(ex_table))]/length(exome_line)),2)
  phg_stats[[glue("{line}_haplotype_match_fraction")]] <- ontype
  ex_table2 <- ex_table[grep(line, names(ex_table), invert=T)]/length(exome_line)
  offtype <- round(ex_table2[ex_table2 > 0.005], 2)
  phg_stats[[glue("{line}_offtype_haplotype_attribution_greater_than_0.005")]] <- offtype
  phg_stats[[glue("{line}_fraction_uncharacterized")]] <- 1-(ontype+sum(offtype))
}
sink(glue("{dir}/output/rPHG/haplotype_comparison_breakdown.txt"))
print(phg_stats)
sink()
