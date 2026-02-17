.libPaths("/home/nicolas.lara/R/x86_64-pc-linux-gnu-library/4.4")
options(java.parameters = "-Xmx100g")
library(rJava)
library(rPHG2)
# library(dplyr)
# library(stringr)
library(glue)
# library(data.table)
library(GenomicRanges)

##set up to be run on-server
#project_dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
project_dir="/project/guedira_seq_map/nico/pangenome"
PHG_dir="/90daydata/guedira_seq_map/nico2/pangenome_multichrom"

##initialize PHG
initPhg(glue("{PHG_dir}/../phgv2_v2.4/lib"))
locCon <- PHGLocalCon(as.character(glue("{PHG_dir}/vcf_dbs/hvcf_files")))
graph <- locCon |> buildHaplotypeGraph()
phgDs <- graph |> readPhgDataSet()

founders <- grep("exm|gbs|UX|GBS", phgDs |> readSamples(), value=T, invert=T)


##generate number of haplotypes per haplotype region, summarize
multimorphic <- phgDs |> numberOfHaplotypes(byRefRange=T) |>
  dplyr::filter(seqnames != 'chrUnknown')
print(table(multimorphic$n_haplo))

##extract only relevant regions
significant_regions <- read.delim(glue("{project_dir}/output/significant_haplotypes.csv"), sep=",")
gr <- GRanges(
  seqnames = significant_regions$seqnames,
  ranges = IRanges(significant_regions$start, significant_regions$end)
)
PHG_filter <-  phgDs |> filterRefRanges(gr)
PHG_filter <- PHG_filter |> filterSamples(c(founders,grep("GBS|UX", phgDs |> readSamples(), value=T)))

##rename haplotypes to equal founder calls for easier human analysis
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
  mapped <- gsub("_G1", "", mapped)
  return(mapped)
}


hapmat_remapped <- apply(PHG_filter@hapIds, 2, remap_haplotype, founders = paste0(founders, "_G1"))
rownames(hapmat_remapped) <- gsub("_G1", "", rownames(PHG_filter@hapIds))
write.csv(hapmat_remapped, glue("{PHG_dir}/output/rPHG/significant_haplotype_calls.csv"), row.names=T)
