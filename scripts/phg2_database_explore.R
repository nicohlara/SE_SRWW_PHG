
# install.packages("pak")
# pak::pak("maize-genetics/rPHG2")

options(java.parameters = "-Xmx12g")
library(rJava)

library(rPHG2)
library(glue)
library(GenomicRanges)
library(dplyr)

dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"
initPhg(glue("{dir}/phg/lib"))

chrom_pos <- read.delim(glue("{dir}/data/chrom_info.txt"))

locCon <- PHGLocalCon(as.character(glue("{dir}/vcf_dbs/hvcf_files")))

graph <- locCon |> buildHaplotypeGraph()

##create dataset, examine
phgDs <- graph |> readPhgDataSet()

phgDs |> numberOfHaplotypes(byRefRange=T) |> dplyr::filter(n_haplo > 1)

##percent monomorphic
multimorphic <- phgDs |> numberOfHaplotypes(byRefRange=T) |> dplyr::filter(n_haplo > 1)
nrow(phgDs |> numberOfHaplotypes(byRefRange=T) |> dplyr::filter(n_haplo <= 1))/(phgDs |> numberOfRefRanges())
##visualize graph
phgDs |> plotHaploCounts(geom="p")


#multimorphic ranges
ranges <- multimorphic |>
  dplyr::filter(seqnames != 'chrUnknown')

gr <- GRanges(
  seqnames = ranges$seqnames,
  ranges = IRanges(ranges$start, ranges$end)
  )

phgDsf <-  phgDs |> filterRefRanges(gr)


compare_overlap <- function(phg, sample1, sample2) {
  samp_haps <- phg |> filterSamples(c(sample1, sample2))
  m <- samp_haps |> readHapIds()
  m2 <- m[,!is.na(m[1,]) & !is.na(m[2,]) ]
  same <- m2[,m2[1,] == m2[2,]]
  different <- m2[,m2[1,] != m2[2,]]
  return(ncol(same)/ncol(m2))
  # total_range <- samp_haps |> numberOfRefRanges()
  # divergent_range <- samp_haps |> numberOfHaplotypes(byRefRange=T) |> dplyr::filter(n_haplo > 1)
  # return(1-(nrow(divergent_range)/total_range))
}

exome <- grep("exm", phgDs |> readSamples(), value=T)
gbs <- grep("gbs", phgDs |> readSamples(), value=T)
founders <- grep("exm|gbs", phgDs |> readSamples(), value=T, invert=T)

##compare
identity_table <- data.frame(seq_sample = c(exome, gbs))
for (sample in founders) {
  print(sample)
  comp_vec <- c()
  for (seqsamp in c(exome, gbs)) {
    print(seqsamp)
    comp_vec <- c(comp_vec, compare_overlap(phgDsf, sample, seqsamp))
    print(tail(comp_vec, 1))
  }
  identity_table[[sample]] = round(comp_vec, 2)
}
# write.table(identity_table, file="clipboard", sep="\t", quote=F, row.names=F)
write.table(identity_table, file=glue("{dir}/output/phg_test1_allregions.tsv"), sep="\t", quote=F, row.names=F)


##characterize haploblocks
haploblocks <- phgDs |> numberOfHaplotypes(byRefRange=T)
haploblocks <- haploblocks |>
  mutate(group = sub("chr([0-9]+)[A-Z]$", "\\1", seqnames),
         subgenome = sub("chr[0-9]+([A-Z])$", "\\1", seqnames))
hist(haploblocks$width)
min(haploblocks$width); mean(haploblocks$width); max(haploblocks$width)

##plot haploblock sizes
# ggplot(haploblocks, aes(x=start, y=width)) +
#   geom_point() +
#   facet_grid(rows=vars(group), cols = vars(subgenome), scales="free")
# ggsave(filename=glue("{dir}/figures/haploblock_sizes.png"), width=8, height=12)
# ##plot number of blocks
# ggplot(haploblocks, aes(x=start, y=n_haplo)) +
#   geom_point() +
#   facet_grid(rows=vars(group), cols = vars(subgenome), scales="free")
# 




ref_ranges <- read.delim(glue("{dir}/ref_ranges.bed"), header = F) |>
  dplyr::rename(start = V2, end = V3) |>
  mutate(start = start + 1)

haplo_regionames <- merge(haploblocks, ref_ranges, by = c("start", "end"))


genic_inter_comp <- haplo_regionames[grep("Traes", haplo_regionames$V4, invert=T),] |>
  group_by(seqnames) |>
  summarise(intergenic = cor(width, n_haplo)) |>
  print(n=22)
genic_inter_comp$genic <- haplo_regionames[grep("Traes", haplo_regionames$V4),] |>
  group_by(seqnames) |>
  summarise(genic = cor(width, n_haplo)) |>
  select(genic)


##compare genic by intergenic accuracy
##genic
genic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4),] |>
  dplyr::filter(seqnames != 'chrUnknown')

intergenic_ranges <- haplo_regionames[grep("Traes", haplo_regionames$V4, invert=T),] |>
  dplyr::filter(seqnames != 'chrUnknown') 

genic_r <- GRanges(
  seqnames = genic_ranges$seqnames,
  ranges = IRanges(genic_ranges$start, genic_ranges$end)
)

intergenic_r <- GRanges(
  seqnames = intergenic_ranges$seqnames,
  ranges = IRanges(intergenic_ranges$start, intergenic_ranges$end)
)

phg_genic <-  phgDs |> filterRefRanges(genic_r)
phg_intergenic <-  phgDs |> filterRefRanges(intergenic_r)


##check correlation of genic haploblock width with number of haplotypes at a block
haploblocks <- phg_intergenic |> numberOfHaplotypes(byRefRange=T)
cor(haploblocks$width, haploblocks$n_haplo)
haploblocks |>
  group_by(seqnames) |>
  summarise(cor(width, n_haplo)) |>
  print(n=22)
##check individual chromosomes
hbc <- haploblocks |> filter(seqnames == "chr5D")
plot(hbc$start, hbc$n_haplo)
points(hbc$start, hbc$width/max(hbc$width)*7, col='red')




plot_size=5; plot_width=3; plot_depth=7
res_val<-45
png(filename = glue("{dir}/figures/intergenic_width_by_haplotypes.png"),  
    width=2*plot_width*plot_size*res_val, height=plot_depth*plot_size*res_val)
par(mar=c(4,4,2,3))
plots <- plot_width*plot_depth
nf <- layout(matrix(c(1:plots), plot_depth, plot_width, byrow=T),
             heights=matrix(c(rep(plot_size, plot_width*plot_depth)), plot_depth, byrow=T),
             widths=matrix(c(rep(plot_size*2, plot_width*plot_depth)), plot_depth, byrow=T))

for (chr in paste0(rep(1:7, each=3), rep(c("A", "B", "D"), rep=7))) {
  hbc <- haploblocks |> filter(seqnames == glue("chr{chr}"))
  plot(hbc$start, hbc$n_haplo, main=chr, xlab="Position", ylab="# of haplotypes")
  points(hbc$start, hbc$width/max(hbc$width)*6+1, col='red')
  axis(4, at=1:7, labels= round(seq(1, max(hbc$width), length.out=7)/1e3, 0), col ='red', col.axis='red')
  mtext("Haploblock length (kb)", side = 4, col = "red", line=2)
  abline(v=chrom_pos[chrom_pos$Chromosome == chr, 'Centromere_start'], col='blue', lwd = 3)
  abline(v=chrom_pos[chrom_pos$Chromosome == chr, 'Centromere_end'], col='blue', lwd=3)  
}
dev.off()



identity_table <- data.frame(seq_sample = c(exome))
for (sample in founders) {
  print(sample)
  comp_vec_genic <- c()
  comp_vec_intergenic <- c()
  for (ex in exome) {
    print(ex)
    comp_vec_genic <- c(comp_vec_genic, compare_overlap(phg_genic, sample, ex))
    comp_vec_intergenic <- c(comp_vec_intergenic, compare_overlap(phg_intergenic, sample, ex))
    print(tail(comp_vec_genic, 1))
    print(tail(comp_vec_intergenic, 1))
    
  }
  identity_table[[glue("{sample}_genic")]] = round(comp_vec_genic, 2)
  identity_table[[glue("{sample}_intergenic")]] = round(comp_vec_intergenic, 2)
}
# write.table(identity_table, file=glue("{dir}/output/phg_test1_genic_intergenic.tsv"), sep="\t", quote=F, row.names=F)

identity_table[,grep("inter", colnames(identity_table), invert=T)]


##haplotype assignation
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
length(hapmat_remapped["exmAGS2000_G1", grep("AGS2000", hapmat_remapped["exmAGS2000_G1",])])/ncol(hapmat_remapped)

# id_table_hap <- data.frame()
# 
# for (ex in exome) {
#   print(ex)
#   ext <- paste0(ex, "_G1")
#   ex_v <- hapmat_remapped[paste0(ex, "_G1"),]
#   values <- c()
#   # for (sample in gsub("_G1", "", founders)) {
#   #   values <-c(values, length(ex_v[grep(sample, ex_v)])/ncol(hapmat_remapped))
#   # }
#   values
#   id_table_hap <- rbind(id_table_hap, values)
# }

# colnames(id_table_hap) <- gsub("_G1", "", founders)
# rownames(id_table_hap) <- exome
# round(id_table_hap, 2)
#  
# ##where are AGS2000 off-type calls going?
# ags <- hapmat_remapped[paste0("exmAGS2000", "_G1"),]
# agst <- table(ags)
# agst2 <- round(agst[grep("AGS2000", names(agst), invert=T)]/length(ags),2)
# agst2[agst2 > 0]
# 
# ##where are Hilliard off-type calls going?
# hil <- hapmat_remapped[paste0("exmHILLIARD", "_G1"),]
# hilt <- table(hil)
# hilt2 <- round(hilt[grep("HILLIARD", names(hilt), invert=T)]/length(hil),2)
# hilt2[hilt2 > 0]


##off call list 
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
sink(glue("{dir}/output/haplotype_comparison_breakdown.txt"))
print(phg_stats)
sink()

