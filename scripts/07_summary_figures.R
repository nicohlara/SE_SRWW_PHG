
library(glue)
library(dplyr)
library(yarrr)
library(stringr)
library(tidyr)
library(multcompView)


dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"

##summarize identity table
id <- read.table(glue("{dir}/output/haplotype_identity_table.csv"), sep = ",", header=T, row.names = 1)
colnames(id) <- gsub("_genic", "", colnames(id))
##exome identity
exid <- id[grep("exm", rownames(id)),]
rownames(exid) <- gsub("exm|-", "", rownames(exid))
exid <- exid[rownames(exid) %in% colnames(exid), colnames(exid) %in% rownames(exid)]
exid <- diag(as.matrix(exid[colnames(exid), colnames(exid)]))
##GBS identity
gbid <- id[grep("gbs", rownames(id), ignore.case = T),]
rownames(gbid) <- gsub("gbs|_SRGBS|-", "", rownames(gbid))
gbid <- gbid[rownames(gbid) %in% colnames(gbid), colnames(gbid) %in% rownames(gbid)]
gbid <- diag(as.matrix(gbid[colnames(gbid), colnames(gbid)]))
 
id_filt <- data.frame(Sequencing = c(rep("Exome_capture", length(exid)), rep("GBS", length(gbid))),
                      Line = c(names(exid), names(gbid)),
                      Percent_agreement = c(exid, gbid))

id_filt2 <- data.frame(Line = names(exid),
           Exome = exid)
id_filt2 <- merge(id_filt2, data.frame(Line = names(gbid),
                           GBS = gbid), by="Line", all=T)
write.table(id_filt2, "clipboard", quote = F,row.names = F)


pirateplot(Percent_agreement ~ Sequencing, data=id_filt, ylim=c(0,1))
# write.csv(accuracies, glue("{dir}/output/genomic_prediction_accuracies.csv"), row.names=F)

##examine prediction accuracies
# acc <- read.delim(glue("{dir}/output/genomic_prediction_accuracies.csv"), sep=",")
files <- list.files(glue("{dir}/output/accuracy_tests"), full.names=T)
tables = lapply(files, read.delim, sep=",")
accuracies <- do.call('rbind', tables)
accuracies <- dplyr::filter(accuracies, relationship_matrix != "Additive_GRM" & 
                              percent_withheld %in% c(0.25, 0.75))

grm_translate <- data.frame(RM = c("GRM", "HRM", "uSNP_GRM"),
                            rename = c("Biparental SNP", "Haplotype", "Unimputed SNP"))
accuracies$relationship_matrix <- factor(grm_translate$rename[match(accuracies$relationship_matrix, grm_translate$RM)],
          levels=grm_translate$rename)



pvals <- accuracies %>%
  group_by(trait, percent_withheld) %>%
  arrange(relationship_matrix) |>
  group_modify(~ {
    Tk <- TukeyHSD(aov(prediction_accuracy ~ relationship_matrix, data=.x))
    tibble(
      comparison = rownames(Tk$relationship_matrix),
      p_val = Tk$relationship_matrix[,'p adj']
    )
  })

make_cld <- function(df, alpha = 0.05, level_order) {
  comparisons <- df |>
    mutate(sig = p_val < alpha)
  sig_vec <- comparisons$sig
  names(sig_vec) <- comparisons$comparison
  ##insert to rearrange letter order
  raw_letters <- multcompLetters(sig_vec)$Letters
  # letters <- multcompLetters(sig_vec)$Letters
  ordered_groups <- intersect(level_order, names(raw_letters))
  new_letters <- setNames(
    letters[seq_along(ordered_groups)],
    ordered_groups
  )
  letter_map <- setNames(
    new_letters[ordered_groups],
    raw_letters[ordered_groups]
  )
  tibble(
    relationship_matrix = names(raw_letters),
    label = letter_map[raw_letters]
  )
}
cld_df <- pvals |>
  group_by(trait, percent_withheld) |>
  group_modify(~ make_cld(.x,level_order= grm_translate$rename)) |>
  ungroup()

heritability <- read.delim(glue("{dir}/output/blue_heritability.csv"), sep=",")
heritability <- heritability[1:4,]

###plot values
p <- ggplot(data=accuracies, aes(x = relationship_matrix, y = prediction_accuracy)) +
  geom_boxplot() +
  facet_grid(trait ~ percent_withheld) +
  geom_hline(data = heritability,
             aes(yintercept= prediction_limit), inherit.aes=F) +
  scale_y_continuous(limits=c(0.2, 0.9)) 
  # title(main="Comparison of GRM accuracy at 25% or 75% withheld data")
p + geom_text(data = cld_df, aes(x = relationship_matrix, y = 0.7, label = label),
    inherit.aes = FALSE, size = 5, vjust = 0)
ggsave(glue("{dir}/figures/accuracy_assessment.png"), width = 6, height = 8)



# 
# res_val=3
# png(filename = glue('{dir}/figures/accuracy_charts.png'),  
#     width=750*res_val, height=700*res_val, res=72*res_val,
#     bg='transparent')
# # pirateplot(prediction_accuracy ~ relationship_matrix + percent_withheld + trait, data=accuracies, ylim = c(0,1))
# dev.off()


##look at GWAS results
files <- list.files(glue("{dir}/output/GWAS"), full.names=T)
# tables = lapply(files, read.delim, sep=",")
# GWAS <- do.call('rbind', tables)
for (f in files) {
  print(f)
  tb <- read.delim(f, sep=",")
  tb <- tb[tb$p <= 0.1,]
  print(dim(tb))
  if (grepl('haplotype', f)) {
    tb <- tb %>%
      mutate(haploregion = id) |>
      separate(haploregion, into = c("chr", "range"), sep = ":") %>%
      separate(range, into = c("pos", "end"), sep = "-") %>%
      mutate(chr = sub("^chr", "", chr))
  } else {
    tb$chr <- str_replace(str_replace(tb$id, "^S", ""), "_\\d*$", "")
    tb$pos <- as.numeric(str_replace(tb$id, "^S\\d[ABD]_", ""))
  }
  pname <- tools::file_path_sans_ext(basename(f))
  ggplot(tb, aes(x=pos, y=-log10(p))) +
    geom_point() +
    facet_grid(cols=vars(chr))
  ggsave(glue("{dir}/figures/GWAS/{pname}.png", length=9, height=5))
}



