
library(glue)
library(dplyr)
library(yarrr)

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
accuracies <- dplyr::filter(accuracies, relationship_matrix != "Additive_GRM")

pvals <- accuracies %>%
  group_by(trait, percent_withheld) %>%
  group_modify(~ {
    tt <- t.test(prediction_accuracy ~ relationship_matrix, data = .x)
    tibble(
      p_value = tt$p.value,
      t_stat  = tt$statistic,
      estimate_diff = diff(tt$estimate)
    )
  }) %>%
  ungroup() |>
  mutate(sig = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

res_val=
png(filename = glue('{dir}/figures/accuracy_charts.png'),  
    width=750*res_val, height=700*res_val, res=72*res_val,
    bg='transparent')
pirateplot(prediction_accuracy ~ relationship_matrix + percent_withheld + trait, data=accuracies, ylim = c(0,1))
dev.off()
