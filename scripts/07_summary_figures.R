
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
acc <- read.delim(glue("{dir}/output/genomic_prediction_accuracies.csv"), sep=",")
pirateplot(prediction_accuracy ~ relationship_matrix + percent_witheld + trait, data=acc, ylim = c(0,1))

