
library(sommer)
library(glue)
library(dplyr)

dir="c:/Users/nalara/Documents/GitHub/SE_SRWW_PHG"

calculate_prediction_accuracy <- function(phenotype_df, relationship_matrix, seed = 12345, withold_percent = 0.2) {
  set.seed(seed)
  y.trn <- phenotype_df
  vv <- sample(rownames(phenotype_df),round(nrow(phenotype_df)*withold_percent))
  y.trn[vv,"pheno"] <- NA
  ##run prediction model
  ans <- mmes(pheno~1,
              random=~vsm(ism(Entry),Gu=relationship_matrix), 
              rcov=~units,nIters=10,
              data=y.trn, verbose = FALSE) # kinship based
  return(cor(ans$u[vv,] ,phenotype_df[vv,"pheno"], use="complete"))
}


seeds <- c(45174796, 69322386, 18488794, 99566509, 17074345, 94020839, 80272836, 86566472, 82166244, 93897769)
families <- c("UX1989", "UX1991", "UX1992", "UX1993", "UX1994",
              "UX1995", "UX1997", "UX2000", "UX2010", "UX2012", 
              "UX2013", "UX2023", "UX2026", "UX2029", "UX2031")
# seeds <- c(45174796)

accuracies <- data.frame()
drop_one_family <- data.frame()
for (mat in c("GRM", "Additive_GRM", "HRM")) {
  ##read in and format Kinship matrix
  K <- read.delim(glue("{dir}/output/{mat}.csv"), sep=",", row.names = 1)
  colnames(K) <- rownames(K)
  K <- as.matrix(K)
  ###read in phenotypic values
  blues <- read.delim(glue("{dir}/data/blues.csv"), sep=",")
  row.names(blues) <- blues$Entry
  ##line up datasets
  blues <- blues[blues$Entry %in% row.names(K),]
  K <- K[blues$Entry, blues$Entry]
  ##perform traitwise GP
  for (trait in c("WDR", "Height", "HD", "PM")) {
    phenotype_df <- blues |> dplyr::select(Entry, !!as.symbol(trait))
    colnames(phenotype_df) <- c("Entry", "pheno")
    for (wp in c(0.1, 0.25, 0.5, 0.75)) {
      print(glue("Predicting {trait} with {wp} portion witheld data"))
      for (se in seeds) {
        pred_acc <- calculate_prediction_accuracy(phenotype_df, K, seed= se, withold_percent=wp)
        # print(pred_acc)
        df <- data.frame(trait = trait, percent_withheld = wp, prediction_accuracy = pred_acc)
        df$relationship_matrix <- mat
        accuracies <- rbind(accuracies, df)
      }
    }
    for (fam in families) {
      y.trn <- phenotype_df
      vv <- grep(fam, rownames(phenotype_df), value=T)
      y.trn[vv,"pheno"] <- NA
      ##run prediction model
      ans <- mmes(pheno~1,
                  random=~vsm(ism(Entry),Gu=K), 
                  rcov=~units,nIters=10,
                  data=y.trn, verbose = FALSE) # kinship based
      pred_acc <- cor(ans$u[vv,] ,phenotype_df[vv,"pheno"], use="complete")
      df <- data.frame(trait = trait, 
                       family_withhheld = fam, 
                       prediction_accuracy = pred_acc)
      drop_one_family <- rbind(drop_one_family, df)
    }
  }
}

write.csv(accuracies, glue("{dir}/output/genomic_prediction_accuracies.csv"), row.names=F)
write.csv(drop_one_family, glue("{dir}/output/genomic_prediction_accuracies_familywise.csv"), row.names=F)

# pirateplot(prediction_accuracy ~ percent_witheld + trait, data=accuracies)

