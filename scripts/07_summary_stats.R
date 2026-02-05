

library(vegan)


###test similarity of Haplotype and SNP-based GRM
GRM <- read.delim(glue("{PHG_dir}/GRM.csv"), sep=",", row.names = 1)
colnames(GRM) <- rownames(GRM)
HRM <- read.delim(glue("{PHG_dir}/HRM.csv"), sep=",", row.names = 1)
colnames(HRM) <- rownames(HRM)
print(mantel(GRM, HRM))


