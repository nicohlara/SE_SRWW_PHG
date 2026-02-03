library(sommer)
library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)
library(rrBLUP)

haplotypes <- phg@hapIds
# haplotypes <- haplotypes[grep("Ref", rownames(haplotypes), invert=T),]
# haplotypes <- haplotypes[,apply(haplotypes, 2, function(x) !all(is.na(x)))]

SNPs <- as.matrix(genotype)
SNPs_subset <- genotype@snps[c("chr", "pos", "id")]
markers <- c()
for (ch in unique(SNPs_subset$chr)) {
  a <- SNPs_subset[SNPs_subset$chr == ch,]
  markers <- c(markers, sample(a$id, 50))
}
SNPs <- SNPs[,markers]

bl <- filter(blues, Entry %in% row.names(SNPs))
rownames(bl) <- bl$Entry
bl <- bl[row.names(SNPs),]
SNPs <- SNPs[row.names(bl),]
# bl <- as.matrix(bl)

# GRM <- A.mat(SNPs)
# 
# null <- mmer(
#   PM ~ 1,
#   random = ~ vs(Entry, Gu = GRM),
#   data = bl
# )
# Vinv <- solve(null$Vi)
# C <- matrix(1, nrow = length(y), ncol = 1)
# CtVinv <- crossprod(C, Vinv)
# CtVinvC <- CtVinv %*% C
# P <- Vinv - Vinv %*% C %*% solve(CtVinvC) %*% CtVinv


# assoc_test <- function(y, X, Vinv) {
#   # X: n x p design matrix
#   XtVinv <- crossprod(X, Vinv)
#   XtVinvX <- XtVinv %*% X
#   XtVinvY <- XtVinv %*% y
#   
#   beta_hat <- solve(XtVinvX, XtVinvY)
#   stat <- drop(t(beta_hat) %*% XtVinvX %*% beta_hat)
#   
#   pchisq(stat, df = ncol(X), lower.tail = FALSE)
# }
assoc_test <- function(y, X, P) {
  XtP <- crossprod(X, P)
  XtPX <- XtP %*% X
  XtPy <- XtP %*% y
  
  beta <- solve(XtPX, XtPy)
  stat <- drop(t(beta) %*% XtPX %*% beta)
  
  pchisq(stat, df = ncol(X), lower.tail = FALSE)
}

# 
# mh <- data.frame()
# for (m in markers) {
#   SNPset <- SNPs[,m]
#   p <- assoc_test(bl[,"PM"], matrix(SNPs[,m], ncol=1) , P)
#   mh <- rbind(mh, data.frame(marker = m, pval = p))
# }
# 
# mh[,c("chr", "pos")] <- str_split_fixed(mh$marker, "_", 2)
# mh <- mh %>%
#   mutate(pos = as.numeric(pos))
# ggplot(mh, aes(x = pos, y = -log10(pval))) +
#   geom_point() +
#   facet_grid(cols = vars(chr))
# 
# 
# 
# ##sommer GWAS
# 
# mix2 <- GWAS(color~1,
#              random=~vsm(ism(id), Gu=A) + Rowf + Colf,
#              rcov=~units, M=GT, gTerm = "u:id",
#              verbose = FALSE, nIters=10,
#              data=DT)
# 
# 
# 
# ##slow dumb GWAS loop to compare
# len <- ncol(SNPs)
# count = 0
# print(Sys.time())
# plot_df <- data.frame()
# j <- 5
# for (i in c(1:len)) {
#   if (round((i/len)*100) > count) {
#     count <- round((i/len)*100)
#     print(paste(count, "%", sep=""))
#   }
#   b <- data.frame(Entry = rownames(SNPs), Marker = SNPs[,i])
#   c <- colnames(SNPs)[i]
#   if (length(unique(b$Marker)) > 1) {
#     d <- merge(blues,b, by="Entry")
#     p_vals <- data.frame(id = c)
#     mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Cross_ID)))
#     p <- 2 * (1 - pnorm(abs(data.frame(coef(summary(mm)))['Marker','t.value'])))
#     p_vals[paste("p_", names(blues[j]), sep="")] <- p
#         mm <- suppressMessages(lmer(data = d, d[,j] ~ Marker + (1|Cross_ID)))
#     }
#     plot_df <- rbind(plot_df, p_vals)
# }
# 
# 
# plot_df$chr <- str_replace(str_replace(plot_df$id, "^S", ""), "_\\d*$", "")
# plot_df$pos <- as.numeric(str_replace(plot_df$id, "^S\\d[ABD]_", ""))
# ggplot(plot_df, aes(x = pos, y = -log10(p_PM))) +
#   geom_point() +
#   facet_grid(cols = vars(chr))


###reimplementation of basic model

y <- bl[["PM"]]
##fit the family effect. I suspect I could replace this with a GRM?
# null <- lmer(y ~ 1 + (1 | Cross_ID), data = bl, REML = TRUE)
null <- mmer( PM ~ 1,random = ~ vs(Entry, Gu = GRM),
              data = bl)
# null <- mixed.solve(y=y, K=GRM, SE=F)
# y_resid <- resid(null)
# y_resid <- y - c(null$beta) - null$u
y_resid <- y - c(null$Beta$Estimate) - null$U$`u:Entry`$PM
G <- as.matrix(SNPs)
fast_assoc <- function(y, G) {
  pvals <- numeric(ncol(G))
  for (i in seq_len(ncol(G))) {
    g <- G[, i]
    if (length(unique(g)) < 2) {
      pvals[i] <- NA
      next
    }
    fit <- lm(y ~ g)
    pvals[i] <- summary(fit)$coefficients[2, 4]
  }
  pvals
}
p <- fast_assoc(y_resid, G)
plot_df <- data.frame(id = colnames(G), p = p)
plot_df$chr <- str_replace(str_replace(plot_df$id, "^S", ""), "_\\d*$", "")
plot_df$pos <- as.numeric(str_replace(plot_df$id, "^S\\d[ABD]_", ""))
ggplot(plot_df, aes(x = pos, y = -log10(p))) +
  geom_point() +
  facet_grid(cols = vars(chr))

