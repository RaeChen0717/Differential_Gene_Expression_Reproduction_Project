rm(list=ls())
library(edgeR)
library(limma)

setwd('/Users/raechen/Downloads/exp15-16/DGE_THIP')
featurecnts <- read.delim("counts/gene_counts.txt", comment.char="#", check.names = FALSE)

count_cols <- grep("\\.sorted\\.bam$", colnames(featurecnts))
cnt <- as.matrix(featurecnts[, 7:ncol(featurecnts)])
rownames(cnt) <- featurecnts$Geneid
colnames(cnt) <- sub("^.*/|\\.sorted\\.bam$", "", colnames(cnt))

samples <- colnames(cnt)
group <- ifelse(grepl("^experiment15_", samples), "plusATR", "minusATR")
group <- factor(group, levels = c("minusATR", "plusATR"))

y <- DGEList(counts = cnt, group = group)

#keep <- filterByExpr(y, group = group)

cpm_mat <- cpm(y)

mean_cpm_minus <- rowMeans(cpm_mat[, group == "minusATR", drop = FALSE])
mean_cpm_plus  <- rowMeans(cpm_mat[, group == "plusATR",  drop = FALSE])

keep <- (mean_cpm_minus >= 5) | (mean_cpm_plus >= 5)
y <- y[keep,, keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method="TMM")

design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)

tab <- topTags(qlf, n = Inf)$table

fc_cutoff <- 0.58    # ~1.5x
sig_level <- 0.05

tab$flag <- "NotSig"
tab$flag[ tab$FDR < sig_level & tab$logFC >  fc_cutoff] <- "UP"
tab$flag[ tab$FDR < sig_level & tab$logFC < -fc_cutoff] <- "DOWN"

cat("Summary (FDR<0.05 & |log2FC|>0.58):\n")
print(table(tab$flag))

write.csv(tab, "result/edgeR_exp15_vs_exp16_results.csv", row.names = TRUE)
