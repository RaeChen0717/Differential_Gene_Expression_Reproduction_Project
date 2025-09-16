rm(list=ls())
library(sleuth)
library(biomaRt)
library(devtools)
library('ggplot2')

setwd('/Users/raechen/Downloads/exp15-16/DGE_TranscriptLvl')
base_dir <- '/Users/raechen/Downloads/exp15-16/DGE_TranscriptLvl/Kallisto_out_ZT10'
sample_id <- dir(file.path(base_dir))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, id))
s2c <- read.table("SC_minusATR_vs_plusATR.txt", header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

s2c

mart <- useEnsembl(biomart = "ensembl", dataset = "dmelanogaster_gene_ensembl")
t2g <-biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so<- sleuth_prep(s2c, ~ condition, target_mapping = t2g,aggregation_column = 'ext_gene', gene_mode = TRUE, extra_bootstrap_summary=T,
                 read_bootstrap_tpm = TRUE)
so <- sleuth_fit(so, ~condition, "full")
models(so)

#so <- sleuth_fit(so, ~condition1+condition2, "full")
#so <- sleuth_fit(so, ~condition1, "reduced")

beta='conditionplusATR'
so <- sleuth_wt(so, which_beta = beta)

test_table <- sleuth_results(so, beta)
test_table

test_table$raw_b <- exp(test_table$b)
test_table$log2_b <- log2(test_table$raw_b)
test_table$neg_log10_qval<- -log10(test_table$qval)

sleuth_gene_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
sleuth_gene_matrix

DGE_table<-cbind(row.names(sleuth_gene_matrix), sleuth_gene_matrix)# Add row names as one colum (gene names
p_val_sort<-(test_table[order(test_table$target_id),]) # sort according to gene name
p_val_reduce<-p_val_sort[c("target_id","pval","qval", "log2_b")] # extract only genename,p and q and log2()
DGE_all<-merge(DGE_table, p_val_reduce, by.x=1 , by.y="target_id") # megrge with gene expression and p, q values
DGE_all<-na.omit(DGE_all)
DGE_all<-(DGE_all[order(-DGE_all$log2_b),])
write.csv(DGE_all, 'DGE_table_minusATR_vs_plusATR.csv')

sig_level=0.1
cutoff = 0.58

test_table$diffexpressed <- "Not sig"
test_table$diffexpressed[test_table$log2_b > cutoff & test_table$qval < sig_level] <- "UP"
test_table$diffexpressed[test_table$log2_b < -cutoff & test_table$qval < sig_level] <- "DOWN"

N_significant<-length(test_table$diffexpressed[test_table$diffexpressed !="Not sig"])
N_UP<-length(test_table$diffexpressed[test_table$diffexpressed =="UP"])
N_DOWN<-length(test_table$diffexpressed[test_table$diffexpressed =="DOWN"])

N_significant
N_UP
N_DOWN

library('ggplot2')
ggplot(test_table) + geom_point(aes(x = log2_b, y = neg_log10_qval, col=diffexpressed))+
  geom_vline(xintercept=0, col="black",  linetype="dashed")+
  scale_color_manual(values=c("blue", "black", "red"))+
  labs(x = "log2(FC)", y="-log10(q)", colour="DEG")

sleuth_live(so)

