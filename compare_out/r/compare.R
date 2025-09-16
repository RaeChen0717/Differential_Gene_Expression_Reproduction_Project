# compare_all.R
# ------------------------------------------------------------------
# Purpose: Compare the following DGE
#   1) Tutorial  vs Published
#   2) Tutorial  vs edgeR
#   3) edgeR     vs Published
# Inputs:
#   - Tutorial CSV:    DGE_table_minusATR_vs_plusATR.csv
#   - edgeR CSV: edgeR_exp15_vs_exp16_results.csv
#   - Published XLSX:  elife-88198-fig7-data1-v1.xlsx
#   - Background set:  background_fbgn.txt
# Outputs (./compare_out):
#   - fisher_results.csv, matched_pairs.csv
#   - venn_up.png, venn_down.png, scatter_logFC.png, volcano_markOverlap.png
# ------------------------------------------------------------------

rm(list = ls())
suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(stringr)
  library(ggplot2); library(biomaRt); library(VennDiagram); library(grid)
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagramLogger")
    futile.logger::flog.appender(function(...) invisible(NULL), name = "VennDiagramLogger")
  }
})

# ========== CONFIG ==========
setwd('/Users/raechen/Downloads/exp15-16')

tutorial_csv   <- "DGE_Tutorial/DGE_table_minusATR_vs_plusATR.csv"     
edgeR_csv      <- "DGE_THIP/edgeR_exp15_vs_exp16_results.csv"          
paper_xlsx     <- "elife-88198-fig7-data1-v1.xlsx"                  

# Threshold：FDR<0.05 & |log2FC|>0.58
qval_cut  <- 0.05
logfc_cut <- 0.58

# outdir
outdir <- "compare_out"

# ========== HELPERS ==========

# 1. (if not FBgn) gene symbol → FBgn
map_to_fbgn <- local({
  mart <- NULL   # store the BioMart connection later
  function(x) {  # input vector x: gene identifiers
    x <- as.character(x)  # Ensure x is character
    
    # connect to Ensembl BioMart (Drosophila dataset)
    if (is.null(mart)) {
      mart <<- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
    }
    
    # detect which are already FlyBase IDs
    is_fbgn <- str_detect(x, "^FBgn\\d+")
    
    # initialize output as NA, then copy over values that are already FBgn IDs
    out <- rep(NA_character_, length(x))
    out[is_fbgn] <- x[is_fbgn]
    
    # collect unique gene symbols (non-FBgn inputs) to be mapped
    need <- unique(na.omit(x[!is_fbgn]))
    
    # if there are any symbols that need conversion
    if (length(need)) {
      # Query BioMart for mapping between gene symbol and FlyBase ID
      mp <- getBM(attributes = c("external_gene_name","flybase_gene_id"),
                  filters = "external_gene_name", values = need, mart = mart)
      
      # lookup table: symbol -> FBgn ID
      lut <- setNames(mp$flybase_gene_id, mp$external_gene_name)
      
      # fill in the missing entries
      out[!is_fbgn] <- lut[x[!is_fbgn]]
    }
    
    out
  }
})


# 2. read inputs: FBgn, logFC, qval
read_tutorial <- function(path) {
  df <- suppressMessages(read_csv(path, show_col_types = FALSE))
  tibble(
    FBgn  = map_to_fbgn(df$V1),
    logFC = as.numeric(df$log2_b),
    qval  = as.numeric(df$qval)
  ) %>% filter(!is.na(FBgn))
}

read_edgeR <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(FBgn = 1) %>%           
    mutate(FBgn = as.character(FBgn),
           logFC = as.numeric(logFC),
           qval  = as.numeric(FDR)) %>%   
    select(FBgn, logFC, qval) %>%
    filter(!is.na(FBgn))
}

read_published <- function(path) {
    read_xlsx(path, sheet = 2) %>%
      transmute(FBgn  = map_to_fbgn(Row.names),
                logFC = as.numeric(logFC),
                qval  = as.numeric(FDR)) %>%
      filter(!is.na(FBgn))
}

# 3. split UP/DOWN 
split_up_down <- function(df, logfc_cut, qval_cut) {
  sig <- df %>% filter(qval < qval_cut)
  list(
    up   = sig %>% filter(logFC >  logfc_cut) %>% pull(FBgn) %>% unique(), #filtered the repeats (different transcripts)
    down = sig %>% filter(logFC < -logfc_cut) %>% pull(FBgn) %>% unique()
  )
}

# 4. Venn（raw + percent；sigdigs = 3）
save_venn <- function(A, B, labels=c("A","B"), path="venn.png") {
  gr <- VennDiagram::venn.diagram(
    x = list(A=A, B=B), category = labels,
    filename = NULL, imagetype = "png",
    fill = c("#CFE8F3","#FDE0DD"), alpha = 0.7,
    cat.cex = 1.2, cex = 1.4, lwd = 2, quiet = TRUE,
    print.mode = c("raw","percent"), sigdigs = 0
  )
  png(path, width=1600, height=1200, res=200)
  grid.newpage(); grid.draw(gr); dev.off()
}

# 5. Fisher 2×2（CSV，OddsRatio，qFDR）
fisher_csv <- function(A, B, background,
                       direction = "UP",
                       outdir = ".",
                       A_name = "A",
                       B_name = "B") {
  A <- intersect(A, background); B <- intersect(B, background)
  a <- length(intersect(A,B)); b <- length(setdiff(A,B))
  c <- length(setdiff(B,A));  d <- length(setdiff(background, union(A,B)))
  ft <- fisher.test(matrix(c(a,b,c,d), 2), alternative = "greater")
  
  df <- data.frame(
    row_label                = c(paste("In",  A_name), paste("Not in", A_name)),
    `In B`                   = c(a, c),
    `Not in B`               = c(b, d),
    check.names = FALSE
  )
  names(df)[2:3] <- c(paste("In",  B_name), paste("Not in", B_name))
  
  stats <- data.frame(
    row_label                = c("odds_ratio", "p_value (BH adj)"),
    `In B`                   = c(sprintf("%.3f", unname(ft$estimate)),
                                 sprintf("%.3g", p.adjust(ft$p.value, method="BH"))),
    `Not in B`               = c("", ""),
    check.names = FALSE
  )
  names(stats)[2:3] <- c(paste("In",  B_name), paste("Not in", B_name))
  
  out <- rbind(df, stats)
  fp <- file.path(outdir, sprintf("fisher_table_%s.csv", tolower(direction)))
  write.csv(out, fp, row.names = FALSE); message("Saved: ", fp)
  invisible(list(a=a,b=b,c=c,d=d, table=out))
}

# 6. Jaccard / only 
jaccard <- function(A,B) { u <- union(A,B); if (length(u)==0) NA_real_ else length(intersect(A,B))/length(u) }

dump_overlap_lists <- function(A, B, tag, outdir) {
  writeLines(sort(intersect(A,B)), file.path(outdir, paste0("intersect_", tag, ".txt")))
  writeLines(sort(setdiff(A,B)),   file.path(outdir, paste0("only_A_",   tag, ".txt")))
  writeLines(sort(setdiff(B,A)),   file.path(outdir, paste0("only_B_",   tag, ".txt")))
}

# 7. logFC correlation：scatter plot + regression + Pearson/Spearman
plot_logFC_correlation <- function(dfA, dfB, labels=c("A","B"), outdir=".", prefix="A_vs_B") {
  matched <- inner_join(
    dfA %>% select(FBgn, logFC_A = logFC),
    dfB %>% select(FBgn, logFC_B = logFC),
    by = "FBgn"
  )
  if (!nrow(matched)) stop("No overlapping genes for correlation plot.")
  r_pear  <- suppressWarnings(cor(matched$logFC_A, matched$logFC_B, method="pearson"))
  r_spear <- suppressWarnings(cor(matched$logFC_A, matched$logFC_B, method="spearman"))
  pt <- cor.test(matched$logFC_A, matched$logFC_B, method="pearson")
  st <- cor.test(matched$logFC_A, matched$logFC_B, method="spearman")
  sink(file.path(outdir, paste0(prefix, "_correlation_stats.txt")))
  cat("### Correlation between ", labels[1], " and ", labels[2], " logFC ###\n", sep="")
  cat("\n-- Pearson --\n");  print(pt)
  cat("\n-- Spearman --\n"); print(st)
  sink()
  p <- ggplot(matched, aes(logFC_A, logFC_B)) +
    geom_point(alpha=.6) +
    geom_smooth(method="lm", se=FALSE, color="red") +
    labs(x=paste0(labels[1], " log2FC"), y=paste0(labels[2], " log2FC"),
         title=sprintf("Pearson r=%.2f, Spearman ρ=%.2f (n=%d)", r_pear, r_spear, nrow(matched))) +
    theme_bw()
  ggsave(file.path(outdir, paste0(prefix, "_logFC_correlation.png")), p, width=7, height=5, dpi=200)
  invisible(list(n=nrow(matched), pearson=r_pear, spearman=r_spear))
}

# 8. Volcano Plot
volcano_mark_overlap <- function(dfA, setsB, qcut=0.05, lfc=0.58, labelA="A", labelB="B", outpng="volcano.png") {
  dfA$in_B <- dfA$FBgn %in% union(setsB$up, setsB$down)
  volcano <- dfA %>% mutate(sig = ifelse(qval < qcut & abs(logFC) > lfc, "sig", "ns"))
  p <- ggplot(volcano, aes(logFC, -log10(pmax(qval, 1e-300)), alpha=sig)) +
    geom_point(aes(shape=in_B)) +
    scale_alpha_manual(values=c(ns=0.3, sig=0.9)) +
    labs(x=paste0("log2FC (", labelA, ")"),
         y="-log10(FDR)", shape=paste0("In ", labelB, " DEG"),
         title=paste0(labelA, " volcano (", labelB, " overlap highlighted)")) +
    theme_bw()
  ggsave(outpng, p, width=7, height=5, dpi=200)
}

# ========== read inputs ==========
tut <- read_tutorial(tutorial_csv)
edg <- read_edgeR(edgeR_csv)
pub <- read_published(paper_xlsx)

# background set: union of all FBgn lists
background_all <- unique(c(tut$FBgn, edg$FBgn, pub$FBgn))

# ========== comparison ==========
comparisons <- list(
  list(A = tut, A_name="Tutorial", B = pub, B_name="Published", out = file.path(outdir, "Tutorial_vs_Published")),
  list(A = tut, A_name="Tutorial", B = edg, B_name="edgeR",     out = file.path(outdir, "Tutorial_vs_edgeR")),
  list(A = edg, A_name="edgeR",     B = pub, B_name="Published", out = file.path(outdir, "edgeR_vs_Published"))
)

for (cmp in comparisons) {
  A <- cmp$A; B <- cmp$B; A_name <- cmp$A_name; B_name <- cmp$B_name; outdir <- cmp$out
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  message("=== Running: ", A_name, " vs ", B_name, " ===")
  
  # 1) split based on threshold
  setsA <- split_up_down(A, logfc_cut, qval_cut)
  setsB <- split_up_down(B, logfc_cut, qval_cut)
  
  # 2) Venn（UP/DOWN)
  save_venn(setsA$up,   setsB$up,   c(paste(A_name,"UP"),   paste(B_name,"UP")),
            file.path(outdir, "venn_up.png"))
  save_venn(setsA$down, setsB$down, c(paste(A_name,"DOWN"), paste(B_name,"DOWN")),
            file.path(outdir, "venn_down.png"))
  
  # 3) Fisher 2×2（UP/DOWN.CSV）
  fisher_csv(setsA$up,   setsB$up,   background_all, direction="UP",
             outdir=outdir, A_name=A_name, B_name=B_name)
  fisher_csv(setsA$down, setsB$down, background_all, direction="DOWN",
             outdir=outdir, A_name=A_name, B_name=B_name)
  
  # 4) Jaccard（CSV）
  j_up   <- jaccard(setsA$up,   setsB$up)
  j_down <- jaccard(setsA$down, setsB$down)
  size_tbl <- tibble(
    set = c(paste0(A_name,"_UP"), paste0(B_name,"_UP"), "Overlap_UP",
            paste0(A_name,"_DOWN"), paste0(B_name,"_DOWN"), "Overlap_DOWN"),
    size = c(length(setsA$up), length(setsB$up), length(intersect(setsA$up, setsB$up)),
             length(setsA$down), length(setsB$down), length(intersect(setsA$down, setsB$down)))
  )
  write_csv(size_tbl, file.path(outdir, "overlap_sizes.csv"))
  write_csv(tibble(jaccard_up=j_up, jaccard_down=j_down), file.path(outdir, "jaccard.csv"))
  
  # 5) only
  dump_overlap_lists(setsA$up,   setsB$up,   "UP",   outdir)
  dump_overlap_lists(setsA$down, setsB$down, "DOWN", outdir)
  
  # 6) logFC
  plot_logFC_correlation(A, B, labels=c(A_name, B_name),
                         outdir=outdir, prefix=paste0(A_name, "_vs_", B_name))
  
  # 7) volcano plot
  volcano_mark_overlap(A, setsB, qcut=qval_cut, lfc=logfc_cut,
                       labelA=A_name, labelB=B_name,
                       outpng=file.path(outdir, paste0("volcano_", A_name, "_mark_", B_name, ".png")))
}

message("All done. See results under: ", normalizePath(outdir))
