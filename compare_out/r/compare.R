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
  library(readr); library(dplyr); library(stringr)
  library(ggplot2); library(biomaRt); library(VennDiagram); library(grid); library(openxlsx) 
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(futile.logger::FATAL, name = "VennDiagramLogger")
    futile.logger::flog.appender(function(...) invisible(NULL), name = "VennDiagramLogger")
  }
})

# ========== CONFIG ==========
setwd('/Users/raechen/Downloads/exp15-16')

tutorial_csv   <- "DGE_Tutorial/result/DGE_table_minusATR_vs_plusATR.csv"     
edgeR_csv      <- "DGE_THIP/result/edgeR_exp15_vs_exp16_results.csv"          
paper_csv      <- "published_data.csv"                  

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
    FBgn  = as.character(df$flybase_gene_id),
    logFC = as.numeric(df$log2_b),
    qval  = as.numeric(df$qval)
  ) %>% filter(!is.na(FBgn), FBgn != "-")
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
  read_csv(path, show_col_types = FALSE) %>%
    transmute(FBgn  = as.character(flybase_gene_id),
              logFC = as.numeric(logFC),
              qval  = as.numeric(FDR)) %>%
    filter(!is.na(FBgn), FBgn != "-")
}

# 3. split UP/DOWN 
split_up_down <- function(df, logfc_cut, qval_cut) {
  sig <- df %>% filter(qval < qval_cut)
  list(
    up   = sig %>% filter(logFC >  logfc_cut) %>% pull(FBgn) %>% unique(), #filtered the repeats (different transcripts)
    down = sig %>% filter(logFC < -logfc_cut) %>% pull(FBgn) %>% unique()
  )
}

# 4. Venn（raw + percent；sigdigs = 0）
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
                       A_name = "A",
                       B_name = "B") {
  A <- intersect(A, background); B <- intersect(B, background)
  a <- length(intersect(A,B)); b <- length(setdiff(A,B))
  c <- length(setdiff(B,A));  d <- length(setdiff(background, union(A,B)))
  ft <- fisher.test(matrix(c(a,b,c,d), 2), alternative = "greater")
  
  df <- data.frame(
    row_label = c(paste("In",  A_name), paste("Not in", A_name)),
    `In B`    = c(a, c),
    `Not in B`= c(b, d),
    check.names = FALSE
  )
  names(df)[2:3] <- c(paste("In",  B_name), paste("Not in", B_name))
  
  stats <- data.frame(
    row_label = c("odds_ratio", "p_value (BH adj)"),
    `In B`    = c(sprintf("%.3f", unname(ft$estimate)),
                  sprintf("%.3g", p.adjust(ft$p.value, method="BH"))),
    `Not in B`= c("", ""),
    check.names = FALSE
  )
  names(stats)[2:3] <- c(paste("In",  B_name), paste("Not in", B_name))
  
  out <- rbind(df, stats)
  return(list(a=a,b=b,c=c,d=d, table=out))
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
volcano_mark_overlap <- function(dfA, B_or_setsB, qcut=0.05, lfc=0.58,
                                 labelA="A", labelB="B", outpng="volcano.png") {
  # B data.frame/setsB（list: up/down
  get_sigB <- function(dfA_FBgn, B_or_setsB) {
    # list(up=,down=)
    if (is.list(B_or_setsB) && all(c("up","down") %in% names(B_or_setsB))) {
      return(dfA_FBgn %in% union(B_or_setsB$up, B_or_setsB$down))
    }
    # data.frame （make sure data type
    B <- B_or_setsB
    if (!("FBgn" %in% names(B))) stop("dfB/setsB needs FBgn column or up/down lists.")
    B <- B %>%
      mutate(
        logFC = suppressWarnings(as.numeric(logFC)),
        qval  = suppressWarnings(as.numeric(qval))
      )
    sigB_ids <- B %>%
      filter(!is.na(qval), !is.na(logFC)) %>%
      filter(qval < qcut & abs(logFC) > lfc) %>%
      pull(FBgn) %>% unique()
    dfA_FBgn %in% sigB_ids
  }
  
  dfA <- dfA %>%
    mutate(
      logFC = suppressWarnings(as.numeric(logFC)),
      qval  = suppressWarnings(as.numeric(qval))
    )
  
  is_sigA <- (dfA$qval < qcut) & (abs(dfA$logFC) > lfc)
  sigB_for_A <- get_sigB(dfA$FBgn, B_or_setsB)
  
  cls <- ifelse(!is_sigA, "ns",
                ifelse(sigB_for_A, "sig_both", "sig_only"))
  
  volcano <- dfA %>%
    mutate(cls = factor(cls, levels = c("ns","sig_only","sig_both")))
  
  n_sig    <- sum(volcano$cls != "ns", na.rm = TRUE)
  n_both   <- sum(volcano$cls == "sig_both", na.rm = TRUE)
  pct_both <- if (n_sig > 0) round(100 * n_both / n_sig, 1) else 0
  
  p <- ggplot(volcano, aes(logFC, -log10(pmax(qval, 1e-300)))) +
    geom_hline(yintercept = -log10(qcut), linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(color = cls), alpha = 0.75, size = 1.2) +
    scale_color_manual(
      name = paste0("DEG categories (", labelA, " vs ", labelB, ")"),
      values = c(ns="grey80", sig_only="#D62728", sig_both="#2CA02C"),
      labels = c(
        ns        = paste0("Not significant in ", labelA),
        sig_only  = paste0("Significant only in ", labelA, " or ", labelB),
        sig_both  = paste0("Significant in both ", labelA, " and ", labelB)
      )
    ) +
    labs(
      title    = paste0("Volcano plot of ", labelA, " (", labelB, " overlap highlighted)"),
      subtitle = sprintf("Significant in %s: %d | Overlap with %s: %d (%.1f%% of sig)",
                         labelA, n_sig, labelB, n_both, pct_both),
      x = paste0("log2FC (", labelA, ")"),
      y = "-log10(FDR)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right")
  
  ggsave(outpng, p, width = 7, height = 5, dpi = 200)
}

# 9. output genes with discordant direction
pairwise_table_disagree <- function(dfA, dfB,
                                    A_name = "A", B_name = "B",
                                    qcut = 0.05, lfc = 0.58,
                                    outdir = ".",
                                    csv_prefix = "pairs_disagree") {
  tbl <- dplyr::inner_join(
    dfA %>% dplyr::select(FBgn,
                          logFC_A = logFC,
                          q_A     = qval),
    dfB %>% dplyr::select(FBgn,
                          logFC_B = logFC,
                          q_B     = qval),
    by = "FBgn"
  ) %>%
    dplyr::mutate(
      direction_agree = sign(logFC_A) == sign(logFC_B),
      sig_A = (q_A < qcut) & (abs(logFC_A) > lfc),
      sig_B = (q_B < qcut) & (abs(logFC_B) > lfc),
      category = dplyr::case_when(
        sig_A & sig_B ~ "both_sig",
        sig_A & !sig_B ~ paste0("sig_only_", A_name),
        !sig_A & sig_B ~ paste0("sig_only_", B_name),
        TRUE ~ "ns_both"
      ),
      delta_logFC = logFC_B - logFC_A,
      abs_delta   = abs(delta_logFC)
    ) %>%
    dplyr::rename(
      !!paste0("logFC_", A_name) := logFC_A,
      !!paste0("logFC_", B_name) := logFC_B,
      !!paste0("qval_",  A_name) := q_A,
      !!paste0("qval_",  B_name) := q_B
    ) %>%
    dplyr::filter(!direction_agree) %>%    
    dplyr::arrange(dplyr::desc(abs_delta))
    fp <- file.path(outdir, sprintf("%s_%s_vs_%s.csv", csv_prefix, A_name, B_name))
    readr::write_csv(tbl, fp)
    message("Saved disagree-only table -> ", fp)
  
    return(tbl)
}

# 10.  histogram & CDF
# dfA/dfB: FBgn, logFC, qval
plot_sigonly_hist_and_cdf <- function(dfA, dfB,
                                      labelA="A", labelB="B",
                                      qcut=0.05, lfc=0.58,
                                      outdir=".") {
  stopifnot(all(c("FBgn","logFC","qval") %in% names(dfA)),
            all(c("FBgn","logFC","qval") %in% names(dfB)))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  sigA_ids <- dfA %>% 
    mutate(qval = as.numeric(qval), logFC = as.numeric(logFC)) %>%
    filter(qval < qcut, abs(logFC) > lfc) %>% pull(FBgn) %>% unique()
  sigB_ids <- dfB %>% 
    mutate(qval = as.numeric(qval), logFC = as.numeric(logFC)) %>%
    filter(qval < qcut, abs(logFC) > lfc) %>% pull(FBgn) %>% unique()
  
  A_only_ids <- setdiff(sigA_ids, sigB_ids)
  B_only_ids <- setdiff(sigB_ids, sigA_ids)
  
  message(sprintf("[%s vs %s] sigA=%d, sigB=%d, A-only=%d, B-only=%d, overlap=%d",
                  labelA, labelB, length(sigA_ids), length(sigB_ids),
                  length(A_only_ids), length(B_only_ids),
                  length(intersect(sigA_ids, sigB_ids))))
  
  q_in_B_for_Aonly <- dfB %>% 
    select(FBgn, qB = qval) %>% mutate(qB = as.numeric(qB)) %>%
    filter(FBgn %in% A_only_ids, is.finite(qB)) %>% pull(qB)
  
  q_in_A_for_Bonly <- dfA %>% 
    select(FBgn, qA = qval) %>% mutate(qA = as.numeric(qA)) %>%
    filter(FBgn %in% B_only_ids, is.finite(qA)) %>% pull(qA)
  
  .safe_nlog10 <- function(p) -log10(pmax(as.numeric(p), 1e-300))
  
  plot_hist <- function(vals, title, xlab, out_png) {
    df <- data.frame(val = vals)
    p <- ggplot(df, aes(val)) +
      geom_histogram(bins = 40, alpha = 0.9) +
      geom_vline(xintercept = -log10(qcut), linetype = "dashed") +
      labs(title = title,
           subtitle = sprintf("n=%d | median=%.2f | cutoff line at %.2f",
                              length(vals), stats::median(vals), -log10(qcut)),
           x = xlab, y = "Count") +
      theme_bw()
    ggsave(out_png, p, width = 7, height = 5, dpi = 200)
  }
  
  plot_cdf <- function(qvals, title, xlab, out_png) {
    df <- data.frame(q = as.numeric(qvals))
    p <- ggplot(df, aes(q)) +
      stat_ecdf(geom = "step") +
      geom_vline(xintercept = qcut, linetype = "dashed") +
      scale_x_continuous(trans = "log10",
                         breaks = c(1, 0.5, 0.1, 0.05, 0.01, 0.001),
                         labels = c("1","0.5","0.1","0.05","1e-2","1e-3")) +
      labs(title = title,
           subtitle = sprintf("n=%d | share < 0.1: %.1f%%",
                              sum(is.finite(qvals)),
                              100*mean(qvals < 0.1, na.rm = TRUE)),
           x = xlab, y = "Empirical CDF") +
      theme_bw()
    ggsave(out_png, p, width = 7, height = 5, dpi = 200)
  }
  
  # A-only in B
  if (length(q_in_B_for_Aonly) > 0) {
    plot_hist(
      vals   = .safe_nlog10(q_in_B_for_Aonly),
      title  = sprintf("“%s-only” genes in %s: -log10(q)", labelA, labelB),
      xlab   = sprintf("-log10(q) in %s (for %s-only)", labelB, labelA),
      out_png= file.path(outdir, sprintf("sigOnly_%s_in_%s_hist.png", labelA, labelB))
    )
    plot_cdf(
      qvals  = q_in_B_for_Aonly,
      title  = sprintf("“%s-only” genes in %s: q-value CDF", labelA, labelB),
      xlab   = sprintf("q in %s (for %s-only)", labelB, labelA),
      out_png= file.path(outdir, sprintf("sigOnly_%s_in_%s_cdf.png", labelA, labelB))
    )
  } else {
    message("No ", labelA, "-only genes (relative to ", labelB, ").")
  }
  
  # B-only in A
  if (length(q_in_A_for_Bonly) > 0) {
    plot_hist(
      vals   = .safe_nlog10(q_in_A_for_Bonly),
      title  = sprintf("“%s-only” genes in %s: -log10(q)", labelB, labelA),
      xlab   = sprintf("-log10(q) in %s (for %s-only)", labelA, labelB),
      out_png= file.path(outdir, sprintf("sigOnly_%s_in_%s_hist.png", labelB, labelA))
    )
    plot_cdf(
      qvals  = q_in_A_for_Bonly,
      title  = sprintf("“%s-only” genes in %s: q-value CDF", labelB, labelA),
      xlab   = sprintf("q in %s (for %s-only)", labelA, labelB),
      out_png= file.path(outdir, sprintf("sigOnly_%s_in_%s_cdf.png", labelB, labelA))
    )
  } else {
    message("No ", labelB, "-only genes (relative to ", labelA, ").")
  }
  
  invisible(list(
    n_Aonly = length(A_only_ids),
    n_Bonly = length(B_only_ids),
    frac_Aonly_inB_q_lt_0.1 = mean(q_in_B_for_Aonly < 0.1, na.rm = TRUE),
    frac_Bonly_inA_q_lt_0.1 = mean(q_in_A_for_Bonly < 0.1, na.rm = TRUE)
  ))}

# 11. 

diagnostics_summary <- function(dfA, dfB, setsA, setsB,
                                labels = c("A","B"),
                                qcut = 0.05, lfc = 0.58) {
  # Find overlapping genes between A and B
  mat <- dplyr::inner_join(
    dfA %>% dplyr::select(FBgn, logFC_A = logFC, q_A = qval),
    dfB %>% dplyr::select(FBgn, logFC_B = logFC, q_B = qval),
    by = "FBgn"
  )
  
  # Error metrics on overlapping genes (logFC differences)
  dif <- mat$logFC_A - mat$logFC_B
  mae <- mean(abs(dif), na.rm = TRUE)
  rmse <- sqrt(mean(dif^2, na.rm = TRUE))
  bias <- mean(dif, na.rm = TRUE)
  
  # Sign consistency (direction agreement)
  disagree <- sum(sign(mat$logFC_A) != sign(mat$logFC_B), na.rm = TRUE)
  disagree_rate <- disagree / nrow(mat)
  
  # “Sig-only” analysis: look at q-value distribution in the other dataset
  is_sig_A <- (dfA$qval < qcut) & (abs(dfA$logFC) > lfc)
  is_sig_B <- (dfB$qval < qcut) & (abs(dfB$logFC) > lfc)
  sigA_only <- setdiff(dfA$FBgn[is_sig_A], dfB$FBgn[is_sig_B])
  sigB_only <- setdiff(dfB$FBgn[is_sig_B], dfA$FBgn[is_sig_A])
  
  # A-only genes: their q-values in B
  q_in_B_for_Aonly <- mat$q_B[ match(sigA_only, mat$FBgn) ]
  shareAonly_q_005_01 <- mean(q_in_B_for_Aonly > 0.05 & q_in_B_for_Aonly <= 0.1, na.rm = TRUE)
  shareAonly_q_lt_01  <- mean(q_in_B_for_Aonly <= 0.1, na.rm = TRUE)
  
  # B-only genes: their q-values in A
  q_in_A_for_Bonly <- mat$q_A[ match(sigB_only, mat$FBgn) ]
  shareBonly_q_005_01 <- mean(q_in_A_for_Bonly > 0.05 & q_in_A_for_Bonly <= 0.1, na.rm = TRUE)
  shareBonly_q_lt_01  <- mean(q_in_A_for_Bonly <= 0.1, na.rm = TRUE)
  
  # Jaccard overlap for UP/DOWN sets
  A_up   <- setsA$up;   B_up   <- setsB$up
  A_down <- setsA$down; B_down <- setsB$down
  j_up   <- {u <- union(A_up, B_up); if (length(u)==0) NA_real_ else length(intersect(A_up, B_up))/length(u)}
  j_down <- {u <- union(A_down, B_down); if (length(u)==0) NA_real_ else length(intersect(A_down, B_down))/length(u)}
  
  tibble::tibble(
    pair            = paste0(labels[1], "_vs_", labels[2]),
    n_overlap_genes = nrow(mat),
    mae_logFC       = mae,
    rmse_logFC      = rmse,
    bias_logFC      = bias,
    sign_disagree_n = disagree,
    sign_disagree_rate = round(disagree_rate, 4),
    
    A_sig_n         = sum(is_sig_A, na.rm = TRUE),
    B_sig_n         = sum(is_sig_B, na.rm = TRUE),
    A_only_n        = length(sigA_only),
    B_only_n        = length(sigB_only),
    
    `Aonly_q_in_B_share_0.05_0.1` = round(shareAonly_q_005_01, 4),
    `Aonly_q_in_B_share_<=0.1`    = round(shareAonly_q_lt_01, 4),
    `Bonly_q_in_A_share_0.05_0.1` = round(shareBonly_q_005_01, 4),
    `Bonly_q_in_A_share_<=0.1`    = round(shareBonly_q_lt_01, 4),
    
    jaccard_up   = round(j_up, 4),
    jaccard_down = round(j_down, 4)
  )
}

# ========== read inputs ==========
tut <- read_tutorial(tutorial_csv)
edg <- read_edgeR(edgeR_csv)
pub <- read_published(paper_csv)

# background set: union of all FBgn lists
background_all <- unique(c(tut$FBgn, edg$FBgn, pub$FBgn))
length(background_all)

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
  
  # 1) Split significant genes based on thresholds
  setsA <- split_up_down(A, logfc_cut, qval_cut)
  setsB <- split_up_down(B, logfc_cut, qval_cut)
  
  # 2) Venn diagrams (UP/DOWN)
  save_venn(setsA$up,   setsB$up,   c(paste(A_name,"UP"),   paste(B_name,"UP")),
            file.path(outdir, "venn_up.png"))
  save_venn(setsA$down, setsB$down, c(paste(A_name,"DOWN"), paste(B_name,"DOWN")),
            file.path(outdir, "venn_down.png"))
  
  # 3) Fisher’s exact test (return objects, not CSVs)
  fres_up   <- fisher_csv(setsA$up,   setsB$up,   background_all,
                          direction="UP", A_name=A_name, B_name=B_name)
  fres_down <- fisher_csv(setsA$down, setsB$down, background_all,
                          direction="DOWN", A_name=A_name, B_name=B_name)
  
  # 4) Jaccard index and overlap sizes (for Excel)
  j_up   <- jaccard(setsA$up,   setsB$up)
  j_down <- jaccard(setsA$down, setsB$down)
  size_tbl <- tibble(
    set = c(paste0(A_name,"_UP"), paste0(B_name,"_UP"), "Overlap_UP",
            paste0(A_name,"_DOWN"), paste0(B_name,"_DOWN"), "Overlap_DOWN"),
    size = c(length(setsA$up), length(setsB$up), length(intersect(setsA$up, setsB$up)),
             length(setsA$down), length(setsB$down), length(intersect(setsA$down, setsB$down)))
  )
  
  # 5) Pairwise table of genes with inconsistent directions
  pairs_disagree <- pairwise_table_disagree(A, B,
                                            A_name = A_name, B_name = B_name,
                                            qcut = qval_cut, lfc = logfc_cut,
                                            outdir = outdir)
  
  # 6) logFC correlation plots
  plot_logFC_correlation(A, B, labels=c(A_name, B_name),
                         outdir=outdir, prefix=paste0(A_name, "_vs_", B_name))
  
  # 7) Volcano plot (B can be passed directly as data.frame)
  volcano_mark_overlap(A, B, qcut=qval_cut, lfc=logfc_cut,
                       labelA=A_name, labelB=B_name,
                       outpng=file.path(outdir, paste0("volcano_", A_name, "_mark_", B_name, ".png")))
  
  # 8) Histograms and CDFs for "sig-only" genes
  plot_sigonly_hist_and_cdf(dfA=A, dfB=B,
                            labelA=A_name, labelB=B_name,
                            qcut=qval_cut, lfc=logfc_cut,
                            outdir=outdir)
  
  # 9) Diagnostics summary (tibble)
  diag_tbl <- diagnostics_summary(A, B, setsA, setsB,
                                  labels = c(A_name, B_name),
                                  qcut = qval_cut, lfc = logfc_cut)

  
  # ---------- Excel summary: create workbook first, then write ----------
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(wb, "fisher_up")
  openxlsx::writeData(wb, "fisher_up", fres_up$table)
  
  openxlsx::addWorksheet(wb, "fisher_down")
  openxlsx::writeData(wb, "fisher_down", fres_down$table)
  
  openxlsx::addWorksheet(wb, "overlap_sizes")
  openxlsx::writeData(wb, "overlap_sizes", size_tbl)
  
  openxlsx::addWorksheet(wb, "jaccard")
  openxlsx::writeData(wb, "jaccard",
                      tibble(jaccard_up = j_up, jaccard_down = j_down))
  
  # Intersection and "only" gene lists
  inter_up   <- sort(intersect(setsA$up,   setsB$up))
  onlyA_up   <- sort(setdiff(setsA$up,     setsB$up))
  onlyB_up   <- sort(setdiff(setsB$up,     setsA$up))
  inter_down <- sort(intersect(setsA$down, setsB$down))
  onlyA_down <- sort(setdiff(setsA$down,   setsB$down))
  onlyB_down <- sort(setdiff(setsB$down,   setsA$down))
  
  openxlsx::addWorksheet(wb, "intersect_UP")
  openxlsx::writeData(wb, "intersect_UP",   tibble(FBgn=inter_up))
  openxlsx::addWorksheet(wb, "only_A_UP")
  openxlsx::writeData(wb, "only_A_UP",      tibble(FBgn=onlyA_up))
  openxlsx::addWorksheet(wb, "only_B_UP")
  openxlsx::writeData(wb, "only_B_UP",      tibble(FBgn=onlyB_up))
  openxlsx::addWorksheet(wb, "intersect_DOWN")
  openxlsx::writeData(wb, "intersect_DOWN", tibble(FBgn=inter_down))
  openxlsx::addWorksheet(wb, "only_A_DOWN")
  openxlsx::writeData(wb, "only_A_DOWN",    tibble(FBgn=onlyA_down))
  openxlsx::addWorksheet(wb, "only_B_DOWN")
  openxlsx::writeData(wb, "only_B_DOWN",    tibble(FBgn=onlyB_down))
  
  # Disagree-only table (if exists)
  if (nrow(pairs_disagree) > 0) {
    openxlsx::addWorksheet(wb, "pairs_disagree")
    openxlsx::writeData(wb, "pairs_disagree", pairs_disagree)
  }
  
  # Diagnostics
  openxlsx::addWorksheet(wb, "Diagnostics")
  openxlsx::writeData(wb, "Diagnostics", diag_tbl)
  
  # Save Excel workbook
  xlsx_path <- file.path(outdir, paste0(A_name, "_vs_", B_name, "_Summary.xlsx"))
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  message("Saved workbook -> ", xlsx_path)
  rm(wb)
}
