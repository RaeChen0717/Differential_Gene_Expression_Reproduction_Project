# DGE Comparison Summary

This report evaluates the consistency of differential gene expression (DGE) analysis across three datasets/pipelines:

- **Tutorial pipeline (sleuth-based)**
- **edgeR re-analysis**
- **Published dataset (Van Swinderen et al., eLife 2023)**

Thresholds used: **FDR < 0.05** and **|log2FC| > 0.58**.  
Focus is on DEG overlap, correlation of log2FC values, and analysis of “sig-only” genes.

---

## 1. Tutorial vs Published

**DEG overlap (from Venn diagrams):**
- UP: **133 shared**, 93 Tutorial-only, 83 Published-only → ~59% overlap  
- DOWN: **98 shared**, 60 Tutorial-only, 73 Published-only → ~57% overlap  

**Log2FC correlation (Tutorial_vs_Published_correlation_stats.txt):**
- Pearson r = **0.91**, 95% CI [0.902–0.910]  
- Spearman ρ = **0.95**  
- Both p < 2.2e-16 → strong agreement  

**Sig-only diagnostics (CDF/Hist plots):**
- Tutorial-only DEGs: many q-values in Published ≈ 0.05–0.1  
- Published-only DEGs: similar, cluster near cutoff  
  Suggests differences come mainly from **borderline q-values**, not direction disagreement.  

**Volcano plot (volcano_Tutorial_mark_Published.png):**  
- Shared DEGs (green) dominate in both tails  
- Sig-only genes (red) cluster near cutoff thresholds  

---

## 2. Tutorial vs edgeR

**DEG overlap (Venn diagrams):**
- UP: **89 shared**, 75 edgeR-only, 24 Tutorial-only → ~50% overlap  
- DOWN: **69 shared**, 44 Tutorial-only, 15 edgeR-only → ~50% overlap  

**Log2FC correlation (Tutorial_vs_edgeR_correlation_stats.txt):**
- Pearson r = **0.91**, 95% CI [0.902–0.910]  
- Spearman ρ = **0.96**  
- p < 2.2e-16 → high consistency  

**Sig-only diagnostics (Hist/CDF):**
- EdgeR-only DEGs: n=67, **73.1%** have q < 0.1 in Tutorial  
- Tutorial-only DEGs: n=38, **92.1%** have q < 0.1 in edgeR  
  Again, differences are **borderline significance calls**, not direction flips.  

**Volcano plot (volcano_Tutorial_mark_edgeR.png):**  
- Shared DEGs clearly beyond thresholds  
- Sig-only DEGs cluster at edges of cutoff lines  

---

## 3. edgeR vs Published

**DEG overlap (Venn diagrams):**
- UP: **164 shared**, 84 edgeR-only, 104 Published-only → ~61% overlap  
- DOWN: **144 shared**, 55 edgeR-only, 56 Published-only → ~72% overlap  

**Log2FC correlation (edgeR_vs_Published_correlation_stats.txt):**
- Pearson r = **0.994**, 95% CI [0.9939–0.9945]  
- Spearman ρ = **0.998**  
- p < 2.2e-16 → **nearly perfect correlation**  

**Sig-only diagnostics (Hist/CDF):**
- EdgeR-only DEGs: majority have Published q between 0.05–0.1  
- Published-only DEGs: similar pattern  
  Most “sig-only” genes are **borderline threshold artifacts**.  

**Volcano plot (volcano_edgeR_mark_Published.png):**  
- Shared DEGs dominate  
- Sig-only DEGs appear near FDR/logFC cutoffs  

---

## Diagnostics & Interpretation

1. **Strong log2FC correlations** across all comparisons:
   - Tutorial vs Published: r ≈ 0.91  
   - Tutorial vs edgeR: r ≈ 0.91  
   - edgeR vs Published: r ≈ 0.994  
   → All pipelines estimate fold changes consistently.  

2. **Overlap rates (50–72%) are below the expected 85%**, but:
   - Volcano, histograms, and CDF plots show that “sig-only” genes usually lie near cutoffs.  
   - These discrepancies arise from small differences in filtering, dispersion modeling, or FDR adjustments.  

3. **Direction consistency:**  
   - Almost no genes flip sign between methods.  
   - Disagreement comes from **significance thresholding**, not opposite biological conclusions.  

---
# Tutorial vs edgeR Comparison Outputs

This folder contains all results from the **Tutorial vs edgeR differential expression analysis (DGE) comparison**.  
Both statistical summaries and visualizations are included.  

---

## 📊 Summary File
- **Tutorial_vs_edgeR_Summary.xlsx**  
  Main Excel summary containing:
  - Fisher’s exact test results (UP/DOWN sets)  
  - Overlap sizes and Jaccard indices  
  - Intersect / only gene lists  
  - Diagnostics summary (MAE, RMSE, bias, sign disagreement rate, sig-only q-value analysis, etc.)  

---

## 🔍 Volcano Plot
- **volcano_Tutorial_mark_edgeR.png**  
  Volcano plot based on Tutorial results, highlighting overlap with edgeR:  
  - **Grey** = Not significant in Tutorial  
  - **Red** = Significant only in Tutorial or edgeR  
  - **Green** = Significant in both Tutorial and edgeR  

---

## 📈 Correlation Analysis
- **Tutorial_vs_edgeR_logFC_correlation.png**  
  Scatter plot of log2 fold changes (logFC) between Tutorial and edgeR.  
  Includes linear regression fit and correlation statistics.  

- **Tutorial_vs_edgeR_correlation_stats.txt**  
  Text file with detailed correlation metrics:  
  - Pearson correlation (r, CI, p-value)  
  - Spearman correlation (ρ, p-value)  

---

## ⚠️ Disagreement Table
- **pairs_disagree_Tutorial_vs_edgeR.csv**  
  List of genes where **Tutorial and edgeR disagree in logFC direction**.  
  Columns include FBgn, logFC/qval from both datasets, category, and delta logFC.  

---

## 🔵 Venn Diagrams
- **venn_up.png**  
  Venn diagram showing overlap of **up-regulated genes**.  
- **venn_down.png**  
  Venn diagram showing overlap of **down-regulated genes**.  

---

## 📊 Sig-only Gene Diagnostics
To investigate why overlap is below 100%, we analyzed “sig-only” genes.  

- **sigOnly_edgeR_in_Tutorial_cdf.png**  
  CDF plot of q-values for edgeR-only genes in Tutorial.  
- **sigOnly_edgeR_in_Tutorial_hist.png**  
  Histogram of –log10(q) values for edgeR-only genes in Tutorial.  

- **sigOnly_Tutorial_in_edgeR_cdf.png**  
  CDF plot of q-values for Tutorial-only genes in edgeR.  
- **sigOnly_Tutorial_in_edgeR_hist.png**  
  Histogram of –log10(q) values for Tutorial-only genes in edgeR.  

> These plots show that most sig-only genes fall **near the FDR cutoff boundary (q ≈ 0.05)**, suggesting that discrepancies arise primarily from threshold effects rather than large differences in effect size or direction.  

---

✅ Together, these results provide a comprehensive comparison between Tutorial and edgeR pipelines in identifying differentially expressed genes (DEGs).
## Conclusion

- Despite moderate Venn overlap (<85%), the **biological conclusions are consistent**.  
- Differences are mainly **statistical cutoff artifacts**, not conflicting fold-change directions.  
- **edgeR vs Published shows the highest agreement (r=0.994, >70% overlap)**, confirming that the pipeline reproduces the published results.  
- Tutorial pipeline shows slightly lower overlap with both edgeR and Published, reflecting differences in statistical modeling.  

---
