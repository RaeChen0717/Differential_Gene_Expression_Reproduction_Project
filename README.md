# Differential Gene Expression (DGE) Comparison Summary

This project compares **three pipelines** for RNA-seq DGE analysis in *Drosophila*:  

1. **Tutorial pipeline (kallisto + sleuth)**  
2. **edgeR pipeline (HISAT2 + featureCounts + edgeR)**  
3. **Published paper pipeline (HISAT2 + HTSeq + edgeR)**  

---

## Results Overview

| Comparison               | Overlap of Significant DEGs | Pearson r | Spearman ρ | Interpretation |
|--------------------------|-----------------------------|-----------|------------|----------------|
| **Tutorial vs edgeR**    | ~70%                        | 0.91      | 0.96       | Reasonably high correlation. Most DEGs consistent, but ~30% differ due to filtering and model choice. |
| **edgeR vs Published**   | ~66%                        | 0.99      | 1.00       | Almost perfect correlation. Confirms edgeR analysis reproduces published results faithfully. |
| **Tutorial vs Published**| ~59%                        | 0.76      | 0.90       | Moderate overlap and lower correlation. Tutorial pipeline diverges more because it uses pseudo-alignment and sleuth’s Wald test. |

---

## Explanations

### 1. **Tutorial vs edgeR**
- **Correlation**: High (Pearson 0.91, Spearman 0.96).  
- **Overlap**: ~70% of DEGs shared.  
- **Differences**:  
  - Sleuth uses bootstrap-based variance estimation, edgeR uses negative binomial GLM.  
  - Different filtering of low-expression genes (sleuth more lenient, edgeR stricter).  
  - Some borderline genes pass threshold in one but not the other.  

**Volcano Plot**: shows many shared DEGs, but also "sig only in Tutorial" or "sig only in edgeR" points. These reflect model-dependent borderline calls.

---

### 2. **edgeR vs Published**
- **Correlation**: Almost perfect (Pearson 0.99, Spearman 1.00).  
- **Overlap**: ~66% DEGs identical.  
- **Interpretation**:  
  - Confirms your re-analysis with edgeR essentially replicated the paper.  
  - Small differences likely due to version changes (edgeR 3.16.5 in paper vs newer in your run), or filtering (mean CPM threshold).  

**Venn Diagrams**: show ~2/3 overlap, which is strong given minor pipeline differences.

---

### 3. **Tutorial vs Published**
- **Correlation**: Moderate (Pearson 0.76, Spearman 0.90).  
- **Overlap**: ~59% DEGs shared.  
- **Why lower?**  
  - Paper used alignment-based counting (HISAT2 + HTSeq).  
  - Tutorial used pseudo-alignment (kallisto + sleuth).  
  - Statistical models differ (Wald test vs quasi-likelihood F-test).  
  - Filtering thresholds differ, especially for low CPM genes.  

**Scatter Plot**: shows general trend agreement (logFC direction consistent), but more variance at extremes.

---

## Biological Meaning

- **High Pearson/Spearman r** across all comparisons means the **direction and magnitude of log2FC are consistent** across pipelines.  
- The **lower Jaccard overlap** (40–70%) is expected — DE calls depend on statistical thresholding and filtering, not just raw effect size.  
- **Key conclusion**:  
  - edgeR best reproduces the published pipeline.  
  - Tutorial pipeline is a reasonable approximation but not identical.  
  - Divergence highlights the sensitivity of DEG lists to **analysis choices**.

---

## Figures (examples)

- **Venn Diagrams** → show DEG overlaps between methods.  
- **Scatter Plots** → demonstrate strong linear correlation of log2FC values.  
- **Volcano Plots** → highlight how many significant DEGs are shared vs unique to each method.  

---

## Take-Home Messages
1. **Consistency**: All pipelines agree on overall direction of regulation.  
2. **edgeR vs Published**: Nearly identical — validates your re-analysis.  
3. **Tutorial differences**: Mainly from pseudo-alignment vs alignment and statistical frameworks.  
4. **Biological robustness**: Core biological findings are reproduced despite pipeline differences.  
