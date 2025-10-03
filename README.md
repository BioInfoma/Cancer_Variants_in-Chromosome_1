# Cancer-linked Chromosome 1 Variants: African Ancestry-associated Frequency Shifts in 1000 Genomes

**Authors:**  
Temitope Kolawole¹*, Grace Adeloye¹, Ifedayo Fashola¹, Esther Orji¹  
¹ Department of Biological Sciences, Redeemer’s University, Nigeria  

*Corresponding author: kolawole17791@run.edu.ng  

---

## 📌 Background
Most cancer genetics studies rely heavily on allele-frequency data from European populations. This creates a bias where variants common in non-European groups, particularly Africans, may be overlooked or mis-ranked.  

This project investigates whether germline variants in **cancer-associated regions of Chromosome 1** show ancestry-related allele frequency patterns, with a focus on African populations, using the **1000 Genomes Project Phase 3 dataset**.

---

## 🧪 Methods
- **Data source:** 1000 Genomes Phase 3, 2,504 individuals across AFR, AMR, EAS, EUR, SAS populations.  
- **Intervals studied:** 8 cancer-linked regions on Chromosome 1.  
- **Variant processing:**  
  - Quality control retained **264 biallelic SNVs**  
  - Annotation performed with **Ensembl VEP**  
  - 12 variants (4.5%) classified as **HIGH-impact**  
- **Population analysis:**  
  - PCA to assess global structure  
  - χ² tests for population-level differences  
  - Fisher’s exact tests (AFR vs non-AFR) with FDR correction  
- **Effect sizes:**  
  - ΔAF (absolute allele frequency change)  
  - Odds ratios  
  - Cohen’s h  
- **Enrichment definition:** ΔAF ≥ 0.10 (sensitivity threshold at ΔAF ≥ 0.05)

---

## 📊 Results
- **27 variants (10.2%)** showed significant AFR vs non-AFR allele frequency differences (*FDR < 0.05*).  
- Of these, **3 were HIGH-impact variants**.  
- **3 variants** reached ΔAF ≥ 0.10, though none were HIGH-impact.  
- No AFR-private variants (≥1% AFR, <0.1% elsewhere) were identified → differences reflect **frequency skews** not exclusive alleles.  
- PCA showed clear AFR vs non-AFR separation, with small but statistically significant global structure (*χ² = 33,833; p < 2.2×10⁻¹⁶; Cramér’s V = 0.017*).  

---

### 🖼️ Example Figures
PCA plot of Chromosome 1 cancer-linked variants across global populations:  
![PCA of AFR vs non-AFR](results/pca_chr1.png)  

Allele frequency shifts of significant variants:  
![Allele Frequency Bar Plot](results/allele_freq_shifts.png)  

---

## 📌 Conclusion
- African ancestry is associated with **27 statistically significant allele frequency shifts** within cancer-related regions of Chromosome 1, including **3 HIGH-impact variants**.  
- No population-exclusive alleles were detected, but **Eurocentric reference panels can mis-rank African variants**.  
- This study demonstrates a **scalable, public-data pipeline** for ancestry-aware variant interpretation.  
- Highlights the need for **inclusive global population panels** in cancer genetics.  

---
