
---

```markdown
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

Allele frequency shifts of significant variants:  
![Allele Frequency Bar Plot](results/allele_freq_shifts.png)  

Variant impact summary:  
![Impact Distribution](results/variant_impact.png)  

---


---

## 📌 Conclusion
- African ancestry is associated with **27 statistically significant allele frequency shifts** within cancer-related regions of Chromosome 1, including **3 HIGH-impact variants**.  
- No population-exclusive alleles were detected, but **Eurocentric reference panels can mis-rank African variants**.  
- This study demonstrates a **scalable, public-data pipeline** for ancestry-aware variant interpretation.  
- Highlights the need for **inclusive global population panels** in cancer genetics.  

---

## 📂 Repository Structure
```

cancer-chr1-analysis/
│── data/          # Processed datasets, frequency files, panel info (raw excluded)
│── scripts/       # R and Python scripts for QC, PCA, stats, and plots
│── results/       # Output figures, PCA plots, clustering, and tables
│── docs/          # Abstract, notes, supplementary files
│── backup/        # Temporary backup files (not tracked in GitHub)
│── ensembl-vep/   # External tool folder (not tracked in GitHub)
│── README.md      # Project overview and usage instructions
│── .gitignore     # Excluded large or unnecessary files

````

---

## ⚙️ Reproducibility
### Requirements
- **R (≥4.0)** with:
  - `tidyverse`, `data.table`, `ggplot2`, `cluster`, `factoextra`
- **Python (optional, for plots/scripts):**
  - `pandas`, `numpy`, `matplotlib`, `seaborn`
- **Other tools:**
  - [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) for variant annotation  
  - [PLINK](https://www.cog-genomics.org/plink/) for population genetics analyses  

### Steps
1. Clone this repository:
   ```bash
   git clone https://github.com/BioInfoma/Cancer_Variants_in-Chromosome_1.git
   cd Cancer_Variants_in-Chromosome_1
````

2. Place processed or subset data into `data/` (raw 1000 Genomes files excluded).

3. Run the main R script:

   ```bash
   Rscript scripts/Chr1_cancer_analysis.R
   ```

4. Some of the Outputs (plots, tables, PCA results) will be saved in the `results/` folder.

---

## 📎 Notes

* **Raw VCFs from 1000 Genomes are not uploaded** due to size.
* This repository contains **processed data, scripts, and reproducible workflows**.
* For full raw datasets, see the [1000 Genomes Project](https://www.internationalgenome.org/).

---

## 🔑 Keywords

allele frequency · cancer genetics · African ancestry · population genomics · Chromosome 1 · 1000 Genomes Project

```

---
 

👉 Do you also want me to draft a **short project description (2–3 sentences)** you can use for your GitHub repo tagline and LinkedIn post?
```
