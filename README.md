# OAS3 in Breast Cancer: A TCGA-BRCA Bioinformatics Analysis

## Overview

This project investigates the expression and biological relevance of **oligoadenylate synthetase 3 (OAS3)** in breast cancer using publicly available data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA).

The project was inspired by the pan-cancer study:

> Zhang X, Zhang H, Zhong L, et al. *Pan-cancer multi-omics profiling of OAS3 reveals its immunological and prognostic associations across human cancers.* PeerJ. 2026. DOI: 10.7717/peerj.20805.

The original study investigated OAS3 across multiple cancer types and reported associations with cancer expression, prognosis, immune regulation, interferon response, and other molecular features.

Rather than attempting to reproduce the entire pan-cancer analysis, this project performs a **breast cancer-specific computational extension** focused on TCGA-BRCA.

---

## Research Question

> **What is the expression pattern, prognostic relevance, and molecular context of OAS3 in breast cancer?**

The analysis addresses five main questions:

1. Is OAS3 differentially expressed between breast tumor and normal tissue?
2. Is OAS3 expression associated with overall survival?
3. Which genes are most strongly co-expressed with OAS3?
4. What biological processes and pathways are associated with OAS3?
5. What genes and pathways differ between OAS3-high and OAS3-low breast tumors?

---

## Key Findings

### 1. OAS3 is elevated in breast tumor tissue

OAS3 expression was significantly higher in Primary Tumor samples than in Solid Tissue Normal samples.

| Sample type         |     n | Median OAS3 |
| ------------------- | ----: | ----------: |
| Primary Tumor       | 1,106 |        5.46 |
| Solid Tissue Normal |   113 |        4.41 |

Wilcoxon rank-sum test:

**P < 2.2 × 10⁻¹⁶**

---

### 2. OAS3 expression was not significantly associated with overall survival

The survival analysis included 1,203 patients and 197 observed events.

Cox proportional hazards analysis:

**HR = 1.20 (95% CI: 0.90–1.59), P = 0.21**

Therefore, this analysis did not identify a statistically significant overall survival difference between the OAS3-high and OAS3-low groups.

This contrasts with the unfavorable prognostic associations reported for OAS3 in selected cancer types in the original pan-cancer study.

---

### 3. OAS3 is strongly associated with interferon/antiviral genes

The strongest positively co-expressed genes included:

* OAS2
* OAS1
* PARP9
* CMPK2
* STAT1
* IFIT3
* DDX60
* IFIH1
* PARP14
* IFIT1
* DDX58
* RSAD2
* IFIT2
* OASL
* USP18
* MX1

This suggests that OAS3 expression is embedded within a broader interferon-stimulated and antiviral transcriptional program.

---

### 4. Functional enrichment supports antiviral and immune biology

GO enrichment identified strong enrichment for:

* Defense response to virus
* Response to virus
* Negative regulation of viral process
* Regulation of viral process
* Regulation of viral genome replication

KEGG enrichment also identified several viral disease pathways.

These KEGG results are interpreted as enrichment of **shared host antiviral and interferon-related genes**, rather than evidence of specific viral infection in breast cancer.

---

### 5. OAS3-high and OAS3-low tumors show distinct transcriptional profiles

Differential expression analysis identified:

**116 significantly differentially expressed genes**

using:

```text
Adjusted P-value < 0.05
Absolute log2 fold change > 1
```

GO enrichment of these genes again highlighted antiviral and interferon-associated biological processes.

---

## Research Workflow

```text
TCGA-BRCA Data
       │
       ├── RNA-seq Expression
       │
       ├── Phenotype
       │
       └── Survival
       │
       ▼
Data Preprocessing
       │
       ▼
OAS3 Expression
       │
       ├── Tumor vs Normal
       │
       ├── Survival Analysis
       │
       ├── Co-expression
       │       │
       │       └── GO / KEGG
       │
       └── OAS3 High vs Low
               │
               └── Differential Expression
                       │
                       └── GO / KEGG
```

---

## Data

The analysis uses publicly available TCGA-BRCA data accessed through **UCSC Xena**.

### Expression

**TCGA-BRCA STAR-TPM**

* 60,660 gene identifiers
* 1,226 samples
* Log2-transformed TPM expression

### Phenotype

TCGA-BRCA clinical and sample-level phenotype information.

### Survival

TCGA-BRCA overall survival information.

### Gene Annotation

GENCODE v36 gene annotation/probemap.

Raw datasets are not included in this repository because of their large file sizes.

---

## Repository Structure

```text
OAS3-BRCA-Bioinformatics/
│
├── data/
│   ├── raw/
│   │   └── README.md
│   │
│   └── processed/
│       └── README.md
│
├── scripts/
│   ├── 01_data_import.R
│   ├── 02_data_preprocessing.R
│   ├── 03_OAS3_expression.R
│   ├── 04_survival_analysis.R
│   ├── 05_coexpression_analysis.R
│   ├── 06_functional_enrichment.R
│   └── 07_differential_expression.R
│
├── figures/
│   ├── Figure1_OAS3_tumor_vs_normal.png
│   ├── Figure2_Kaplan_Meier.pdf
│   ├── Figure3_OAS3_coexpression.png
│   ├── Figure4_GO_Dotplot.pdf
│   ├── Figure5_KEGG_Dotplot.pdf
│   ├── Figure6_volcanoplot_OAS3_High_vs_Low.png
│   ├── Figure7_GO_DEG_Dotplot.pdf
│   └── Figure8_KEGG_DEG_Dotplot.pdf
│
├── results/
│   ├── top30_OAS3_coexpressed_genes.csv
│   ├── GO_results.csv
│   ├── KEGG_results.csv
│   ├── GO_DEG_results.csv
│   └── KEGG_DEG_results.csv
│
├── METHODS.md
├── DISCUSSION.md
├── README.md
└── .gitignore
```

*File names may differ slightly depending on the final organization of the analysis scripts and figures.*

---

## Methods

A detailed description of the computational workflow, data preprocessing, statistical analyses, differential expression analysis, and functional enrichment methods is available in:

**[METHODS.md](METHODS.md)**

---

## Results and Discussion

Interpretation of the findings, comparison with the original pan-cancer study, limitations, and future directions are provided in:

**[DISCUSSION.md](DISCUSSION.md)**

---

## Technologies and R Packages

The analysis was conducted using:

* R 4.6.1
* RStudio
* TCGA / UCSC Xena
* `readr`
* `dplyr`
* `tibble`
* `ggplot2`
* `survival`
* `survminer`
* `limma`
* `clusterProfiler`
* `org.Hs.eg.db`
* `TCGAbiolinks`

---

## Reproducibility

To reproduce the analysis:

1. Clone or download this repository.
2. Obtain the TCGA-BRCA datasets from UCSC Xena.
3. Place the downloaded raw datasets in `data/raw/`.
4. Install the required R packages.
5. Run the scripts sequentially.
6. Figures will be generated in `figures/`.
7. Analysis tables will be generated in `results/`.

Raw TCGA datasets are excluded from Git tracking because of their large size.

---

## Relationship to the Original Study

This project is a **breast cancer-specific computational extension** of the published pan-cancer study rather than a full replication.

The original study examined OAS3 across multiple cancer types and incorporated expression, survival, immune infiltration, immune checkpoints, stemness, mutation, genomic, and experimental analyses.

This project focuses on a subset of those biological questions within TCGA-BRCA:

* OAS3 expression
* Tumor-versus-normal comparison
* Overall survival
* Co-expression
* Differential expression
* GO enrichment
* KEGG enrichment

The purpose is to determine whether selected observations from the pan-cancer analysis are also evident specifically in breast cancer.

---

## Limitations

The main limitations of this analysis include:

* Unequal numbers of tumor and normal samples.
* Reliance on retrospective TCGA data.
* Lack of comprehensive clinical adjustment in the survival analysis.
* Correlation-based co-expression analysis cannot establish causality.
* KEGG viral pathway enrichment does not indicate viral infection.
* Findings require validation in independent breast cancer cohorts.

---

## Future Directions

Potential extensions include:

* Breast cancer molecular subtype analysis.
* OAS3 expression across pathological stages.
* Multivariable survival analysis.
* Independent cohort validation.
* Immune infiltration analysis.
* Immune checkpoint analysis.
* Treatment-response analysis.
* Gene Set Enrichment Analysis (GSEA).
* Protein-level validation.

---

## Citation

If you use this analysis or repository, please cite the original study that inspired the project:

**Zhang X, Zhang H, Zhong L, Yang W, Yang C, Pan Y, Liu B.**
*Pan-cancer multi-omics profiling of OAS3 reveals its immunological and prognostic associations across human cancers.* PeerJ. 2026;14:e20805.
DOI: 10.7717/peerj.20805.

The original article is published under a **Creative Commons CC BY 4.0** license.

This repository contains an independent computational analysis and should not be considered a reproduction of the original authors' complete analysis.

---

## Disclaimer

This project is intended for research, educational, and bioinformatics portfolio purposes. The findings should not be interpreted as clinical evidence or used for patient diagnosis, prognosis, or treatment decisions.
