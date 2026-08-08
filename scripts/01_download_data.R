install.packages(c(
  "tidyverse",
  "survival",
  "survminer",
  "ggplot2",
  "SummarizedExperiment"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "limma",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db"
))
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(limma)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)

# Prepare the downloaded data
brca <- GDCprepare(query)

class(brca)
