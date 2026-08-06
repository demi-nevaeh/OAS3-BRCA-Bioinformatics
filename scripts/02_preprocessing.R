library(readr)

expr <- read_tsv("data/raw/TCGA-BRCA.star_tpm.tsv.gz")
dim(expr)
head(expr[,1:5])
colnames(expr)[1:10]
#rename the gene column
colnames(expr)[1] <- "gene_id"
#set gene_id column as rows
library(dplyr)
library(tibble)
expr_mat <- expr %>%
  column_to_rownames("gene_id")
dim(expr_mat)
dim(expr)
head(expr_mat[, 1:5])

#ANNOTATE GENES

#to download probemap to annotate gene
download.file(
  "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
  destfile = "data/raw/gencode.v36.annotation.gtf.gene.probemap",
  mode = "wb"
)
library(readr)

probemap <- read_tsv(
  "data/raw/gencode.v36.annotation.gtf.gene.probemap"
)
head(probemap)
colnames(probemap)
#to check if the respective genes names match
head(rownames(expr_mat))
head(probemap$id)
#to create a copy before merging
expr_df <- data.frame(
  gene_id = rownames(expr_mat),
  expr_mat,
  check.names = FALSE
)
dim(expr_df)
#merging with the annotation
library(dplyr)

expr_annotated <- left_join(
  expr_df,
  probemap,
  by = c("gene_id" = "id")
)
dim(expr_annotated)
sum(is.na(expr_annotated$gene))
head(expr_annotated[, c("gene_id", "gene")])

oas3 <- expr_annotated %>%
  filter(gene == "OAS3")

oas3
expr_annotated %>%
  filter(gene %in% c("OAS1", "OAS2", "OAS3"))

dir.create("data/processed", showWarnings = FALSE)

saveRDS(
  expr_annotated,
  "data/processed/expr_annotated.rds"
)

oas_family <- expr_annotated %>%
  filter(gene %in% c(
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL"
  ))

oas3 <- expr_annotated %>%
  filter(gene == "OAS3")

oas3[, 1:10]
#UPLOADING THE PHENOTYPE
library(readr)

pheno <- read_tsv("data/raw/TCGA-BRCA.clinical.tsv.gz")
dim(pheno)
head(pheno)
colnames(pheno)
#UNDERSTAND THE PHENOTYPE DATA
table(pheno$sample_type.samples)
head(pheno$sample)

#TO CREATE A PHENOTYPE TABLE
library(dplyr)

pheno_clean <- pheno %>%
  select(
    sample,
    sample_type.samples,
    gender.demographic,
    age_at_index.demographic,
    vital_status.demographic,
    ajcc_pathologic_stage.diagnoses
  ) %>%
  distinct()

dim(pheno_clean)
head(pheno_clean)

saveRDS(expr_annotated,
        "data/processed/expr_annotated.rds")

saveRDS(pheno_clean,
        "data/processed/pheno_clean.rds")
saveRDS(expr_mat,
        "data/processed/expr_mat.rds")
