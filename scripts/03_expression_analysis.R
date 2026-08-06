library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Load processed data
expr_annotated <- readRDS(
  "data/processed/expr_annotated.rds"
)

pheno_clean <- readRDS(
  "data/processed/pheno_clean.rds"
)

expr_mat <- readRDS(
  "data/processed/expr_mat.rds"
)

#MATCH expression samples with phenotype
expr_samples <- colnames(expr_mat)

length(intersect(expr_samples, pheno_clean$sample))

#CREATE A SIMPLE METADATA
sample_info <- pheno_clean %>%
  filter(sample %in% colnames(expr_mat))
#extract OAS3 expression
oas3 <- expr_annotated %>%
  filter(gene == "OAS3")
oas3[,1:10]
#CONVERT TO A LONG TABLE FORMAT
library(tidyr)

oas3_long <- oas3 %>%
  select(gene, all_of(colnames(expr_mat))) %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  )
head(oas3_long)
dim(oas3_long)
#MERGE WITH THE PHENOTYPE
oas3_data <- left_join(
  oas3_long,
  sample_info,
  by = "sample"
)
head(oas3_data)
dim(oas3_data)
#use only tumor and normal samples
oas3_tn <- oas3_data %>%
  filter(sample_type.samples %in% c("Primary Tumor",
                                    "Solid Tissue Normal"))
table(oas3_tn$sample_type.samples)

#do the wilcox test
wilcox.test(
  expression ~ sample_type.samples,
  data = oas3_tn
)
#to calculate the median expression, to know if OAS3 gene is downregulated or upregulated
oas3_tn %>%
  group_by(sample_type.samples) %>%
  summarise(
    n = n(),
    median = median(expression),
    mean = mean(expression),
    sd = sd(expression)
  )
#to create my plot
library(ggplot2)

ggplot(oas3_tn,
       aes(x = sample_type.samples,
           y = expression,
           fill = sample_type.samples)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2,
              alpha = 0.25,
              size = 0.8) +
  labs(
    title = "OAS3 Expression in TCGA-BRCA",
    x = "",
    y = "log2(TPM + 1)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
ggsave(
  "figures/Figure1_OAS3_Tumour_vs_Normal.png",
  width = 6,
  height = 5,
  dpi = 300
)
