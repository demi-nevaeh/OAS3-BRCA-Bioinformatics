oas3_surv <- readRDS(
  "data/processed/oas3_surv.rds"
)

expr_mat <- readRDS(
  "data/processed/expr_mat.rds"
)


metadata <- oas3_surv %>%
  dplyr::select(sample, OAS3_group)

head(metadata)

common_samples <- intersect(colnames(expr_mat), metadata$sample)

length(common_samples)

expr_dge <- expr_mat[, common_samples]

metadata <- metadata %>%
  filter(sample %in% common_samples)

metadata <- metadata[
  match(colnames(expr_dge), metadata$sample),
]
all(colnames(expr_dge) == metadata$sample)
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma")

library(limma)

group <- factor(metadata$OAS3_group)

design <- model.matrix(~ group)

design

fit <- lmFit(as.matrix(expr_dge), design)

fit <- eBayes(fit)

deg <- topTable(
  fit,
  coef = 2,
  number = Inf
)

head(deg)

deg$gene_id <- rownames(deg)

annotation <- expr_annotated %>%
  dplyr::select(gene_id, gene) %>%
  distinct()

deg <- left_join(deg, annotation, by = "gene_id")
deg_sig <- deg %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )
nrow(deg_sig)
write.csv(
  deg,
  "results/OAS3_High_vs_Low_DEGs.csv",
  row.names = FALSE
)
library(ggplot2)

ggplot(
  deg,
  aes(logFC, -log10(adj.P.Val))
) +
  geom_point(alpha = 0.5) +
  theme_classic()

ggsave(
  "results/Figure6_volcanoplotOAS3_High_vs_Low.png",
  width = 8,
  height = 6,
  dpi = 300
)

deg$Significance <- "Not Significant"

deg$Significance[
  deg$adj.P.Val < 0.05 &
    deg$logFC > 1
] <- "Upregulated"

deg$Significance[
  deg$adj.P.Val < 0.05 &
    deg$logFC < -1
] <- "Downregulated"

volcano_plot <- ggplot(deg,
                       aes(logFC,
                           -log10(adj.P.Val),
                           color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not Significant" = "grey70"
  )) +
  theme_classic() +
  labs(
    title = "Differentially Expressed Genes between OAS3-High and OAS3-Low Tumors",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value"
  )
volcano_plot

ggsave(
  filename = "figures/Figure13_Improved_Volcano_Plot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

head(deg)
deg_sig <- deg %>%
  filter(
    adj.P.Val < 0.05,
    abs(logFC) > 1
  )
nrow(deg_sig)
saveRDS(deg_sig,
        "data/processed/deg_sig.rds")
