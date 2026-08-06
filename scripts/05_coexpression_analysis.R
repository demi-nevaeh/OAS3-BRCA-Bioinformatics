expr_annotated <- readRDS(
  "data/processed/expr_annotated.rds"
)

oas3_surv <- readRDS(
  "data/processed/oas3_surv.rds"
)

table(oas3_surv$OAS3_group)

which(expr_annotated$gene == "OAS3")

oas3_vector <- as.numeric(
  expr_numeric[expr_annotated$gene == "OAS3", ]
)

cor_values <- apply(
  expr_numeric,
  1,
  function(x)
    cor(
      x,
      oas3_vector,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
)


cor_results <- cor_results %>%
  arrange(desc(abs(correlation)))

head(cor_results, 20)

#preparing for GO/KEGG enrichment
#remove oas3 itself
cor_results_unique <- cor_results_unique %>%
  filter(gene != "OAS3")
#Select strongly positively correlated genes
top_genes <- cor_results_unique %>%
  filter(correlation >= 0.5)
nrow(top_genes)

#to perform GO Enrichment
BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Hs.eg.db)
gene_ids <- bitr(
  top_genes$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
#run GO enrichment
go_results <- enrichGO(
  gene          = gene_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
head(as.data.frame(go_results))
kegg_results <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(as.data.frame(kegg_results))

write.csv(
  as.data.frame(go_results),
  "results/OAS3_GO_BP.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(kegg_results),
  "results/OAS3_KEGG.csv",
  row.names = FALSE
)

library(clusterProfiler)

dotplot(go_results,
        showCategory = 10,
        font.size = 12,
        title = "GO Biological Process Enrichment of OAS3 Co-expressed Genes")
ggsave(
  "results/Figure4_GO_Dotplot.png",
  width = 8,
  height = 6,
  dpi = 300
)
dotplot(kegg_results,
        showCategory = 10,
        font.size = 12,
        title = "KEGG Pathway Enrichment of OAS3 Co-expressed Genes")
ggsave(
  "results/Figure5_KEGG_Dotplot.png",
  width = 8,
  height = 6,
  dpi = 300
)
install.packages("pheatmap")
library(pheatmap)

top30 <- cor_results %>%
  filter(gene != "OAS3") %>%
  arrange(desc(correlation))

top30 <- top30[1:30, ]

heat_data <- expr_annotated %>%
  filter(gene %in% top30$gene) %>%
  dplyr::select(-gene_id, -gene, -chrom, -chromStart, -chromEnd, -strand)

rownames(heat_data) <- top30$gene

pheatmap(
  as.matrix(heat_data),
  scale = "row",
  show_colnames = FALSE,
  fontsize_row = 8
)
write.csv(
  top30,
  "results/Top30_OAS3_Coexpressed_Genes_BRCA.csv",
  row.names = FALSE
)
saveRDS(
  top30,
  "results/Top30_OAS3_Coexpressed_Genes_BRCA.rds"
)


