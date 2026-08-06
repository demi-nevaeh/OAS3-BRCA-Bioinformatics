deg_sig <- readRDS("data/processed/deg_sig.rds")

gene_ids <- bitr(
  deg_sig$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
head(gene_ids)

dim(gene_ids)

go_deg <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(as.data.frame(go_deg))

kegg_deg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)
head(as.data.frame(kegg_deg))

write.csv(
  as.data.frame(go_deg),
  "results/GO_DEG_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(kegg_deg),
  "results/KEGG_DEG_results.csv",
  row.names = FALSE
)

dotplot(go_deg,
        showCategory = 15) +
  ggtitle("GO Biological Processes Enriched in OAS3-associated DEGs")
ggsave(
  "figures/Figure9_GO_Dotplot.png",
  width=8,
  height=6,
  dpi=300
)

kegg_plot <- dotplot(
  kegg_deg,
  showCategory = 15
) +
  ggtitle("KEGG Pathways Enriched in OAS3-associated DEGs")

kegg_plot

ggsave(
  "figures/Figure10_KEGG_Dotplot.png",
  width=8,
  height=6,
  dpi=300
)

barplot(go_deg, showCategory = 15)
ggsave(
  "figures/Figure11_GO_BARplot.png",
  width=8,
  height=6,
  dpi=300
)

barplot(kegg_deg, showCategory = 15)
ggsave(
  "figures/Figure12_KEGG_BARplot.png",
  width=8,
  height=6,
  dpi=300
)
