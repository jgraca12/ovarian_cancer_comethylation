# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 08 - GO enrichment analysis
# -----------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 1. Load WGCNA results
load("results/platinum_association.RData")
load("results/WGCNA_network.RData")

# Genes in the turquoise module
turquoise_genes <- unique(
  gene_results$Gene[
    gene_results$Module == "turquoise"
  ]
)

# Background: all genes included in main WGCNA analysis
background_genes <- colnames(datMeth_var)

stopifnot(
  length(turquoise_genes) == 363,
  length(background_genes) == 2339,
  all(turquoise_genes %in% background_genes),
  anyDuplicated(background_genes) == 0
)

# -----------------------------------------------------------------------------
# 2. GO Biological Process enrichment
# -----------------------------------------------------------------------------

go_results <- enrichGO(
  gene = turquoise_genes,
  universe = background_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# View enrichment results
go_table <- as.data.frame(go_results)

dim(go_table)
head(go_table)

nrow(go_table)

go_table[, c(
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust"
)] |> head(10)

# -----------------------------------------------------------------------------
# 3. Visualize GO enrichment
# -----------------------------------------------------------------------------

p_go <- dotplot(
  go_results,
  showCategory = 15
) +
  ggtitle("GO enrichment: turquoise module")

print(p_go)

ggsave(
  "figures/turquoise_GO_enrichment.png",
  plot = p_go,
  width = 10,
  height = 7,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 4. Simplify redundant GO terms
# -----------------------------------------------------------------------------

go_simplified <- simplify(
  go_results,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

# Compare number of enriched terms
nrow(as.data.frame(go_results))
nrow(as.data.frame(go_simplified))

# -----------------------------------------------------------------------------
# 5. Final GO enrichment plot
# -----------------------------------------------------------------------------

p_go_simplified <- dotplot(
  go_simplified,
  showCategory = 19,
  label_format = 40
) +
  ggtitle("GO enrichment: turquoise module") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 13)
  )

print(p_go_simplified)

ggsave(
  "figures/turquoise_GO_enrichment_simplified.png",
  plot = p_go_simplified,
  width = 12,
  height = 9,
  dpi = 300
)

# -----------------------------------------------------------------------------
# 6. Save GO enrichment results
# -----------------------------------------------------------------------------

write.csv(
  as.data.frame(go_results),
  "results/turquoise_GO_enrichment.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(go_simplified),
  "results/turquoise_GO_enrichment_simplified.csv",
  row.names = FALSE
)

save(
  go_results,
  go_simplified,
  turquoise_genes,
  background_genes,
  file = "results/GO_enrichment.RData"
)
