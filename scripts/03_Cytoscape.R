# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 03 - Export WGCNA network to Cytoscape
# -----------------------------------------------------------------------------

library(WGCNA)

# ------------------------------------------------------------
# 1. Load WGCNA results
# ------------------------------------------------------------

load("results/WGCNA_network.RData")

dim(datMeth_var)
table(moduleColors)

# ------------------------------------------------------------
# 2. Calculate the Topological Overlap Matrix (TOM)
# ------------------------------------------------------------

# Calculate signed adjacency using the selected soft-thresholding power
adjacency_matrix <- adjacency(
  datMeth_var,
  power = softPower,
  type = "signed"
)

# Calculate TOM to measure similarity in network connections
TOM <- TOMsimilarity(
  adjacency_matrix,
  TOMType = "signed"
)

# Add gene names
rownames(TOM) <- colnames(datMeth_var)
colnames(TOM) <- colnames(datMeth_var)

# ------------------------------------------------------------
# 3. Prepare node attributes
# ------------------------------------------------------------

gene_info <- read.csv("results/module_genes_kME.csv")
top_hubs <- read.csv("results/top10_hub_genes_per_module.csv")

# Keep all genes assigned to a WGCNA module (exclude grey)
nodes_all <- gene_info[
  gene_info$Module != "grey",
]

# Identify the top 10 genes by kME in each module
nodes_all$TopHub <- nodes_all$Gene %in% top_hubs$Gene

dim(nodes_all)
table(nodes_all$Module)
table(nodes_all$TopHub)

# ------------------------------------------------------------
# 4. Examine within-module TOM connections
# ------------------------------------------------------------

# Calculate the strongest TOM connection for each gene
strongest_TOM <- numeric()

for (mod in c("blue", "brown", "turquoise", "yellow")) {
  
  genes <- colnames(datMeth_var)[moduleColors == mod]
  tom_mod <- TOM[genes, genes]
  
  # Exclude self-connections
  diag(tom_mod) <- NA
  
  strongest_TOM <- c(
    strongest_TOM,
    apply(tom_mod, 1, max, na.rm = TRUE)
  )
}

summary(strongest_TOM)

quantile(
  strongest_TOM,
  c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
)

# ------------------------------------------------------------
# 5. Select the top 3 connections per gene
# ------------------------------------------------------------

# For each gene, select its three strongest connections
# to other genes in the same WGCNA module

edge_list <- data.frame()

for (mod in c("blue", "brown", "turquoise", "yellow")) {
  
  genes <- colnames(datMeth_var)[moduleColors == mod]
  tom_mod <- TOM[genes, genes]
  
  # Exclude self-connections
  diag(tom_mod) <- NA
  
  for (gene in genes) {
    
    weights <- tom_mod[gene, ]
    
    # Select the three strongest TOM connections
    top_genes <- names(sort(weights, decreasing = TRUE))[1:3]
    
    temp <- data.frame(
      Source = gene,
      Target = top_genes,
      TOM = weights[top_genes],
      Module = mod
    )
    
    edge_list <- rbind(edge_list, temp)
  }
}

# ------------------------------------------------------------
# 6. Remove duplicated undirected edges
# ------------------------------------------------------------

# An interaction A-B is equivalent to B-A
# Keep each gene pair only once

edge_list$pair <- apply(
  edge_list[, c("Source", "Target")],
  1,
  function(x) paste(sort(x), collapse = "_")
)

edges_top3 <- edge_list[!duplicated(edge_list$pair), ]

edges_top3$pair <- NULL

# Check final network
nrow(edges_top3)
length(unique(c(edges_top3$Source, edges_top3$Target)))
table(edges_top3$Module)

# Check that every module-assigned gene appears in the network
setequal(
  nodes_all$Gene,
  unique(c(edges_top3$Source, edges_top3$Target))
)

# ------------------------------------------------------------
# 7. Export node and edge tables for Cytoscape
# ------------------------------------------------------------

# Rename node attributes for clarity in Cytoscape
nodes_cytoscape_all <- nodes_all

names(nodes_cytoscape_all)[
  names(nodes_cytoscape_all) == "Module"
] <- "WGCNA_Module"

names(nodes_cytoscape_all)[
  names(nodes_cytoscape_all) == "kME"
] <- "WGCNA_kME"

names(nodes_cytoscape_all)[
  names(nodes_cytoscape_all) == "TopHub"
] <- "WGCNA_TopHub"

head(nodes_cytoscape_all)

write.csv(
  edges_top3,
  "results/cytoscape_edges_top3.csv",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  nodes_cytoscape_all,
  "results/cytoscape_nodes_top3.csv",
  row.names = FALSE,
  quote = FALSE
)

