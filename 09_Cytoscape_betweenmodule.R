# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 09 - Cytoscape network with between-module connections
# -----------------------------------------------------------------------------

library(WGCNA)

# 1. Load WGCNA results
load("results/WGCNA_network.RData")

# Load existing Cytoscape tables
nodes <- read.csv("results/cytoscape_nodes_top3.csv")
edges_within <- read.csv("results/cytoscape_edges_top3.csv")

# Check data
dim(nodes)
dim(edges_within)
table(nodes$WGCNA_Module)

# -----------------------------------------------------------------------------
# 2. Calculate TOM
# -----------------------------------------------------------------------------

adjacency_matrix <- adjacency(
  datMeth_var,
  power = softPower,
  type = "signed"
)

TOM <- TOMsimilarity(
  adjacency_matrix,
  TOMType = "signed"
)

rownames(TOM) <- colnames(datMeth_var)
colnames(TOM) <- colnames(datMeth_var)

# -----------------------------------------------------------------------------
# 3. Examine between-module TOM values
# -----------------------------------------------------------------------------

# Keep only genes assigned to detected modules
genes <- nodes$Gene

tom_selected <- TOM[genes, genes]

# Module assignment for each gene
module <- setNames(nodes$WGCNA_Module, nodes$Gene)

# Extract unique gene pairs
upper <- which(
  upper.tri(tom_selected),
  arr.ind = TRUE
)

between_edges <- data.frame(
  Source = rownames(tom_selected)[upper[, 1]],
  Target = colnames(tom_selected)[upper[, 2]],
  TOM = tom_selected[upper]
)

# Identify between-module edges
between_edges$SourceModule <- module[between_edges$Source]
between_edges$TargetModule <- module[between_edges$Target]

between_edges <- between_edges[
  between_edges$SourceModule != between_edges$TargetModule,
]

# Examine TOM values
summary(between_edges$TOM)

quantile(
  between_edges$TOM,
  c(0.50, 0.75, 0.90, 0.95, 0.99, 0.999)
)

# -----------------------------------------------------------------------------
# 4. Select strongest between-module connections
# -----------------------------------------------------------------------------

# Select the top 0.1% of between-module TOM values
threshold <- quantile(
  between_edges$TOM,
  0.999,
  na.rm = TRUE
)

edges_between <- between_edges[
  between_edges$TOM >= threshold,
]

# Mark edge types
edges_between$EdgeType <- "Between-module"
edges_within$EdgeType <- "Within-module"

# Keep the same columns in both edge tables
edges_between <- edges_between[
  , c("Source", "Target", "TOM", "EdgeType")
]

edges_within <- edges_within[
  , c("Source", "Target", "TOM", "EdgeType")
]

# Combine the two edge sets
edges_combined <- rbind(
  edges_within,
  edges_between
)

# Check results
nrow(edges_within)
nrow(edges_between)
nrow(edges_combined)

table(edges_combined$EdgeType)

# -----------------------------------------------------------------------------
# 5. Export combined network
# -----------------------------------------------------------------------------

write.csv(
  edges_combined,
  "results/cytoscape_edges_inter_module.csv",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  nodes,
  "results/cytoscape_nodes_inter_module.csv",
  row.names = FALSE,
  quote = FALSE
)

# Count selected between-module edges by module pair
table(
  paste(
    pmin(between_edges$SourceModule[
      between_edges$TOM >= threshold
    ],
    between_edges$TargetModule[
      between_edges$TOM >= threshold
    ]),
    pmax(between_edges$SourceModule[
      between_edges$TOM >= threshold
    ],
    between_edges$TargetModule[
      between_edges$TOM >= threshold
    ]),
    sep = " - "
  )
)
