# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 07 - Gene-specific co-methylation changes and GPR97 analysis
# -----------------------------------------------------------------------------

library(WGCNA)

# -----------------------------------------------------------------------------
# 1. Load previous results
# -----------------------------------------------------------------------------

load("results/platinum_association.RData")
load("results/differential_comethylation.RData")

# Check data
dim(turquoise_data)
table(status)


# -----------------------------------------------------------------------------
# 2. Identify gene-gene connections with the largest changes
# -----------------------------------------------------------------------------

# Extract unique gene pairs (exclude diagonal)
upper <- upper.tri(cor_sens)

# Create dataframe of gene-pair correlations
edge_changes <- data.frame(
  Gene1 = rownames(cor_sens)[row(cor_sens)[upper]],
  Gene2 = colnames(cor_sens)[col(cor_sens)[upper]],
  Cor_Sensitive = cor_sens[upper],
  Cor_Resistant = cor_res[upper]
)

# Calculate correlation differences
edge_changes$DeltaR <- with(
  edge_changes,
  Cor_Resistant - Cor_Sensitive
)

edge_changes$AbsDeltaR <- abs(edge_changes$DeltaR)

# Sort by largest absolute change
edge_changes <- edge_changes[
  order(-edge_changes$AbsDeltaR),
]

# Examine the largest changes
head(edge_changes, 10)
summary(edge_changes$DeltaR)

# Proportion of gene pairs with increased signed correlation
mean(edge_changes$DeltaR > 0)

# Number of pairs with large changes
sum(edge_changes$AbsDeltaR > 0.3)


# -----------------------------------------------------------------------------
# 3. Identify genes involved in large correlation changes
# -----------------------------------------------------------------------------

# Select connections with absolute correlation changes > 0.3
large_changes <- edge_changes[
  edge_changes$AbsDeltaR > 0.3,
]

# Count how often each gene appears
gene_counts <- table(c(
  large_changes$Gene1,
  large_changes$Gene2
))

# Sort from highest to lowest
gene_counts <- sort(gene_counts, decreasing = TRUE)

# Top 15 genes
head(gene_counts, 15)

# Load WGCNA hub information
hub_data <- read.csv(
  "results/top10_hub_genes_per_module.csv"
)

# Select turquoise hubs
turquoise_hubs <- hub_data[
  hub_data$Module == "turquoise",
]

# Compare genes with large changes to WGCNA hubs
top15_genes <- names(head(gene_counts, 15))

intersect(
  top15_genes,
  turquoise_hubs$Gene
)

comparison <- data.frame(
  Gene = top15_genes,
  LargeChangeConnections = as.integer(head(gene_counts, 15))
)

comparison$WGCNA_Top10Hub <- comparison$Gene %in%
  turquoise_hubs$Gene

comparison


# -----------------------------------------------------------------------------
# 4. Compare gene-level mean correlations
# -----------------------------------------------------------------------------

# Exclude self-correlations
diag(cor_sens) <- NA
diag(cor_res) <- NA

# Mean signed correlation with other turquoise genes
connectivity_results <- data.frame(
  Gene = rownames(cor_sens),
  MeanCor_Sensitive = rowMeans(cor_sens, na.rm = TRUE),
  MeanCor_Resistant = rowMeans(cor_res, na.rm = TRUE)
)

# Difference between platinum groups
connectivity_results$DeltaConnectivity <-
  connectivity_results$MeanCor_Resistant -
  connectivity_results$MeanCor_Sensitive

# Sort by largest increase
connectivity_results <- connectivity_results[
  order(-connectivity_results$DeltaConnectivity),
]

head(connectivity_results, 15)


# -----------------------------------------------------------------------------
# 5. Permutation test for the largest gene-level increase
# -----------------------------------------------------------------------------

set.seed(123)
n_perm <- 1000

# Function to calculate mean signed correlations
mean_connectivity <- function(data) {
  
  cor_mat <- cor(
    data,
    use = "pairwise.complete.obs"
  )
  
  diag(cor_mat) <- NA
  
  rowMeans(cor_mat, na.rm = TRUE)
}

# Function to calculate differences between groups
connectivity_diff <- function(group_labels) {
  
  sens <- turquoise_data[
    group_labels == "Sensitive", , drop = FALSE
  ]
  
  res <- turquoise_data[
    group_labels == "Resistant", , drop = FALSE
  ]
  
  mean_connectivity(res) - mean_connectivity(sens)
}

# Observed gene-level differences
observed_gene_diff <- connectivity_diff(status)

# Largest observed increase
observed_max <- max(observed_gene_diff, na.rm = TRUE)

# Store permutation maxima
perm_max <- numeric(n_perm)

# Shuffle platinum labels
for (b in seq_len(n_perm)) {
  
  shuffled_status <- sample(status)
  
  perm_gene_diff <- connectivity_diff(shuffled_status)
  
  perm_max[b] <- max(perm_gene_diff, na.rm = TRUE)
  
  if (b %% 100 == 0) {
    cat("Completed", b, "permutations\n")
  }
}

# Permutation p-value for the largest increase
p_gene_max <- (
  1 + sum(perm_max >= observed_max)
) / (n_perm + 1)

observed_max
p_gene_max


# -----------------------------------------------------------------------------
# 6. Check GPR97 using Spearman correlations
# -----------------------------------------------------------------------------

# Function to calculate mean signed Spearman correlations
mean_connectivity_spearman <- function(data) {
  
  cor_mat <- cor(
    data,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  diag(cor_mat) <- NA
  
  rowMeans(cor_mat, na.rm = TRUE)
}

sens_spearman <- mean_connectivity_spearman(
  turquoise_data[status == "Sensitive", , drop = FALSE]
)

res_spearman <- mean_connectivity_spearman(
  turquoise_data[status == "Resistant", , drop = FALSE]
)

# Examine GPR97
gpr97_spearman <- c(
  Sensitive = sens_spearman["GPR97"],
  Resistant = res_spearman["GPR97"],
  Difference = res_spearman["GPR97"] -
    sens_spearman["GPR97"]
)

gpr97_spearman


# -----------------------------------------------------------------------------
# 7. Select GPR97 connections for Cytoscape
# -----------------------------------------------------------------------------

# Select connections involving GPR97
gpr97_edges <- edge_changes[
  edge_changes$Gene1 == "GPR97" |
    edge_changes$Gene2 == "GPR97",
]

# Select the 15 largest absolute correlation changes
gpr97_top15 <- head(
  gpr97_edges[order(-gpr97_edges$AbsDeltaR), ],
  15
)

# Identify the other gene in each pair
gpr97_top15$Partner <- ifelse(
  gpr97_top15$Gene1 == "GPR97",
  gpr97_top15$Gene2,
  gpr97_top15$Gene1
)

# Inspect selected connections
gpr97_top15[
  , c("Partner", "Cor_Sensitive", "Cor_Resistant", "DeltaR")
]


# -----------------------------------------------------------------------------
# 8. Export networks for Cytoscape
# -----------------------------------------------------------------------------

# Platinum-sensitive network
edges_sensitive <- data.frame(
  source = "GPR97",
  target = gpr97_top15$Partner,
  Correlation = gpr97_top15$Cor_Sensitive,
  AbsCorrelation = abs(gpr97_top15$Cor_Sensitive),
  DeltaR = gpr97_top15$DeltaR
)

# Platinum-resistant network
edges_resistant <- data.frame(
  source = "GPR97",
  target = gpr97_top15$Partner,
  Correlation = gpr97_top15$Cor_Resistant,
  AbsCorrelation = abs(gpr97_top15$Cor_Resistant),
  DeltaR = gpr97_top15$DeltaR
)

write.csv(
  edges_sensitive,
  "results/GPR97_sensitive_edges.csv",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  edges_resistant,
  "results/GPR97_resistant_edges.csv",
  row.names = FALSE,
  quote = FALSE
)

# Combined edge table for comparing both conditions
edges_combined <- data.frame(
  source = "GPR97",
  target = gpr97_top15$Partner,
  Sensitive = abs(gpr97_top15$Cor_Sensitive),
  Resistant = abs(gpr97_top15$Cor_Resistant),
  DeltaR = gpr97_top15$DeltaR
)

write.csv(
  edges_combined,
  "results/GPR97_comparison_edges.csv",
  row.names = FALSE,
  quote = FALSE
)


# -----------------------------------------------------------------------------
# 9. Save results
# -----------------------------------------------------------------------------

write.csv(
  connectivity_results,
  "results/GPR97_gene_connectivity_changes.csv",
  row.names = FALSE
)

write.csv(
  gpr97_top15,
  "results/GPR97_top15_edges.csv",
  row.names = FALSE
)

save(
  edge_changes,
  connectivity_results,
  observed_gene_diff,
  observed_max,
  perm_max,
  p_gene_max,
  gpr97_spearman,
  gpr97_top15,
  file = "results/GPR97_analysis.RData"
)
