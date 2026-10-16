# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 02 - WGCNA network construction
# -----------------------------------------------------------------------------

library(WGCNA)

# ------------------------------------------------------------
# 1. Load preprocessed methylation data
# ------------------------------------------------------------

load("results/preprocessed_methylation.RData")

dim(datMeth_var)

# ------------------------------------------------------------
# 2. Select soft-thresholding power
# ------------------------------------------------------------

# Candidate soft-thresholding powers
powers <- c(1:10, seq(12, 20, by = 2))

# Test powers
sft <- pickSoftThreshold(
  datMeth_var,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)

png("figures/soft_threshold.png", width = 1400, height = 650)

par(mfrow = c(1, 2))

# Scale-free topology fit
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft-thresholding power",
  ylab = "Scale-free topology fit (signed R²)",
  type = "n",
  main = "Scale independence"
)

text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = 0.8
)

abline(h = 0.8, lty = 2)


# Mean connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft-thresholding power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = 0.8
)

dev.off()

par(mfrow = c(1, 1))


# Choose beta = 14, the lowest tested power reaching R² >= 0.8
# while maintaining reasonable mean connectivity
softPower <- 14

# ------------------------------------------------------------
# 3. Construct WGCNA network and detect modules
# ------------------------------------------------------------

net <- blockwiseModules(
  datMeth_var,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)

# Number of genes in each module
table(moduleColors)

# ------------------------------------------------------------
# 4. Visualize detected modules
# ------------------------------------------------------------

png("figures/WGCNA_modules.png", width = 1400, height = 800)

plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# ------------------------------------------------------------
# 5. Calculate module eigengenes
# ------------------------------------------------------------

MEs <- orderMEs(net$MEs)

dim(MEs)
colnames(MEs)

# ------------------------------------------------------------
# 6. Examine correlations between module eigengenes
# ------------------------------------------------------------

colnames(MEs)[colnames(MEs) == "ME0"] <- "Grey"
colnames(MEs)[colnames(MEs) == "ME1"] <- "Turquoise"
colnames(MEs)[colnames(MEs) == "ME2"] <- "Blue"
colnames(MEs)[colnames(MEs) == "ME3"] <- "Brown"
colnames(MEs)[colnames(MEs) == "ME4"] <- "Yellow"

# Exclude the grey eigengene (unassigned genes)
MEs_modules <- MEs[, colnames(MEs) != "Grey"]

MEcor_modules <- cor(
  MEs_modules,
  use = "pairwise.complete.obs"
)

round(MEcor_modules, 2)

png("figures/module_eigengene_correlations.png",
    width = 1000, height = 900)

par(mar = c(7, 7, 4, 2))

labeledHeatmap(
  Matrix = MEcor_modules,
  xLabels = colnames(MEs_modules),
  yLabels = colnames(MEs_modules),
  ySymbols = colnames(MEs_modules),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = round(MEcor_modules, 2),
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Co-methylation module correlations"
)

dev.off()

# ------------------------------------------------------------
# 7. Calculate module membership (kME)
# ------------------------------------------------------------

# Module eigengenes using the original numeric labels
MEs_numeric <- orderMEs(net$MEs)

# Correlation of every gene with every module eigengene
kME <- cor(
  datMeth_var,
  MEs_numeric,
  use = "pairwise.complete.obs"
)

dim(kME)
colnames(kME)

# ------------------------------------------------------------
# 8. Assign module membership to each gene
# ------------------------------------------------------------

# Create table containing gene and module
gene_info <- data.frame(
  Gene = colnames(datMeth_var),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

# Map module colors to their eigengene columns
module_to_ME <- c(
  turquoise = "ME1",
  blue = "ME2",
  brown = "ME3",
  yellow = "ME4",
  grey = "ME0"
)

# Get each gene's kME for its assigned module
gene_info$kME <- sapply(
  seq_len(nrow(gene_info)),
  function(i) {
    kME[i, module_to_ME[gene_info$Module[i]]]
  }
)

head(gene_info)
summary(gene_info$kME)

# ------------------------------------------------------------
# 9. Identify the top 10 hub genes per module
# ------------------------------------------------------------

module_genes <- gene_info[gene_info$Module != "grey", ]

top_hubs <- do.call(
  rbind,
  lapply(
    split(module_genes, module_genes$Module),
    function(x) {
      head(x[order(x$kME, decreasing = TRUE), ], 10)
    }
  )
)

# Select the 10 genes with the highest kME in each module
# High kME indicates strong correlation with the module eigengene

top_hubs

write.csv(
  top_hubs,
  "results/top10_hub_genes_per_module.csv",
  row.names = FALSE
)

write.csv(
  module_genes,
  "results/module_genes_kME.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Save results
# ------------------------------------------------------------

save(
  net,
  moduleColors,
  datMeth_var,
  softPower,
  MEs_modules,
  gene_info,
  file = "results/WGCNA_network.RData"
)
