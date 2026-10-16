# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 04 - WGCNA without variance filtering
# -----------------------------------------------------------------------------

library(WGCNA)

# ------------------------------------------------------------
# 1. Load methylation data after missingness QC
# ------------------------------------------------------------

load("results/preprocessed_methylation.RData")

# Use all genes remaining after missingness filtering
dim(datMeth)

# ------------------------------------------------------------
# 2. Check data quality
# ------------------------------------------------------------

gsg_all <- goodSamplesGenes(
  datMeth,
  verbose = 3
)

gsg_all$allOK

# ------------------------------------------------------------
# 3. Select soft-thresholding power
# ------------------------------------------------------------

# Test an extended range of powers
powers_extended <- c(1:10, seq(12, 30, by = 2))

sft_all_extended <- pickSoftThreshold(
  datMeth,
  powerVector = powers_extended,
  networkType = "signed",
  verbose = 5
)

sft_all_extended$fitIndices

png("figures/soft_threshold_all_genes.png",
    width = 1400, height = 650)

par(mfrow = c(1, 2))

# Scale-free topology fit
plot(
  sft_all_extended$fitIndices[, 1],
  -sign(sft_all_extended$fitIndices[, 3]) *
    sft_all_extended$fitIndices[, 2],
  xlab = "Soft-thresholding power",
  ylab = "Scale-free topology fit (signed R²)",
  type = "n",
  main = "Scale independence"
)

text(
  sft_all_extended$fitIndices[, 1],
  -sign(sft_all_extended$fitIndices[, 3]) *
    sft_all_extended$fitIndices[, 2],
  labels = powers_extended,
  cex = 0.8
)

abline(h = 0.8, lty = 2)

# Mean connectivity
plot(
  sft_all_extended$fitIndices[, 1],
  sft_all_extended$fitIndices[, 5],
  xlab = "Soft-thresholding power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  sft_all_extended$fitIndices[, 1],
  sft_all_extended$fitIndices[, 5],
  labels = powers_extended,
  cex = 0.8
)

dev.off()

par(mfrow = c(1, 1))

# Selected beta = 28 for the exploratory analysis
# Scale-free topology fit was below 0.8 for the tested powers
softPower_all <- 28

# ------------------------------------------------------------
# 4. Construct WGCNA network using all genes
# ------------------------------------------------------------

net_all <- blockwiseModules(
  datMeth,
  power = softPower_all,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)

moduleColors_all <- labels2colors(net_all$colors)

# ------------------------------------------------------------
# 5. Examine detected modules
# ------------------------------------------------------------

table(moduleColors_all)

# Number of genes in each block
sapply(net_all$blockGenes, length)

# ------------------------------------------------------------
# 6. Examine module distribution across blocks
# ------------------------------------------------------------

block_module_table <- sapply(
  seq_along(net_all$blockGenes),
  function(i) {
    table(
      factor(
        moduleColors_all[net_all$blockGenes[[i]]],
        levels = c("blue", "brown", "grey", "turquoise", "yellow")
      )
    )
  }
)

colnames(block_module_table) <- paste0(
  "Block_", seq_len(ncol(block_module_table))
)

block_module_table

# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

save(
  net_all,
  moduleColors_all,
  softPower_all,
  sft_all_extended,
  file = "results/WGCNA_all_genes.RData"
)
