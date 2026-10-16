# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 05 - Association with platinum status
# -----------------------------------------------------------------------------

library(WGCNA)
library(ggplot2)
library(reshape2)

# -----------------------------------------------------------------------------
# 1. Load WGCNA and clinical data
# -----------------------------------------------------------------------------

load("results/WGCNA_network.RData")

clinical <- read.delim(
  "data/data_clinical_patient.txt",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Module eigengenes
MEs <- as.data.frame(orderMEs(net$MEs))

# Add sample and patient IDs
MEs$SampleID <- rownames(MEs)
MEs$PatientID <- substr(MEs$SampleID, 1, 12)

# Check patient IDs
sum(MEs$PatientID %in% clinical$`#Patient Identifier`)
anyDuplicated(MEs$PatientID)
anyDuplicated(clinical$`#Patient Identifier`)

# Examine platinum status categories
table(clinical$`Platinum Status`, useNA = "ifany")


# -----------------------------------------------------------------------------
# 2. Match methylation samples with platinum status
# -----------------------------------------------------------------------------

# Keep patients classified as sensitive or resistant
clinical_platinum <- clinical[
  clinical$`Platinum Status` %in% c("Sensitive", "Resistant"),
  c("#Patient Identifier", "Platinum Status")
]

colnames(clinical_platinum) <- c("PatientID", "PlatinumStatus")

# Check that each patient appears only once
stopifnot(
  !anyDuplicated(MEs$PatientID),
  !anyDuplicated(clinical_platinum$PatientID)
)

# Match clinical information to module eigengenes
ME_platinum <- merge(
  MEs,
  clinical_platinum,
  by = "PatientID"
)

# Rename modules
names(ME_platinum)[names(ME_platinum) == "ME1"] <- "Turquoise"
names(ME_platinum)[names(ME_platinum) == "ME2"] <- "Blue"
names(ME_platinum)[names(ME_platinum) == "ME3"] <- "Brown"
names(ME_platinum)[names(ME_platinum) == "ME4"] <- "Yellow"
names(ME_platinum)[names(ME_platinum) == "ME0"] <- "Grey"

# Check matched patients
dim(ME_platinum)
table(ME_platinum$PlatinumStatus)


# -----------------------------------------------------------------------------
# 3. Compare module eigengenes between platinum groups
# -----------------------------------------------------------------------------

modules <- c("Blue", "Brown", "Turquoise", "Yellow")

# Convert data to long format for plotting
ME_long <- melt(
  ME_platinum,
  id.vars = c("PatientID", "PlatinumStatus"),
  measure.vars = modules,
  variable.name = "Module",
  value.name = "Eigengene"
)

# Boxplots
p <- ggplot(
  ME_long,
  aes(x = PlatinumStatus, y = Eigengene, fill = PlatinumStatus)
) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.4) +
  facet_wrap(~ Module, scales = "free_y") +
  labs(
    x = "Platinum status",
    y = "Module eigengene"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(p)

ggsave(
  "figures/platinum_module_eigengenes.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 4. Test module eigengene differences
# -----------------------------------------------------------------------------

module_results <- data.frame(
  Module = modules,
  P_value = NA_real_
)

for (i in seq_along(modules)) {
  
  mod <- modules[i]
  
  test <- wilcox.test(
    ME_platinum[[mod]] ~ ME_platinum$PlatinumStatus,
    exact = FALSE
  )
  
  module_results$P_value[i] <- test$p.value
}

# Correct for testing four modules
module_results$FDR <- p.adjust(
  module_results$P_value,
  method = "BH"
)

module_results

# Compare median eigengene values
module_medians <- aggregate(
  ME_platinum[, modules],
  by = list(PlatinumStatus = ME_platinum$PlatinumStatus),
  FUN = median
)

module_medians


# -----------------------------------------------------------------------------
# 5. Match gene methylation profiles with platinum status
# -----------------------------------------------------------------------------

# Select genes assigned to WGCNA modules
module_genes <- colnames(datMeth_var)[moduleColors != "grey"]

length(module_genes)

# Extract methylation profiles
meth_modules <- datMeth_var[, module_genes, drop = FALSE]

# Convert sample IDs to patient IDs
patient_ids <- substr(rownames(meth_modules), 1, 12)

# Match patients in the same order as ME_platinum
match_idx <- match(ME_platinum$PatientID, patient_ids)

# Check for missing matches
sum(is.na(match_idx))
stopifnot(!anyNA(match_idx))

# Extract matched methylation data
meth_platinum <- meth_modules[match_idx, , drop = FALSE]

# Platinum status for the matched patients
status <- ME_platinum$PlatinumStatus

dim(meth_platinum)
table(status)


# -----------------------------------------------------------------------------
# 6. Gene-level differential methylation analysis
# -----------------------------------------------------------------------------

resistant <- status == "Resistant"
sensitive <- status == "Sensitive"

# Function to compare methylation between groups
test_gene <- function(x) {
  
  x_res <- x[resistant]
  x_sens <- x[sensitive]
  
  # Remove missing values
  x_res <- x_res[!is.na(x_res)]
  x_sens <- x_sens[!is.na(x_sens)]
  
  # Difference in median methylation
  delta_beta <- median(x_res) - median(x_sens)
  
  # Wilcoxon rank-sum test
  p_value <- wilcox.test(
    x_res,
    x_sens,
    exact = FALSE
  )$p.value
  
  c(
    DeltaBeta = delta_beta,
    P_value = p_value,
    N_Resistant = length(x_res),
    N_Sensitive = length(x_sens)
  )
}

# Apply to all module genes
gene_tests <- t(apply(meth_platinum, 2, test_gene))

gene_results <- data.frame(
  Gene = rownames(gene_tests),
  gene_tests,
  row.names = NULL
)

# Correct for multiple testing
gene_results$FDR <- p.adjust(
  gene_results$P_value,
  method = "BH"
)

# Add WGCNA module assignments
gene_results$Module <- moduleColors[
  match(gene_results$Gene, colnames(datMeth_var))
]

# Sort by adjusted p-value
gene_results <- gene_results[
  order(gene_results$FDR, gene_results$P_value),
]

head(gene_results, 10)

# Number of FDR-significant genes
sum(gene_results$FDR < 0.05)

# Significant genes per module
table(gene_results$Module[gene_results$FDR < 0.05])


# -----------------------------------------------------------------------------
# 7. Examine nominal associations across modules
# -----------------------------------------------------------------------------

# Number of genes with nominal p < 0.05 per module
table(
  Module = gene_results$Module,
  Nominal_P_less_0.05 = gene_results$P_value < 0.05
)

# Percentage of nominally associated genes per module
prop.table(
  table(
    gene_results$Module,
    gene_results$P_value < 0.05
  ),
  margin = 1
) * 100

# Compare turquoise with all other modules
enrichment_table <- table(
  Turquoise = gene_results$Module == "turquoise",
  Nominal = gene_results$P_value < 0.05
)

enrichment_table

# Fisher's exact test (exploratory)
fisher.test(enrichment_table)


# -----------------------------------------------------------------------------
# 8. Permutation test for module enrichment
# -----------------------------------------------------------------------------

set.seed(123)
n_perm <- 1000

# Keep genes in the same order as gene_results
meth_test <- meth_platinum[, gene_results$Gene, drop = FALSE]
gene_modules <- gene_results$Module

# Difference in nominal association rate:
# one module versus all other modules
enrichment_score <- function(p_values, modules) {
  
  associated <- p_values < 0.05
  
  sapply(unique(modules), function(mod) {
    mean(associated[modules == mod]) -
      mean(associated[modules != mod])
  })
}

# Observed enrichment scores
observed_score <- enrichment_score(
  gene_results$P_value,
  gene_modules
)

# Store permutation results
perm_scores <- matrix(
  NA_real_,
  nrow = n_perm,
  ncol = length(observed_score),
  dimnames = list(NULL, names(observed_score))
)

# Shuffle platinum labels while keeping methylation data unchanged
for (b in seq_len(n_perm)) {
  
  shuffled_status <- sample(status)
  
  perm_p <- apply(meth_test, 2, function(x) {
    
    x_res <- x[shuffled_status == "Resistant"]
    x_sens <- x[shuffled_status == "Sensitive"]
    
    # Use the same missing-value handling as the original tests
    x_res <- x_res[!is.na(x_res)]
    x_sens <- x_sens[!is.na(x_sens)]
    
    wilcox.test(
      x_res,
      x_sens,
      exact = FALSE
    )$p.value
  })
  
  perm_scores[b, ] <- enrichment_score(
    perm_p,
    gene_modules
  )
  
  if (b %% 100 == 0) {
    cat("Completed", b, "permutations\n")
  }
}

# Turquoise-specific permutation p-value
p_turquoise <- (
  1 + sum(
    perm_scores[, "turquoise"] >= observed_score["turquoise"]
  )
) / (n_perm + 1)

# Account for selecting the strongest module
p_max <- (
  1 + sum(
    apply(perm_scores, 1, max) >= max(observed_score)
  )
) / (n_perm + 1)

round(observed_score, 4)
p_turquoise
p_max


# -----------------------------------------------------------------------------
# 9. Save results
# -----------------------------------------------------------------------------

write.csv(
  module_results,
  "results/platinum_module_results.csv",
  row.names = FALSE
)

write.csv(
  gene_results,
  "results/platinum_gene_results.csv",
  row.names = FALSE
)

save(
  ME_platinum,
  meth_platinum,
  status,
  gene_results,
  module_results,
  module_medians,
  observed_score,
  perm_scores,
  p_turquoise,
  p_max,
  file = "results/platinum_association.RData"
)
