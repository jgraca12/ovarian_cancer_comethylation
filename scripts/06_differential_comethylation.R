# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 06 - Differential co-methylation analysis
# -----------------------------------------------------------------------------

library(WGCNA)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1. Load platinum association results
# -----------------------------------------------------------------------------

load("results/platinum_association.RData")

# Modules identified in the WGCNA network
modules <- c("blue", "brown", "turquoise", "yellow")

# Check data
dim(meth_platinum)
table(status)
table(gene_results$Module)


# -----------------------------------------------------------------------------
# 2. Compare gene-gene correlations within modules
# -----------------------------------------------------------------------------

# Function to summarize correlations for one module
module_correlations <- function(module_name) {
  
  # Select genes belonging to the module
  genes <- gene_results$Gene[
    gene_results$Module == module_name
  ]
  
  # Extract methylation data
  data <- meth_platinum[, genes, drop = FALSE]
  
  # Separate platinum groups
  data_sensitive <- data[
    status == "Sensitive", , drop = FALSE
  ]
  
  data_resistant <- data[
    status == "Resistant", , drop = FALSE
  ]
  
  # Calculate Pearson correlations
  cor_sensitive <- cor(
    data_sensitive,
    use = "pairwise.complete.obs"
  )
  
  cor_resistant <- cor(
    data_resistant,
    use = "pairwise.complete.obs"
  )
  
  # Extract unique gene pairs (exclude diagonal)
  upper <- upper.tri(cor_sensitive)
  
  r_sensitive <- cor_sensitive[upper]
  r_resistant <- cor_resistant[upper]
  
  # Summarize median absolute correlations
  median_sensitive <- median(
    abs(r_sensitive), na.rm = TRUE
  )
  
  median_resistant <- median(
    abs(r_resistant), na.rm = TRUE
  )
  
  data.frame(
    Module = module_name,
    Genes = length(genes),
    MedianAbsCor_Sensitive = median_sensitive,
    MedianAbsCor_Resistant = median_resistant,
    Difference = median_resistant - median_sensitive
  )
}

# Apply to all four modules
cor_results <- do.call(
  rbind,
  lapply(modules, module_correlations)
)

cor_results


# -----------------------------------------------------------------------------
# 3. Permutation test for differences in co-methylation
# -----------------------------------------------------------------------------

set.seed(123)
n_perm <- 1000

# Function to calculate median absolute Pearson correlation
median_abs_cor <- function(data) {
  
  cor_matrix <- cor(
    data,
    use = "pairwise.complete.obs"
  )
  
  median(
    abs(cor_matrix[upper.tri(cor_matrix)]),
    na.rm = TRUE
  )
}

# Calculate correlation differences for each module
calculate_differences <- function(group_labels) {
  
  sapply(modules, function(mod) {
    
    genes <- gene_results$Gene[
      gene_results$Module == mod
    ]
    
    data <- meth_platinum[, genes, drop = FALSE]
    
    sensitive_data <- data[
      group_labels == "Sensitive", , drop = FALSE
    ]
    
    resistant_data <- data[
      group_labels == "Resistant", , drop = FALSE
    ]
    
    median_abs_cor(resistant_data) -
      median_abs_cor(sensitive_data)
  })
}

# Observed differences
observed_diff <- calculate_differences(status)

# Store permutation results
perm_diff <- matrix(
  NA_real_,
  nrow = n_perm,
  ncol = length(modules),
  dimnames = list(NULL, modules)
)

# Shuffle platinum labels
for (b in seq_len(n_perm)) {
  
  shuffled_status <- sample(status)
  
  perm_diff[b, ] <- calculate_differences(
    shuffled_status
  )
  
  if (b %% 100 == 0) {
    cat("Completed", b, "permutations\n")
  }
}

# Two-sided permutation p-value for turquoise
p_turquoise_cor <- (
  1 + sum(
    abs(perm_diff[, "turquoise"]) >=
      abs(observed_diff["turquoise"])
  )
) / (n_perm + 1)

# Correct for selecting the module with the largest difference
p_max_cor <- (
  1 + sum(
    apply(abs(perm_diff), 1, max) >=
      max(abs(observed_diff))
  )
) / (n_perm + 1)

round(observed_diff, 4)
p_turquoise_cor
p_max_cor


# -----------------------------------------------------------------------------
# 4. Compare correlations in the turquoise module
# -----------------------------------------------------------------------------

# Select turquoise genes
turquoise_genes <- gene_results$Gene[
  gene_results$Module == "turquoise"
]

turquoise_data <- meth_platinum[
  , turquoise_genes, drop = FALSE
]

# Calculate correlations for each platinum group
cor_sens <- cor(
  turquoise_data[status == "Sensitive", , drop = FALSE],
  use = "pairwise.complete.obs"
)

cor_res <- cor(
  turquoise_data[status == "Resistant", , drop = FALSE],
  use = "pairwise.complete.obs"
)

# Extract unique gene pairs
upper <- upper.tri(cor_sens)

cor_df <- data.frame(
  Correlation = c(
    cor_sens[upper],
    cor_res[upper]
  ),
  PlatinumStatus = rep(
    c("Sensitive", "Resistant"),
    each = sum(upper)
  )
)

# Plot correlation distributions
p <- ggplot(
  cor_df,
  aes(x = Correlation, fill = PlatinumStatus)
) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Turquoise module: co-methylation correlations",
    x = "Gene-gene Pearson correlation",
    y = "Density"
  ) +
  theme_bw()

print(p)

ggsave(
  "figures/turquoise_correlation_distributions.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Boxplot of absolute correlations
# -----------------------------------------------------------------------------

# Convert signed correlations to absolute correlations
cor_df$AbsCorrelation <- abs(cor_df$Correlation)

# Plot distributions by platinum status
p_abs <- ggplot(
  cor_df,
  aes(x = PlatinumStatus,
      y = AbsCorrelation,
      fill = PlatinumStatus)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  labs(
    title = "Turquoise module: co-methylation strength",
    x = "Platinum status",
    y = "Absolute Pearson correlation"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

print(p_abs)

ggsave(
  "figures/turquoise_absolute_correlations_boxplot.png",
  plot = p_abs,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Compare median correlation strength across all modules
# -----------------------------------------------------------------------------

# Convert results into long format for plotting
plot_data <- rbind(
  data.frame(
    Module = cor_results$Module,
    PlatinumStatus = "Sensitive",
    MedianAbsCor = cor_results$MedianAbsCor_Sensitive
  ),
  data.frame(
    Module = cor_results$Module,
    PlatinumStatus = "Resistant",
    MedianAbsCor = cor_results$MedianAbsCor_Resistant
  )
)

# Plot median absolute correlations for each module
p_modules <- ggplot(
  plot_data,
  aes(
    x = Module,
    y = MedianAbsCor,
    fill = PlatinumStatus
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(
    title = "Co-methylation strength across WGCNA modules",
    x = "Module",
    y = "Median absolute Pearson correlation",
    fill = "Platinum status"
  ) +
  theme_bw()

ggsave(
  "figures/module_platinum_correlation_comparison.png",
  plot = p_modules,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_modules)

# -----------------------------------------------------------------------------
# 5. Check the effect of unequal group sizes
# -----------------------------------------------------------------------------

set.seed(123)
n_repeats <- 500

# Separate turquoise methylation data
sens_data <- turquoise_data[
  status == "Sensitive", , drop = FALSE
]

res_data <- turquoise_data[
  status == "Resistant", , drop = FALSE
]

# Median absolute correlation in resistant patients
res_cor <- median_abs_cor(res_data)

# Repeatedly sample sensitive patients to match resistant group size
sens_subsample_cor <- replicate(n_repeats, {
  
  selected <- sample(
    seq_len(nrow(sens_data)),
    size = nrow(res_data),
    replace = FALSE
  )
  
  median_abs_cor(
    sens_data[selected, , drop = FALSE]
  )
})

# Compare results
summary(sens_subsample_cor)
res_cor

# Fraction of sensitive subsamples with weaker correlations
mean(sens_subsample_cor < res_cor)


# -----------------------------------------------------------------------------
# 6. Check correlations using Spearman instead of Pearson
# -----------------------------------------------------------------------------

# Function to calculate median absolute Spearman correlation
median_abs_spearman <- function(data) {
  
  cor_matrix <- cor(
    data,
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  
  median(
    abs(cor_matrix[upper.tri(cor_matrix)]),
    na.rm = TRUE
  )
}

spearman_sens <- median_abs_spearman(
  turquoise_data[status == "Sensitive", , drop = FALSE]
)

spearman_res <- median_abs_spearman(
  turquoise_data[status == "Resistant", , drop = FALSE]
)

spearman_results <- data.frame(
  Sensitive = spearman_sens,
  Resistant = spearman_res,
  Difference = spearman_res - spearman_sens
)

spearman_results


# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------

write.csv(
  cor_results,
  "results/differential_comethylation_modules.csv",
  row.names = FALSE
)

write.csv(
  spearman_results,
  "results/turquoise_spearman_comparison.csv",
  row.names = FALSE
)

save(
  cor_results,
  observed_diff,
  perm_diff,
  p_turquoise_cor,
  p_max_cor,
  turquoise_data,
  cor_sens,
  cor_res,
  sens_subsample_cor,
  res_cor,
  spearman_results,
  file = "results/differential_comethylation.RData"
)
