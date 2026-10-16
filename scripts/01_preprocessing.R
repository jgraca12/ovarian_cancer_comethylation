# -----------------------------------------------------------------------------
# DNA methylation co-methylation network - Ovarian cancer
# 01 - Data preprocessing
# -----------------------------------------------------------------------------
library(WGCNA)

# Create output folders if they do not already exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Load and inspect methylation data
# -----------------------------------------------------------------------------

# rows = genes, columns = tumor samples
# first two columns = gene identifiers
methylation <- read.delim("data/data_methylation_hm27.txt", header = TRUE, check.names = FALSE)

# Dataset dimensions
dim(methylation)
head(methylation[,1:6])
colnames(methylation)[1:10]

#Number of genes = 9452
nrow(methylation)

#Number of sample = 489
ncol(methylation) - 2

#Number of missing methylation values
sum(is.na(methylation[, -(1:2)]))

#Number of genes with at least 1 missing value
sum(rowSums(is.na(methylation[, -(1:2)])) > 0)

#Number of samples containing at least one missing value
sum(colSums(is.na(methylation[, -(1:2)])) > 0)

#Check duplicated gene symbols
sum(duplicated(methylation$Hugo_Symbol))

#Check missing gene symbols
sum(is.na(methylation$Hugo_Symbol) | methylation$Hugo_Symbol == "")

# Check that methylation measurements are numeric
all(sapply(methylation[, -(1:2)], is.numeric))


# ------------------------------------------------------------
# 2. Investigate missing values
# ------------------------------------------------------------

# Number of missing values per gene
missing_per_gene <- rowSums(is.na(methylation[, -(1:2)]))

summary(missing_per_gene)

# Percentage of missing values per gene
missing_percent_gene <- missing_per_gene / (ncol(methylation) - 2) * 100

summary(missing_percent_gene)

# Number of genes above different missingness thresholds
sum(missing_percent_gene > 5)
sum(missing_percent_gene > 10)
sum(missing_percent_gene > 20)

# ------------------------------------------------------------
# 3. Missing values per sample
# ------------------------------------------------------------

missing_per_sample <- colSums(is.na(methylation[, -(1:2)]))

missing_percent_sample <- 
  missing_per_sample / nrow(methylation) * 100

summary(missing_percent_sample)

# Samples above different missingness thresholds
sum(missing_percent_sample > 5)  #chosen missingness cutoff
sum(missing_percent_sample > 10)
sum(missing_percent_sample > 20)

# Worst affected samples
sort(missing_percent_sample, decreasing = TRUE)[1:10]

# ------------------------------------------------------------
# 4. Filter genes and samples with >5% missing values
# ------------------------------------------------------------

# Keep genes with <= 5% missing values
keep_genes <- missing_percent_gene <= 5

methylation_filtered <- methylation[keep_genes, ]

# Recalculate sample missingness after gene filtering
sample_missing_filtered <- colMeans(
  is.na(methylation_filtered[, -(1:2)])
) * 100

# Keep samples with <= 5% missing values
keep_samples <- sample_missing_filtered <= 5

# Keep the two identifier columns + acceptable samples
methylation_filtered <- methylation_filtered[
  ,
  c(TRUE, TRUE, keep_samples)
]

# Check dimensions and remaining missing values after filtering
dim(methylation_filtered)
sum(is.na(methylation_filtered[, -(1:2)]))

# Genes that still contain at least one NA
sum(rowSums(is.na(methylation_filtered[, -(1:2)])) > 0)

# Maximum number of missing values in a remaining gene
max(rowSums(is.na(methylation_filtered[, -(1:2)])))


# ------------------------------------------------------------
# 5. Prepare data for WGCNA
# ------------------------------------------------------------

# Extract methylation values
meth_matrix <- as.matrix(methylation_filtered[, -(1:2)])

# Add gene symbols
rownames(meth_matrix) <- methylation_filtered$Hugo_Symbol

# Transpose: samples become rows, genes become columns
datMeth <- as.data.frame(t(meth_matrix))

dim(datMeth)

gsg <- goodSamplesGenes(datMeth, verbose = 3)

gsg$allOK # TRUE: all genes and samples pass WGCNA quality checks, despite some missing values

if (!gsg$allOK) {
  stop("Some genes or samples did not pass WGCNA quality checks.")
}

# ------------------------------------------------------------
# 6. Calculate methylation variance for each gene
# ------------------------------------------------------------

gene_variance <- apply(datMeth, 2, var, na.rm = TRUE)

summary(gene_variance)

# Check genes with zero variance
sum(gene_variance == 0)

# Look at the most variable genes
head(sort(gene_variance, decreasing = TRUE), 10)

png("figures/methylation_variance.png", width = 1200, height = 800)

hist(
  gene_variance,
  breaks = 50,
  main = "Distribution of gene methylation variance",
  xlab = "Variance"
)

dev.off()

quantile(
  gene_variance,
  probs = c(0.5, 0.75, 0.8, 0.9, 0.95, 0.99),
  na.rm = TRUE
)

# ------------------------------------------------------------
# 7. Select top 25% most variable genes
# ------------------------------------------------------------

#Focus on genes with the greatest methylation variability across samples
#and reduce the number of genes used for network construction

# Variance threshold corresponding to the 75th percentile
variance_cutoff <- quantile(gene_variance, 0.75, na.rm = TRUE)

variance_cutoff

# Keep genes above the 75th percentile
keep_variable <- gene_variance >= variance_cutoff

datMeth_var <- datMeth[, keep_variable]

dim(datMeth_var)

summary(apply(datMeth_var, 2, var, na.rm = TRUE))

# ------------------------------------------------------------
# 8. Sample clustering to detect potential outliers
# ------------------------------------------------------------

# Cluster samples based on their methylation profiles to identify potential outliers
sampleTree <- hclust(
  dist(datMeth_var),
  method = "average"
)

png("figures/sample_clustering.png", width = 1200, height = 800)

plot(
  sampleTree,
  labels = FALSE,
  main = "Sample clustering to detect outliers",
  xlab = "Samples",
  sub = "",
  ylab = "Height"
)

dev.off()

# ------------------------------------------------------------
# 9. Save results
# ------------------------------------------------------------
save(
  datMeth,
  datMeth_var,
  file = "results/preprocessed_methylation.RData"
)
