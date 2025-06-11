# EdgeR PCA Workflow: All Samples for PCA
#Complied by Lisa Ryno, Isabella Moppel, and Jason Kuchtey June 2025

# install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# change the library path as needed, to the first value returned by typing
#   .libPaths()
#   in the console
my_lib_path = .libPaths()[1]

# Load necessary libraries
library(tximport)
library(edgeR)

# Read in sample metadata
# The file 'samples.txt' should have a Filename column (e.g., "LMR-P1-28-0/quant/quant.sf"),
# a Sample column (e.g., "P1-28"), and a Sugar column ("treated" or "untreated").
samples <- read.table(file = "samples.txt", header = TRUE)

# Build file paths to quantification outputs from Salmon
filenames <- paste(samples$Filename)
files <- file.path(filenames)
names(files) <- paste0(1:nrow(samples))  # Adjusted to match number of samples

# Confirm all files exist
stopifnot(all(file.exists(files)))

# Load tx2gene mapping table
# Created with makeTxDbFromGFF() previously; this step assumes it’s pre-generated and available
load("tx2gene.RData")  # tx2gene must be a data.frame with TXNAME and GENEID columns

# Import quantification data from Salmon
txi_all <- tximport(files = files, type = "salmon", tx2gene = tx2gene)

# Extract count matrix
cts <- txi_all$counts

# Normalize transcript lengths
normMat <- txi_all$length
normMat <- normMat / exp(rowMeans(log(normMat)))
normCts <- cts / normMat

# Calculate effective library sizes and normalize
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Create DGEList and attach sample metadata
y <- DGEList(counts = cts)
y$samples <- samples[, c("Sample", "Sugar", "Temp", "BorP", "Media")]

# Explicitly calculate and assign library sizes
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

# Calculate log-transformed counts per million
norm_counts <- cpm(y, log = TRUE)

# Perform PCA on transposed normalized count matrix
pca <- prcomp(t(norm_counts), scale. = TRUE)

# Save PCA results
pc_scores <- as.data.frame(pca$x)
explained_variance <- as.data.frame(pca$sdev^2 / sum(pca$sdev^2))

write.csv(pc_scores, "PCA_scores.csv", row.names = TRUE)
write.csv(explained_variance, "PCA_explained_variance.csv", row.names = FALSE)
