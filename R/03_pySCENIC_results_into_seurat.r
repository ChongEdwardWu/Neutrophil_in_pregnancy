### Section 0: Preparation --------------------------------------------------

# Suppress warnings and messages
suppressMessages(suppressWarnings(source("path_to_.radian_profile")))

# Clear the environment
rm(list = ls())
graphics.off()
gc()

# Set your work directory
workdir <- "path_to_data"
setwd(workdir)

# Create directories for figures and results
if (!file.exists(file.path(workdir, "figures"))) {
    dir.create(file.path(workdir, "figures"))
}
if (!file.exists(file.path(workdir, "results"))) {
    dir.create(file.path(workdir, "results"))
}

### Section 1: Load SCENIC Results --------------------------------------------------------------

# Required packages:
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(data.table)
library(tidyverse)
set.seed(123)

# Load the processed SCENIC results
seu <- readRDS("01_3_seu_with_dimRed_and_clusters.rds")

### Section 2: Filtering Genes for pySCENIC ---------------------------------------------------------------

# Extract the expression matrix from the SCT assay
exprMat <- seu@assays$SCT@data

# Extract cell metadata
cellInfo <- seu@meta.data

# Define a threshold for gene expression
# Keep genes expressed in at least 1% of cells
threshold <- log1p(1 * 0.01 * ncol(exprMat))  # log1p(1% of cells)

# Identify genes expressed above the threshold
loci1 <- which(rowSums(exprMat) > threshold)
table(rowSums(exprMat) > threshold)
exprMat_filter <- exprMat[loci1, ]
dim(exprMat_filter)

### Section 3: Transfer into Loom Files ---------------------------------------------------------------

# Define the SCENIC results directory and file paths
SCENICdir <- "path_to_SCENIC_results"
scenicLoomPath <- file.path(SCENICdir, '02_pySCENIC_seu_filtered.scenic.loom')
motifEnrichmentFile <- file.path(SCENICdir, '02_pySCENIC_seu_filtered.motifs.csv')
regulonAucFile <- file.path(SCENICdir, '02_pySCENIC_seu_filtered.auc.csv')
BinarymatFile <- file.path(SCENICdir, '02_pySCENIC_seu_filtered.bin.csv')
regulonAucThresholdsFile <- file.path(SCENICdir, '02_pySCENIC_seu_filtered.thresholds.csv')

# Check if all required files exist
all(file.exists(scenicLoomPath), file.exists(motifEnrichmentFile), file.exists(regulonAucFile), file.exists(BinarymatFile), file.exists(regulonAucThresholdsFile))

# Retrieve AUC scores per cell from the SCENIC loom file
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
close_loom(loom)

# Read motif enrichment and regulon AUC data
motifEnrichment <- fread(motifEnrichmentFile, header = TRUE, skip = 1)[-1,]
colnames(motifEnrichment)[1:2] <- c("TF", "MotifID")

regulonAUC2 <- fread(regulonAucFile, header = TRUE) %>%
  column_to_rownames("Cell") %>%
  t()
rownames(regulonAUC2) <- gsub("[(+)]", "", rownames(regulonAUC2))

regulonBin <- fread(BinarymatFile, header = TRUE) %>%
  column_to_rownames("Cell") %>%
  t()
rownames(regulonBin) <- gsub("[(+)]", "", rownames(regulonBin))

# Get AUC matrix
AUCmat <- AUCell::getAUC(regulonAUC)
rownames(AUCmat) <- gsub("[(+)]", "", rownames(AUCmat))

# Save the processed SCENIC results
save(regulons_incidMat, regulons, regulonAUC, motifEnrichment, regulonAUC2, regulonBin, AUCmat,
  file = "03_pySCENIC_results.RData"
)

### Section 4: Incorporate Data into Seurat Object --------------------------------------------------------------

# Load required packages
library(Seurat)

# Load the Seurat object
seu <- readRDS("01_3_seu_with_dimRed_and_clusters.rds")

# Merge AUC matrix into the Seurat object
seu[['AUC']] <- CreateAssayObject(data = regulonAUC2)
seu <- ScaleData(seu, assay = 'AUC', features = rownames(regulonAUC2))

# Merge AUC-bin matrix into the Seurat object
seu[['Bin']] <- CreateAssayObject(data = regulonBin)
seu <- ScaleData(seu, assay = 'Bin', features = rownames(regulonBin))

# Save the updated Seurat object
DefaultAssay(seu) <- 'SCT'
saveRDS(seu, file = "03_pySCENIC2seurat.rds")

# Print success message
cat("pySCENIC results have been incorporated into the Seurat object, and Seurat object saved successfully.\n")

# Optional: Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")
