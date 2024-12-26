### Section 0: Preparation --------------------------------------------------

suppressMessages(suppressWarnings(source("path_to_.radian_profile")))

# Clear the environment
rm(list = ls())
graphics.off()
gc()

# Set your working directory
workdir <- "path_to_data"
setwd(workdir)

# Create directories for figures and results
if (!file.exists(file.path(workdir, "figures"))) {
    dir.create(file.path(workdir, "figures"))
}
if (!file.exists(file.path(workdir, "results"))) {
    dir.create(file.path(workdir, "results"))
}

### Section 1: Load Seurat Object ---------------------------------------------------------------

# Load libraries
library(Seurat)
library(tidyverse)
library(magrittr)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)
set.seed(123)

# Load the processed Seurat object after clustering
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

# Define the output loom file path
loom_file <- "02_pySCENIC_seu_filtered.loom"

# Build loom file using SCopeLoomR
loom <- build_loom(loom_file, dgem = exprMat_filter)

# Add cell annotations to the loom file
cell_annotations <- cellInfo %>%
  select(seurat_clusters, Sample_Name)

loom <- add_cell_annotation(loom, cell_annotations)

# Close the loom file to save changes
close_loom(loom)
