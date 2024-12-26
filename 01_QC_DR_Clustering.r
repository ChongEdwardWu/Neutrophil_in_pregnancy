### Section 0: Preparation --------------------------------------------------
source("path_to_.radian_profile")

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

### Section 1: QC --------------------------------------------------

# Load necessary libraries
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(cowplot)
library(httpgd)
library(future)
library(glmGamPoi)
library(clustree)
set.seed(123)
nworkers <- 16

hgd()

# Load the processed Seurat object
seu <- readRDS("path_to_seurat_object/XJY1SMK_Seurat.rds")

# Modify the Sample_Name values based on the new labels
seu$Sample_Name <- dplyr::recode(seu$Sample_Name,
    "SampleTag01_hs" = "GDM_LD",
    "SampleTag02_hs" = "GDM_HD",
    "SampleTag03_hs" = "NP_LD",
    "SampleTag04_hs" = "NP_HD"
)

# Convert Sample_Name to a factor
seu$Sample_Name <- factor(seu$Sample_Name)
table(seu$Sample_Name)

# Step 1: Calculating Quality Metrics (Mitochondrial, Ribosomal, and Custom Genes)
# Calculate percentage of mitochondrial genes
seu <- PercentageFeatureSet(seu, pattern = "^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l", col.name = "percent_mito")

# Calculate percentage of ribosomal genes
seu <- PercentageFeatureSet(seu, pattern = "^Rp[sl]|^RP[SL]", col.name = "percent_ribo")

# Calculate percentage of hemoglobin, heat shock proteins, and platelet genes
seu <- PercentageFeatureSet(seu, pattern = "^Hb[ab]-|^HB[^(P)]", col.name = "percent_hemoglobin")
seu <- PercentageFeatureSet(seu, pattern = "^HSP|^DNAJ|^Hsp|^Dnaj", col.name = "percent_hsp")
seu <- PercentageFeatureSet(seu, pattern = "Pecam1|Pf4|PECAM1|PF4", col.name = "percent_plat")

# Visualize updated metadata
head(seu@meta.data)

# Step 2: Visualizing QC Metrics
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hemoglobin", "percent_hsp", "percent_plat")
VlnPlot(seu, features = qc_features, group.by = "Sample_Name", pt.size = 0.1, ncol = 4) + NoLegend()

# Step 3: Scatter Plots of QC Metrics

# Scatter Plot 1: nCount_RNA vs nFeature_RNA
p1 <- FeatureScatter(seu,
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
    group.by = "Sample_Name", pt.size = 1
) 
p1 + scale_y_continuous(breaks = c(10, 15, 20)) + # Set y-axis ticks
    geom_hline(yintercept = c(200, 300, 500), linetype = "dashed", color = "lightgrey") 
# choose 300

# Scatter Plot 2: nFeature_RNA vs percent_mito (add ticks and dashed lines)
p2 <- FeatureScatter(seu,
    feature1 = "nFeature_RNA", feature2 = "percent_mito",
    group.by = "Sample_Name", pt.size = 1
)
p2 + scale_y_continuous(breaks = c(10, 15, 20)) + # Set y-axis ticks
    geom_hline(yintercept = c(10, 15, 20), linetype = "dashed", color = "lightgrey") # Add dashed lines
# choose 15

# Scatter Plot 3: nFeature_RNA vs percent_ribo
p3 <- FeatureScatter(seu,
    feature1 = "nFeature_RNA", feature2 = "percent_ribo",
    group.by = "Sample_Name", pt.size = 1
)
p3 + geom_hline(yintercept = c(10), linetype = "dashed", color = "lightgrey")
# choose 10

# Scatter Plot 4: percent_hemoglobin vs percent_mito (add ticks and dashed lines)
p4 <- FeatureScatter(seu,
    feature1 = "nFeature_RNA", feature2 = "percent_hemoglobin",
    group.by = "Sample_Name", pt.size = 1
) 
p4 + geom_hline(yintercept = c(1), linetype = "dashed", color = "lightgrey")

# Combine all scatter plots into a grid layout
cowplot::plot_grid(
    p1, p2, p3, p4,
    ncol = 2,
    labels = c("A", "B", "C", "D") # Optional: Add labels to the plots
) + theme(plot.margin = margin(1, 1, 1, 1, "cm")) # Optional: Adjust plot margins

# Step 4: Filtering Low-Quality Cells
# Filter out cells multiplets and undetermined cells
seu <- seu[, !(seu$sample_type_rename %in% c("Multiplet", "Undetermined"))]
seu$sample_type_rename <- factor(seu$sample_type_rename)

# Filter out cells with fewer than 200 detected genes
table(seu$nFeature_RNA > 300)
seu <- subset(seu, subset = nFeature_RNA > 300)

# Filter out cells with high mitochondrial content
table(seu$percent_mito < 15)
seu <- subset(seu, subset = percent_mito < 15)

# Filter out cells with high hemoglobin gene expression (>1%)
table(seu$percent_hemoglobin < 1)
seu <- subset(seu, subset = percent_hemoglobin < 1)

# Save the filtered object after QC
saveRDS(seu, file = "01_1_seu_filtered_after_QC.rds")
# seu <- readRDS("01_1_seu_filtered_after_QC.rds")

# Print success message
cat("QC completed and Seurat object saved successfully.\n")


### Section 2: SCTransform --------------------------------------------------
# Normalization using SCTransform with regression for mitochondrial and ribosomal content

# Load necessary libraries
library(Seurat)
library(sctransform)
library(glmGamPoi)

# Define variables to regress out (mitochondrial and ribosomal content)
vars_to_regress <- c("percent_mito", "percent_ribo")

# Apply SCTransform for normalization using glmGamPoi method
tryCatch(
  {
    seu <- SCTransform(
      seu,
      vars.to.regress = vars_to_regress,
      verbose = TRUE,
      method = "glmGamPoi",
      return.only.var.genes = FALSE  # Deprecated in Seurat v4 and later
    )
    
    # Save the Seurat object after SCTransform with a consistent naming convention
    saveRDS(seu, file = "01_2_seu_after_SCTransform.rds")
    
    # Print success message
    cat("SCTransform completed and Seurat object saved successfully.\n")
  },
  error = function(e) {
    cat("An error occurred during SCTransform:\n")
    print(e)
  }
)

# Optional: Trigger garbage collection to free up memory
gc()

# seu <- readRDS("01_2_seu_after_SCTransform.rds")


### Section 3: Dimension Reduction and Clustering --------------------------------------------------
# Perform dimensionality reduction (PCA, UMAP, t-SNE) and clustering on the Seurat object

# Set the default assay to SCT
DefaultAssay(seu) <- "SCT"

# Set up parallelization using the future package
plan("multicore", workers = nworkers, future.seed = TRUE) # Uncomment if not set earlier
set.seed(123)  # For reproducibility

# Identify and remove undesirable genes before PCA
hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), value = TRUE)
hb_genes <- grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), value = TRUE)
bad_features <- unique(c(
  hist_genes, hb_genes,
  grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
       rownames(seu@assays$SCT@counts),
       value = TRUE)
))

# Select the features for PCA, excluding bad genes
PCA_features <- setdiff(VariableFeatures(seu), bad_features)

# Run PCA with the selected features
seu <- RunPCA(seu, features = PCA_features, npcs = 50, verbose = TRUE)

# Visualize elbow plot to determine the number of PCs to use
ElbowPlot(seu, ndims = 50)

# Run t-SNE and UMAP using the first 20 PCs
# seu <- RunTSNE(seu, reduction = "pca", dims = 1:20)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)

# Find neighbors based on the PCA reduction and the first 20 PCs
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)

# Perform clustering across multiple resolutions to identify the optimal granularity
seu <- FindClusters(seu, resolution = c(seq(0.01, 0.1, 0.01), seq(0.2, 1, 0.1)))

# Relabel clusters for each resolution to ensure they are sequential and start from 1
cluster_res_cols <- grep("SCT_snn_res.", colnames(seu@meta.data), value = TRUE)
for (j in cluster_res_cols) {
  k <- seu@meta.data[[j]]
  levels(k) <- as.character(seq_along(levels(k)))
  seu@meta.data[[j]] <- k
}

# Use clustree to visualize clustering across different resolutions
clustree_plot <- clustree(seu@meta.data, prefix = "SCT_snn_res.", return = "plot")
print(clustree_plot)

# Disable parallelization
plan("sequential")  # Uncomment if not set earlier

# Save the clustree plot for reference
ggsave(
  filename = "figures/01_clustree_plot.png",
  plot = clustree_plot,
  width = 12,
  height = 8
)

# Choose the desired resolution based on clustree and other criteria
res <- 0.08
seu$seurat_clusters <- factor(seu@meta.data[[paste0("SCT_snn_res.", res)]])

# Visualize UMAP clustering with the selected resolution
umap_plot <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "Sample_Name",
  group.by = paste0("SCT_snn_res.", res),
  label = TRUE
) + coord_fixed(ratio = 1)

print(umap_plot)

# Filter out bad clusters if any
# Example: Filter out cluster 4 as outliers or bad cells
discardCl <- seu$seurat_clusters %in% c(4)
table(discardCl)
seu_filt <- seu[, !(discardCl)]
seu <- seu_filt

# Set the chosen resolution for clusters
seu$seurat_clusters <- factor(seu@meta.data[[paste0("SCT_snn_res.", res)]])
Idents(seu) <- seu$seurat_clusters

# Visualize UMAP clustering with the selected resolution
umap_plot <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "Sample_Name",
  group.by = "seurat_clusters",
  label = TRUE
) + coord_fixed(ratio = 1)

print(umap_plot)

# Save the UMAP plot
ggsave(
  filename = "figures/01_UMAP_clusters.png",
  plot = umap_plot,
  width = 12,
  height = 4
)

# Visualize some QC metrics across clusters using violin plots
qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hsp")
vln_plot <- VlnPlot(
  seu,
  features = qc_metrics,
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 3
)

print(vln_plot)

# Save the violin plots
ggsave(
  filename = "figures/01_violin_plots_QC_metrics_clusters.png",
  plot = vln_plot,
  width = 18,
  height = 12
)

# Summarize cluster distribution across groups
source_cluster <- tibble(
  Cluster = seu$seurat_clusters,
  Group = seu$Sample_Name
) %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  select(Cluster, Group, perc)

# Plot cluster distribution using ggplot2
cluster_distribution_plot <- ggplot(source_cluster, aes(x = Group, y = perc, fill = Cluster)) +
  geom_col(colour = "black") +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("Group") +
  ylab("% of cells")

print(cluster_distribution_plot)

# Save the cluster distribution plot
ggsave(
  filename = "figures/01_cluster_distribution_plot.png",
  plot = cluster_distribution_plot,
  width = 12,
  height = 8
)

# Save the Seurat object after dimension reduction and clustering with a consistent naming convention
saveRDS(seu, file = "01_3_seu_with_dimRed_and_clusters.rds")
# seu <- readRDS("01_3_seu_with_dimRed_and_clusters.rds")

# Print success message
cat("Dimension reduction and clustering completed, and Seurat object saved successfully.\n")

# Optional: Save session information for reproducibility
writeLines(capture.output(sessionInfo()), "session_info.txt")