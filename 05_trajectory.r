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

# Load necessary packages
library(tidyverse)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(future)
library(rJava)
library(xlsx)
library(stringr)
library(httpgd)
library(monocle3)
library(SeuratWrappers)
library(plotly)
library(slingshot)
library(tradeSeq)
library(condiments)
library(scater)
library(scran)
library(RColorBrewer)
library(BiocParallel)
set.seed(123)

hgd()

# Set seed for reproducibility and specify the number of workers for parallel computing
# set.seed(123)
# nworkers <- 16
# library(future)
# plan("sequential")
# plan("multicore", workers = nworkers)

# Load the processed Seurat object
seu <- readRDS(file = "04_GDMvsNP_allCellType.rds")

# Define genes that might be noise/irrelevant
hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), value = TRUE)
hb_genes <- grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), value = TRUE)
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
         rownames(seu@assays$SCT@counts),
         value = TRUE)
))

colnames(seu@meta.data)

# Visualize UMAP clustering with the selected resolution
umap_plot <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "group",
  group.by = "CellType_l2",
  label = TRUE
) + coord_fixed(ratio = 1)

print(umap_plot)

### Section 1: Global Trajectory Inference --------------------------------------------------------------

# We use slingshot for trajectory inference,

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seu, assay = "RNA")

colnames(colData(sce))
colData(sce)$slingshot_clusters <- colData(sce)$CellType_l2
table(colData(sce)$slingshot_clusters)
sce <- slingshot(sce,
  reducedDim = "UMAP",
  # start.clus = "LD-DEFA3", end.clus = "HD-CXCR4",
  clusterLabels = colData(sce)$slingshot_clusters,
  omega = TRUE,
  approx_points = 1000
)

saveRDS(sce, file = "05_sce_for_slingshot.rds")
# sce <- readRDS("05_sce_for_slingshot.rds")

# Determine branch assignments
curveAssignments <- slingBranchID(sce)
sce$curveAssignments <- curveAssignments
table(curveAssignments)

# Plot the results
# p1: source of cell samples
# p2: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:8)[factor(colData(sce)$group)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:3, bty = "n", 
       legend = levels(factor(colData(sce)$group)))

plot(reducedDims(sce)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce)$slingshot_clusters)))
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = TRUE)

# Plot the results
# p1: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
# p2: Principal curves
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce)$slingshot_clusters)))
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = TRUE)

plot(reducedDims(sce)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce)$slingshot_clusters)))
lines(SlingshotDataSet(sce), lwd = 2, col = 'black')

# Trajectory Inference and Differential Progression
progressionTest(sce, conditions = sce$group,
                lineages = TRUE)

differentiationTest(sce, conditions = sce$group, 
                    global = FALSE, pairwise = TRUE)

# Draw a UMAP plot where each point is a cell and is colored by the average slingshot pseudotime across paths. 
# The principal curves fitted to each lineage are shown in black.
pseudo.paths <- slingPseudotime(sce)
head(pseudo.paths)

shared.pseudo <- rowMeans(pseudo.paths, na.rm = TRUE)
pseudo.paths <- pseudo.paths %>%
  as.data.frame() %>%
  mutate(shared = shared.pseudo) %>%
  mutate(assignment = curveAssignments)
# View(pseudo.paths)
class(curveAssignments)

gg <- plotUMAP(sce, colour_by = I(shared.pseudo))
embedded <- embedCurves(sce, "UMAP")
embedded <- slingCurves(embedded)
# head(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data = embedded, aes(x = umap_1, y = umap_2), size = 1.2)
}
gg

### Section 2: Partial Trajectory Inference --------------------------------------------------------------

# Convert to SingleCellExperiment
sce_sub <- as.SingleCellExperiment(seu_sub, assay = "RNA")

colnames(colData(sce_sub))
colData(sce_sub)$slingshot_clusters <- colData(sce_sub)$CellType_l2
table(colData(sce_sub)$slingshot_clusters)
sce_sub <- slingshot(sce_sub,
  reducedDim = "UMAP",
  # start.clus = "LD-DEFA3", end.clus = "HD-CXCR4",
  clusterLabels = colData(sce_sub)$slingshot_clusters,
  omega = TRUE,
  approx_points = 1000
)

saveRDS(sce_sub, file = "05_sce_sub_for_slingshot.rds")
# sce_sub <- readRDS("05_sce_sub_for_slingshot.rds")

# Determine branch assignments
curveAssignments <- slingBranchID(sce_sub)
sce_sub$curveAssignments <- curveAssignments
table(curveAssignments)

# Plot the results
# p1: source of cell samples
# p2: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce_sub)$UMAP,
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
     col = alpha(c(1:8)[factor(colData(sce_sub)$group)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:3, bty = "n", 
       legend = levels(factor(colData(sce_sub)$group)))

plot(reducedDims(sce_sub)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce_sub)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce_sub)$slingshot_clusters)))
lines(SlingshotDataSet(sce_sub), lwd = 2, type = 'lineages', col = 'black', show.constraints = TRUE)

# Plot the results
# p1: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
# p2: Principal curves
layout(matrix(1:2, nrow = 1))
par(mar = c(4.5,4,1,1))

plot(reducedDims(sce_sub)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce_sub)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce_sub)$slingshot_clusters)))
lines(SlingshotDataSet(sce_sub), lwd = 2, type = 'lineages', col = 'black', show.constraints = TRUE)

plot(reducedDims(sce_sub)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce_sub)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce_sub)$slingshot_clusters)))
lines(SlingshotDataSet(sce_sub), lwd = 2, col = 'black')

# Trajectory Inference and Differential Progression
progressionTest(sce_sub, conditions = sce_sub$group,
                lineages = FALSE)

# Draw a UMAP plot where each point is a cell and is colored by the average slingshot pseudotime across paths. 
# The principal curves fitted to each lineage are shown in black.
pseudo.paths <- slingPseudotime(sce_sub)
head(pseudo.paths)

shared.pseudo <- rowMeans(pseudo.paths, na.rm = TRUE)
pseudo.paths <- pseudo.paths %>%
  as.data.frame() %>%
  mutate(shared = shared.pseudo) %>%
  mutate(assignment = curveAssignments)
# View(pseudo.paths)
class(curveAssignments)

gg <- plotUMAP(sce_sub, colour_by = I(shared.pseudo))
embedded <- embedCurves(sce_sub, "UMAP")
embedded <- slingCurves(embedded)
# head(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data = embedded, aes(x = umap_1, y = umap_2), size = 1.2)
}
gg

# Determine branch assignments again (if needed)
curveAssignments <- slingBranchID(sce_sub)
sce_sub$curveAssignments <- curveAssignments
table(curveAssignments)
sce_sub.uni <- sce_sub[, sce_sub$curveAssignments %in% c("1")]

layout(matrix(1:1, nrow = 1))
plot(reducedDims(sce_sub.uni)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(2:4)[factor(colData(sce_sub.uni)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 2:4, bty = "n", legend = levels(factor(colData(sce_sub.uni)$slingshot_clusters)))
lines(SlingshotDataSet(sce_sub.uni), lwd = 2, col = 'black')