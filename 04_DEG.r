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
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(future)
library(rJava)
library(xlsx)
library(stringr)
library(httpgd)
library(clusterProfiler)
# library(org.Mm.eg.db) #! For mouse
library(org.Hs.eg.db) #! For Human
library(fgsea)
library(msigdbr)
library(scater)
library(DESeq2)
library(SingleCellExperiment)
library(edgeR)

hgd()

# Set seed for reproducibility and specify the number of workers for parallel computing
set.seed(123)
nworkers <- 4
library(future)
plan("sequential")
# plan("multicore", workers = nworkers)

# Load the processed Seurat object
seu <- readRDS("03_pySCENIC2seurat.rds")

# Define genes that might be noise/irrelevant
hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), value = TRUE)
hb_genes <- grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), value = TRUE)
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
         rownames(seu@assays$SCT@counts),
         value = TRUE)
))

features <- setdiff(rownames(seu@assays$SCT@data), bad_features)

# Reorganize group levels
seu$Sample_Name <- factor(seu$Sample_Name,
    levels = c("NP_HD", "NP_LD", "GDM_HD", "GDM_LD")
)
levels(seu$Sample_Name)

### Section 1: Clustering --------------------------------------------------------------

# Choose the desired resolution based on clustree and other criteria
res <- 0.2
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

### Section 2: Identify Cluster Markers --------------------------------------------------------------

# Load necessary packages for this section
library(future)
plan("multicore", workers = 16)

# Set the default assay and prepare the object for finding markers
colnames(seu@meta.data)

# Idents(seu) <- "seurat_clusters"
Idents(seu) <- "CellType_l2"
DefaultAssay(seu) <- "SCT"
seu <- PrepSCTFindMarkers(object = seu)

# Wilcoxon rank sum test
seu_DEG <- FindAllMarkers(seu,
  only.pos = TRUE,
  min.pct = 0.1,
  features = features,
  logfc.threshold = 0.25,
  densify = TRUE
)
# Arrange the results
library(dplyr)
seu_DEG <- seu_DEG %>%
  group_by(cluster) %>%
  arrange(p_val_adj, .by_group = TRUE)
# View(seu_DEG)

# ROC test
seu_ROC <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = features,
    test.use = "roc",
    densify = TRUE
)

# Regulon assay
seu_regulon <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = NULL,
    test.use = "roc",
    assay = "AUC",
    logfc.threshold = 0.005,
    densify = TRUE
)

# Regulon Bin assay
seu_regulon_bin <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    features = NULL,
    test.use = "roc",
    assay = "Bin",
    logfc.threshold = 0.005,
    densify = TRUE
)

# Save results
save(seu_DEG, seu_ROC, seu_regulon, seu_regulon_bin,
  file = "04_GDM_all_clusters_markers.rds"
)

# For exporting results
library(rJava)
library(xlsx)
library(stringr)
jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
  gc()
}

write.xlsx(seu_DEG %>% as.data.frame(), 
           file = "results/04_clusters_markers.xlsx", 
           row.names = FALSE,
           sheetName = "clusters_DEG", 
           append = FALSE)
jgc()
write.xlsx(seu_ROC %>% as.data.frame(), 
           file = "results/04_clusters_markers.xlsx", 
           row.names = FALSE,
           sheetName = "clusters_ROC", 
           append = TRUE)
jgc()
write.xlsx(seu_regulon %>% as.data.frame(), 
           file = "results/04_clusters_markers.xlsx", 
           row.names = FALSE,
           sheetName = "clusters_regulon", 
           append = TRUE)
jgc()
write.xlsx(seu_regulon_bin %>% as.data.frame(), 
           file = "results/04_clusters_markers.xlsx", 
           row.names = FALSE,
           sheetName = "clusters_regulon_bin", 
           append = TRUE)
jgc()

# Print success message
cat("All cluster markers determined.\n")

# Based on these markers, define clusters
seu$CellType_l2 <- factor(seu$seurat_clusters)
table(seu$CellType_l2)

# Modify the CellType_l2 values based on the new labels
seu$CellType_l2 <- dplyr::recode(seu$seurat_clusters,
    "1" = "HD-NFKBIA",
    "2" = "LD-OLFM4",
    "3" = "LD-OLR1",
    "4" = "LD-DEFA3",
    "5" = "LD-MMP9",
    "6" = "HD-CXCR4",
    "7" = "HD-CXCL1",
    "8" = "LD-S100A4"
)
seu$CellType_l2 <- factor(seu$CellType_l2,
    levels = c("LD-DEFA3", "LD-OLFM4", "LD-OLR1", "LD-MMP9", "LD-S100A4", "HD-NFKBIA", "HD-CXCL1", "HD-CXCR4")
)
levels(seu$CellType_l2)

# Visualize UMAP clustering with the selected resolution
Idents(seu) <- "CellType_l2"
umap_plot <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "Sample_Name",
  group.by = "CellType_l2",
  label = FALSE
) + coord_fixed(ratio = 1)

print(umap_plot)

# Save the UMAP plot
ggsave(
  filename = "figures/03_UMAP_clusters.png",
  plot = umap_plot,
  width = 12,
  height = 4
)

# Summarize cluster distribution across groups
source_cluster <- tibble(
  Cluster = seu$CellType_l2,
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
  filename = "figures/03_cluster_distribution_plot.png",
  plot = cluster_distribution_plot,
  width = 12,
  height = 8
)

### Section 3: Identify Differential Genes Between LD-MMP9 vs LD-S100A4 --------------------------------------------------

# Visualize UMAP clustering with the selected resolution
umap_plot <- DimPlot(
  seu,
  reduction = "umap",
  split.by = "group",
  group.by = "CellType_l2",
  label = FALSE
) + coord_fixed(ratio = 1)

print(umap_plot)

jgc <- function(){
  .jcall("java/lang/System", method = "gc")
  gc()
}

seu_sub <- seu[, seu$CellType_l2 %in% c("LD-MMP9", "LD-S100A4")]

# Extract UMAP embeddings
umap_embeddings <- Embeddings(seu_sub, "umap")
cells_to_keep <- rownames(umap_embeddings)[umap_embeddings[, 2] > 0]
# Subset the Seurat object to include only these cells
seu_sub <- subset(seu_sub, cells = cells_to_keep)
# Check the new Seurat object
seu_sub

features <- setdiff(rownames(seu_sub@assays$SCT@data), bad_features)
Idents(seu_sub) <- "CellType_l2"

## Differential Expression Analysis
resultfile <- "results/04_LD-S100A4_vs_LD-MMP9-DEG.xlsx"  # Set the result file name here

DEG_SCT <- FindMarkers(
        seu_sub,
        assay = "SCT",
        test.use = "wilcox",
        ident.1 = "LD-S100A4",
        ident.2 = "LD-MMP9",
        min.pct = 0.01,
        only.pos = FALSE,
        logfc.threshold = 0,
        recorrect_umi = FALSE,
        densify = FALSE
)
write.xlsx(DEG_SCT,
        file = resultfile,
        row.names = TRUE,
        sheetName = "DEG",
        append = FALSE
)
DEG_AUC <- FindMarkers(
        seu_sub,
        assay = "AUC",
        test.use = "wilcox",
        ident.1 = "LD-S100A4",
        ident.2 = "LD-MMP9",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0.000,
        densify = TRUE
)
write.xlsx(DEG_AUC,
        file = resultfile,
        row.names = TRUE,
        sheetName = "Regulon AUC scores",
        append = TRUE
)
DEG_Bin <- FindMarkers(
        seu_sub,
        assay = "Bin",
        test.use = "wilcox",
        ident.1 = "LD-S100A4",
        ident.2 = "LD-MMP9",
        min.pct = 0.1,
        only.pos = FALSE,
        logfc.threshold = 0.005,
        densify = TRUE
)
write.xlsx(DEG_Bin,
        file = resultfile,
        row.names = TRUE,
        sheetName = "Binary regulon AUC scores",
        append = TRUE
)

# Enriched GO terms
universe_genes <- na.omit(unique(rownames(seu_sub@assays$SCT)))
markersUP <- DEG_SCT %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC >= 0.25) %>%
    arrange(dplyr::desc(avg_log2FC)) %>%
    rownames() %>%
    na.omit() %>%
    unique()
goUP <- enrichGO(
    gene = markersUP,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db, # For human
    # OrgDb = org.Mm.eg.db, # For mouse
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 50,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    readable = TRUE
)
write.xlsx(goUP,
        file = resultfile,
        row.names = TRUE,
        sheetName = "ChAT_LD-S100A4_enriched_GO",
        append = TRUE
)

markersDN <- DEG_SCT %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC <= -0.25) %>%
    arrange(avg_log2FC) %>%
    rownames() %>%
    na.omit() %>%
    unique()
goDN <- enrichGO(
    gene = markersDN,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db, # For human
    # OrgDb = org.Mm.eg.db, # For mouse
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 50,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    readable = TRUE
)
write.xlsx(goDN,
        file = resultfile,
        row.names = TRUE,
        sheetName = "ChAT_LD-MMP9_enriched_GO",
        append = TRUE
)
jgc()

### Section 4: Gene Set Enrichment Analysis for Clusters: LD-S100A4 vs LD-MMP9 --------------------------------------------------

set.seed(123)
m_df = msigdbr(species = "Homo sapiens")
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n = 23)

# Prepare gene list
logFC <- DEG_SCT$avg_log2FC
names(logFC) <- rownames(DEG_SCT)
gsea_list <- sort(logFC, decreasing = TRUE)
gsea_list <- gsea_list[!is.na(names(gsea_list))]
length(gsea_list)

# Hallmark
pathwaysDF <- msigdbr("Homo sapiens", category = "H")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
Hres <- fgseaRes
head(Hres)
write.xlsx(Hres,
        file = resultfile,
        row.names = TRUE,
        sheetName = "GSEA_HallMark",
        append = TRUE
)
jgc()

# Save the Seurat objects
saveRDS(seu, file = "04_GDMvsNP_allCellType.rds")
# seu <- readRDS(file = "04_GDMvsNP_allCellType.rds")

saveRDS(seu_sub, file = "04_LD-S100A4_vs_LD-MMP9.rds")
# seu_sub <- readRDS(file = "04_LD-S100A4_vs_LD-MMP9.rds")
