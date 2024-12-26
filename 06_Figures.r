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

hgd()

# Set seed for reproducibility and specify the number of workers for parallel computing
set.seed(123)

# Load the processed Seurat objects
seu <- readRDS(file = "04_GDMvsNP_allCellType.rds")
seu_sub <- readRDS(file = "04_LD-S100A4_vs_LD-MMP9.rds")

# Define genes that might be noise/irrelevant
hist_genes <- grep("Hist", rownames(seu@assays$SCT@counts), value = TRUE)
hb_genes <- grep("^Hb[ab]-|^HB[^(P)]", rownames(seu@assays$SCT@counts), value = TRUE)
bad_features <- unique(c(
    hist_genes, hb_genes,
    grep("^mt-|^Mtmr|^MT-|MTRNR2L|Mtrnr2l|Rp[sl]|^RP[SL]|^HSP|^DNAJ|^Hsp|^Dnaj|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc",
         rownames(seu@assays$SCT@counts),
         value = TRUE)
))

# seu$CellType_l2 <- dplyr::recode(seu$CellType_l2,
#     "HD-NFKBIA" = "HD-NFKBIA",
#     "LD-OLFM4" = "LD-OLFM4",
#     "LD-PADI4" = "LD-OLR1",
#     "LD-DEFA3" = "LD-DEFA3",
#     "LD-MMP9" = "LD-MMP9",
#     "HD-CXCR4" = "HD-CXCR4",
#     "HD-CXCL1" = "HD-CXCL1",
#     "LD-S100A4" = "LD-S100A4"
# )

# seu$CellType_l2 <- factor(seu$CellType_l2,
#     levels = c("LD-DEFA3", "LD-OLFM4", "LD-OLR1", "LD-MMP9", "LD-S100A4", "HD-NFKBIA", "HD-CXCL1", "HD-CXCR4")
# )
levels(seu$CellType_l2)

# Define color panel
col <- c(
    "LD-DEFA3" = "#FA9645",
    "LD-OLFM4" = "#4fc2dc",
    "LD-OLR1" = "#F684C1",
    "LD-MMP9" = "#AA8ED6",
    "LD-S100A4" = "#EB6F5D",
    "HD-NFKBIA" = "#e0c52e",
    "HD-CXCL1" = "#4cb842",
    "HD-CXCR4" = "#8C564B"
)

### Section 1: Dimension Reduction --------------------------------------------------

# Plot 1: UMAP and clusters 
DimPlot(seu,
    reduction = "umap",
    cols = col,
    #split.by = "group",
    group.by = "CellType_l2", 
    label = TRUE
) +
    coord_fixed(ratio = 1)
ggsave(file = "figures/fig1_umap_cluster.png", width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig1_umap_cluster.pdf", width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')

### Section 2: Dimension Reduction - Split --------------------------------------------------

DimPlot(
  seu,
  reduction = "umap",
  cols = col,
  split.by = "Sample_Name",
  group.by = "CellType_l2",
  label = FALSE
) + coord_fixed(ratio = 1)
ggsave(file = "figures/fig2_umap_cluster_split.png", width = 250, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig2_umap_cluster_split.pdf", width = 250, height = 150, units = "mm", device = "pdf", bg = 'transparent')

### Section 3: Cluster Abundance in Clusters-Bar Chart --------------------------------------------------

source_cluster <- tibble(Cluster = seu$CellType_l2,
                        Group = seu$Sample_Name) %>%
  group_by(Cluster, Group) %>%
  summarise(no.cell = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    total.no = sum(no.cell),
    perc = 100 * no.cell / total.no
  ) %>%
  select(Cluster, Group, perc) 
write.csv(source_cluster %>% spread(Cluster, perc), file = "results/06_Sample_CellType_abundance.csv")

ggplot(
  source_cluster,
  aes(x = Group, y = perc, fill = Cluster)
) +
  geom_col(colour = "black") +
  scale_fill_manual(values = col) +
  coord_fixed(ratio = 1 / 10) +
  theme_bw() +
  xlab("") +
  ylab("%")
ggsave(file = "figures/fig3_cluster_abundance.png", width = 500, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig3_cluster_abundance.pdf", width = 500, height = 300, units = "mm", device = "pdf", bg = 'transparent')

table(seu$Sample_Name, seu$CellType_l2)

### Section 4: Cluster Markers --------------------------------------------------

load("04_GDM_all_clusters_markers.rds")
Idents(seu) <- seu$CellType_l2

# Define the genes of interest based on the updated manuscript description
# Mature subsets (Mat-NFKBIA, Mat-CXCR4, Mat-CXCL1) associated genes:
# B2M, BTG2, CSF3R, CXCL8, FCGR3B, JUNB, NAMPT, NFKBIA, ZFP36, IL1A, IL1B, CCL3, CCL4, CXCL1, CXCR4
#
# Immature subsets (Imm-DEFA3, Imm-OLR1, Imm-OLFM4, Imm-MMP9, Imm-S100A4) associated genes:
# DEFA3, OLR1, OLFM4, MMP9, S100A4, plus CAMP, LTF, LCN2, MMP8, CEACAM8
#
# Bridging subset genes:
# MME (CD10), TLR2, FCGR2B (CD32), OLR1 (LOX-1) (already included OLR1 above)
#
# Combine all mentioned genes (remove duplicates if any)
feature_genes <- c(
  # Immature-related
  "DEFA3", "OLFM4", "OLR1", "CAMP", "LTF", "LCN2", "MMP8", "CEACAM8","MMP9", "S100A4", 
  # Mature-related
  "B2M", "BTG2", "CSF3R", "CXCL8", "FCGR3B", "JUNB", "NAMPT",  "ZFP36", "NFKBIA",
  "CXCL1", "CXCR4"
)

DotPlot(seu,
    assay = "SCT",
    # group.colors = cl,
    c("white", "#AD272B"),
    features = rev(feature_genes)
) +
  # coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title = element_blank(),
        legend.direction = "vertical",
        legend.position = "bottom")

ggsave(file = "figures/fig4_dotmap_clusterDEGs.png", width = 400, height = 150, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig4_dotmap_clusterDEGs.pdf", width = 400, height = 150, units = "mm", device = "pdf", bg = 'transparent')

### Section 5: Cluster Markers in UMAP --------------------------------------------------

library(Nebulosa)

p <- plot_density(
  object = seu, 
  features = c("DEFA3", "OLFM4", "OLR1", "MMP9", "S100A4", "NFKBIA", "CXCL1", "CXCR4"),
  reduction = "umap", 
  method = "wkde",
  adjust = 2
)
p + plot_layout(ncol = 4)

ggsave(file = "figures/fig5_cluster_markers_umap.png", width = 450, height = 200, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig5_cluster_markers_umap.pdf", width = 450, height = 200, units = "mm", device = "pdf", bg = 'transparent')

### Section 6: Differential Genes Between HD vs LD --------------------------------------------------

library(EnhancedVolcano)
library(stringr)
PseudoBulk_DEG <- readRDS("04_PseudobulK_DEG-HD_vs_LD.rds")
PseudoBulk_regulons <- readRDS("04_Pseudobulk_regulons-HD_vs_LD.rds")

# Volcano plot
# Prepare data for volcano plot
DEG_res <- PseudoBulk_DEG %>% as.data.frame()
head(DEG_res)
volcano_DEG <- DEG_res %>%
  rownames_to_column(var = "symbol") %>%
  filter(!str_detect(symbol, "^ENSG|^MIR")) %>% 
  arrange(PValue)

head(volcano_DEG)

p_min_10 <- DEG_res %>% dplyr::filter(abs(logFC) > 2) %>% arrange(PValue) %>% head(5) %>% rownames(.)
logFC_max_5 <- DEG_res %>% dplyr::filter(PValue < 0.001) %>% arrange(desc(logFC)) %>% head(5) %>% rownames(.)
logFC_min_5 <- DEG_res %>% dplyr::filter(PValue < 0.001)  %>% arrange(logFC) %>% head(5) %>% rownames(.)

label_features <- c(label_features, p_min_10, logFC_max_5, logFC_min_5) %>% unique()

DEG_Volcano <- EnhancedVolcano(volcano_DEG,
  lab = volcano_DEG$symbol,
  x = "logFC",
  y = "PValue",
  # xlim = c(-3, 3),
   selectLab = label_features,
  pCutoff = 0.001,
  FCcutoff = 2,
  xlab = bquote(~ Log[2] ~ "fold change"),
  pointSize = 1.0,
  labSize = 6.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  colAlpha = 4 / 5,
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  lengthConnectors = 2,
  arrowheads = FALSE,
  colConnectors = "black"
)
DEG_Volcano

# Prepare data for volcano plot
regulon_res <- PseudoBulk_regulons %>% as.data.frame()
head(regulon_res)
volcano_regulon <- regulon_res %>%
  rownames_to_column(var = "symbol") %>%
  arrange(PValue)

head(volcano_regulon)

p_min_10 <- regulon_res %>% dplyr::filter(abs(logFC) > 2) %>% dplyr::filter(PValue < 0.001) %>% arrange(PValue) %>% head(10) %>% rownames(.)
logFC_max_5 <- regulon_res %>% dplyr::filter(abs(logFC) > 2) %>% dplyr::filter(PValue < 0.001) %>% arrange(desc(logFC)) %>% head(10) %>% rownames(.)
logFC_min_5 <- regulon_res %>% dplyr::filter(abs(logFC) > 2) %>% dplyr::filter(PValue < 0.001)  %>% arrange(logFC) %>% head(10) %>% rownames(.)

label_regulons <- c(p_min_10, logFC_max_5, logFC_min_5) %>% unique()

regulon_Volcano <- EnhancedVolcano(volcano_regulon,
  lab = volcano_regulon$symbol,
  x = "logFC",
  y = "PValue",
  # xlim = c(-3, 3),
   selectLab = label_regulons,
  pCutoff = 0.001,
  FCcutoff = 2,
  xlab = bquote(~ Log[2] ~ "fold change"),
  pointSize = 1.0,
  labSize = 6.0,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  colAlpha = 4 / 5,
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  lengthConnectors = 2,
  arrowheads = FALSE,
  colConnectors = "black"
)
regulon_Volcano

# Combine the two volcano plots into a single plot with two columns
combined_plot <- DEG_Volcano + regulon_Volcano +  plot_layout(ncol = 2)

# Display the combined plot
combined_plot

ggsave(combined_plot, file = "figures/fig6_Volcano_DEGs.png", width = 400, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(combined_plot, file = "figures/fig6_Volcano_DEGs.pdf", width = 400, height = 300, units = "mm", device = "pdf", bg = 'transparent')

### Section 7: GSEA: HD vs LD --------------------------------------------------

library(fgsea)
library(msigdbr)
set.seed(123)
m_df = msigdbr(species = "Homo sapiens")
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat), n = 23)

# Prepare gene list
DEG_res <- PseudoBulk_DEG %>%
  as.data.frame()
head(DEG_res)
logFC <- DEG_res$logFC
names(logFC) <- rownames(DEG_res)
gsea_list <- sort(logFC, decreasing = TRUE)
gsea_list <- gsea_list[!is.na(names(gsea_list))]
length(gsea_list)

# Hallmark
pathwaysDF <- msigdbr("Homo sapiens", category = "H")
pathways <- split(pathwaysDF$gene_symbol, pathwaysDF$gs_name)
fgseaRes <- fgsea(pathways, gsea_list, minSize = 10, maxSize = 500, nproc = 1)
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ", "))
Hres <- fgseaRes %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  arrange(NES)
head(Hres)

NES_max_5 <- Hres %>% arrange(desc(NES)) %>% head(5) %>% .$pathway
NES_min_5 <- Hres %>% arrange(NES) %>% head(5) %>% .$pathway

# Define a general function to create a geneset
create_geneset <- function(pathways_list) {
  term <- c()  # Initialize an empty vector to store pathway names
  gene <- c()  # Initialize an empty vector to store all genes

  # Loop through each pathway in the pathways_list
  for (i in pathways_list) {
    pathway_name <- str_replace(i, "HALLMARK_", "")  # Remove the "HALLMARK_" prefix
    genes <- unlist(pathways[i])  # Unlist the genes in the current pathway
    term <- c(term, rep(pathway_name, length(genes)))  # Repeat the pathway name for each gene
    gene <- c(gene, genes)  # Append the genes to the gene vector
  }

  # Create a dataframe with 'term' and 'gene' columns
  geneset <- data.frame(term = term, gene = gene)
  return(geneset)  # Return the created geneset dataframe
}

# Create the geneset for NES_max_5 pathways
geneset_max <- create_geneset(NES_max_5)

# Create the geneset for NES_min_5 pathways
geneset_min <- create_geneset(NES_min_5)

# View the resulting geneset
head(geneset_max)
head(geneset_min)

# Perform GSEA
library(clusterProfiler)
library(fgsea)
library(enrichplot)
# Running GSEA using the prepared gene list and gene set
HD_GSEA <- GSEA(gsea_list, TERM2GENE = geneset_max, verbose = TRUE)
LD_GSEA <- GSEA(gsea_list, TERM2GENE = geneset_min, verbose = TRUE)

# Plotting GSEA results for the gene sets with adjusted y-axis limits
HD_GSEA_plot <- gseaplot2(HD_GSEA, geneSetID = 1:5, subplots = 1:2)
LD_GSEA_plot <- gseaplot2(LD_GSEA, geneSetID = 1:5, subplots = 1:2)

# Use print() to capture the plots
HD_GSEA_plot <- print(HD_GSEA_plot)
LD_GSEA_plot <- print(LD_GSEA_plot)

ggsave(HD_GSEA_plot, file = "figures/fig7_HD_GSEA.png", width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(HD_GSEA_plot, file = "figures/fig7_HD_GSEA.pdf", width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')

ggsave(LD_GSEA_plot, file = "figures/fig8_LD_GSEA.png", width = 300, height = 300, units = "mm", dpi = 300, device = "png")
ggsave(LD_GSEA_plot, file = "figures/fig8_LD_GSEA.pdf", width = 300, height = 300, units = "mm", device = "pdf", bg = 'transparent')

### Section 8: Global Trajectory --------------------------------------------------

library(slingshot)
library(tradeSeq)
library(condiments)
library(scater)
library(scran)
library(RColorBrewer)
library(BiocParallel)
sce <- readRDS("05_sce_for_slingshot.rds")

# Determine branch assignments
curveAssignments <- slingBranchID(sce)
sce$curveAssignments <- curveAssignments
table(curveAssignments)

# Plot the results
# p1: source of cell samples
# p2: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
layout(matrix(1:1, nrow = 1))
par(mar = c(4.5,4,1,1))
plot(reducedDims(sce)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce)$slingshot_clusters)))
lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = 'black', show.constraints = TRUE)

# Define file names
png_file <- "figures/fig9_global_trajectory.png"
pdf_file <- "figures/fig9_global_trajectory.pdf"

# Set dimensions in millimeters and resolution
width_mm <- 300
height_mm <- 300
dpi <- 300  # Dots per inch

# Convert dimensions to inches
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

# Define the plotting function
plot_global_trajectory <- function() {
  layout(matrix(1:1, nrow = 1))
  par(mar = c(4.5, 4, 1, 1))
  plot(
    reducedDims(sce)$UMAP,
    asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2",
    col = scales::alpha(c(1:8)[factor(colData(sce)$slingshot_clusters)], alpha = 0.5)
  )
  legend(
    "topleft", pch = 16, col = 1:8, bty = "n",
    legend = levels(factor(colData(sce)$slingshot_clusters))
  )
  lines(
    SlingshotDataSet(sce), lwd = 2, type = 'lineages',
    col = 'black', show.constraints = TRUE
  )
}

# Save as PNG
png(
  filename = png_file, width = width_in, height = height_in,
  units = "in", res = dpi
)
plot_global_trajectory()
dev.off()

# Save as PDF
pdf(
  file = pdf_file, width = width_in, height = height_in
)
plot_global_trajectory()
dev.off()

### Section 9: Local Trajectory --------------------------------------------------

# Load required libraries
library(Seurat)
library(slingshot)
library(tradeSeq)
library(ggplot2)
library(patchwork)
library(dittoSeq)

sce_sub <- readRDS("05_sce_sub_for_slingshot.rds")

# Determine branch assignments
curveAssignments <- slingBranchID(sce_sub)
sce_sub$curveAssignments <- curveAssignments
table(curveAssignments)

# Plot the results
# p1: source of cell samples
# p2: the lineage structure estimated by the cluster-based minimum spanning tree by using the type argument.
layout(matrix(1:1, nrow = 1))
par(mar = c(4.5,4,1,1))
plot(reducedDims(sce_sub)$UMAP, 
     asp = 1, pch = 16, xlab = "UMAP-1", ylab = "UMAP-2", 
     col = alpha(c(1:8)[factor(colData(sce_sub)$slingshot_clusters)], alpha = 0.5))
legend("topleft", pch = 16, col = 1:8, bty = "n", legend = levels(factor(colData(sce_sub)$slingshot_clusters)))
lines(SlingshotDataSet(sce_sub), lwd = 2, col = 'black')

# Define file names for Figure 10
png_file_fig10 <- "figures/fig10_local_trajectory.png"
pdf_file_fig10 <- "figures/fig10_local_trajectory.pdf"

# Set dimensions in millimeters and convert to inches
width_mm_fig10 <- 300
height_mm_fig10 <- 300
dpi_fig10 <- 300  # Resolution in dots per inch

width_in_fig10 <- width_mm_fig10 / 25.4
height_in_fig10 <- height_mm_fig10 / 25.4

# Define the plotting function for Figure 10
plot_sub_global_trajectory <- function() {
  layout(matrix(1:1, nrow = 1))
  par(mar = c(4.5, 4, 1, 1))
  
  # Plot UMAP with clusters
  plot(
    reducedDims(sce_sub)$UMAP,
    asp = 1, 
    pch = 16, 
    xlab = "UMAP-1", 
    ylab = "UMAP-2", 
    col = scales::alpha(c(1:8)[factor(colData(sce_sub)$slingshot_clusters)], alpha = 0.5)
  )
  
  # Add legend
  legend(
    "topleft", 
    pch = 16, 
    col = 1:8, 
    bty = "n", 
    legend = levels(factor(colData(sce_sub)$slingshot_clusters))
  )
  
  # Add Slingshot lineages
  lines(
    SlingshotDataSet(sce_sub), 
    lwd = 2, 
    col = 'black'
  )
}

# Save Figure 10 as PNG
png(
  filename = png_file_fig10, 
  width = width_in_fig10, 
  height = height_in_fig10,
  units = "in", 
  res = dpi_fig10
)
plot_sub_global_trajectory()
dev.off()

# Save Figure 10 as PDF
pdf(
  file = pdf_file_fig10, 
  width = width_in_fig10, 
  height = height_in_fig10
)
plot_sub_global_trajectory()
dev.off()

### Section 10: Differential Differentiation Between GDM and NP --------------------------------------------------

# Load necessary libraries
library(Seurat)
library(slingshot)
library(tradeSeq)
library(ggplot2)
library(patchwork)
library(dittoSeq)

# Load the Seurat object containing your subset of cells
seu_sub <- readRDS(file = "04_LD-S100A4_vs_LD-MMP9.rds")

# Load the SingleCellExperiment object prepared for slingshot
sce_sub <- readRDS("05_sce_sub_for_slingshot.rds")

colnames(colData(sce_sub))

# Extract pseudotime values from the slingshot object and convert to a data frame
pseudotime <- slingPseudotime(sce_sub) %>% as.data.frame()

# Get the names of the lineages (e.g., "Lineage1", "Lineage2", etc.)
Lineages = colnames(pseudotime)

# Add pseudotime information to the Seurat object's metadata
# Here, only the first lineage is being added (i = 1)
for(i in 1:1){
  # Extract pseudotime values for the current lineage
  pseudotime_sub <- pseudotime[,i]
  
  # Add the pseudotime values as a new metadata column in the Seurat object
  seu_sub <- AddMetaData(
    object = seu_sub,
    metadata = pseudotime_sub,
    col.name = Lineages[i]
  )
}

# Display the column names of the Seurat object's metadata to verify addition
colnames(seu_sub@meta.data)

# Define a custom function to plot the density of pseudotime distributions
pseudotime_density <- function(seurat_obj,
                               Lineage,
                               cluster_label,
                               color,
                               alpha = 0.5){
  # Create a data frame with cell type and pseudotime information
  df <- data.frame(seurat_obj[[cluster_label]], seurat_obj[[Lineage]]) 
  
  # Rename the columns for clarity
  colnames(df) <- c("celltype", "Lineage")
  
  # Remove any rows with missing values
  df <- na.omit(df)
  
  # Generate the density plot using ggplot2
  p <- ggplot(df, aes(x = Lineage, fill = celltype)) +
    geom_density(alpha = alpha) +                       # Use the alpha parameter here
    theme_bw() +                                      # Use a clean theme
    scale_fill_manual(values = color) +               # Manually set fill colors
    labs(
      title = paste("Pseudotime Density for", Lineage),
      x = "Pseudotime",
      y = "Density",
      fill = cluster_label
    )
  
  # Return the plot object
  return(p)
}

# Generate the first density plot:
# - Lineage: "Lineage1"
# - Cluster Label: "CellType_l2"
# - Color: predefined color palette 'col'
p1 <- pseudotime_density(
  seurat_obj = seu_sub,
  Lineage = "Lineage1",
  cluster_label = "CellType_l2",
  color = col
) + 
  theme(axis.title.x = element_blank())  # Remove x-axis title for the first plot

# Generate the second density plot:
# - Lineage: "Lineage1"
# - Cluster Label: "group"
# - Color: specific colors for groups
p2 <- pseudotime_density(
  seurat_obj = seu_sub,
  Lineage = "Lineage1",
  cluster_label = "group",
  color = dittoColors(),  # Replace with 'custom_colors' if preferred
  alpha = 0.2             # Adjusted alpha value
) 
p2

# Combine the two plots vertically using patchwork
combined_plot <- p1 + p2 + 
  plot_layout(
    ncol = 1,          # Arrange plots in one column (stacked vertically)
    heights = c(1, 2), # Equal heights for both plots
    guides = 'collect' # Collect legends into a single combined legend
  )

# Display the combined plot
print(combined_plot)

# Save the combined plot as PNG
ggsave(
  combined_plot,
  filename = "figures/fig11_differential_differentiation.png",
  width = 300,
  height = 300,
  units = "mm",
  dpi = 300,
  device = "png"
)

# Save the combined plot as PDF with transparent background
ggsave(
  combined_plot,
  filename = "figures/fig11_differential_differentiation.pdf",
  width = 300,
  height = 300,
  units = "mm",
  device = "pdf",
  bg = 'transparent'
)

### Section 11: Local Feature Genes --------------------------------------------------

library(Nebulosa)

features <- c("MME", "FCGR2B", "IL13RA1", "CCR1","TREM1","CXCR2", "ITGAX", "CD14","TLR2", "CR1", "CLEC2B", "OLR1", "CD101","TACSTD2","CD24","CD200R1")

p <- plot_density(
  object = seu_sub, 
  features = features,
  reduction = "umap", 
  method = "wkde",
  adjust = 2
)
p + plot_layout(ncol = 4)

ggsave(file = "figures/fig12_local_feature_genes_umap.png", width = 450, height = 200, units = "mm", dpi = 300, device = "png")
ggsave(file = "figures/fig12_local_feature_genes_umap.pdf", width = 450, height = 200, units = "mm", device = "pdf", bg = 'transparent')

### Section 12: Split Trajectory & Feature Gene Expression --------------------------------------------------

# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Ensure that the default assay is properly set
DefaultAssay(seu_sub) <- "SCT"

# Plot 1: UMAP and clusters (split by group)
# Note: Ensure that 'group' and 'CellType_l2' columns exist in seu_sub@meta.data
# 'col' is commented out; if you have a predefined color palette, you can uncomment and define it.
p1 <- DimPlot(
  object = seu_sub,
  reduction = "umap",
  # cols = col,  # Uncomment and define 'col' if you want custom colors
  split.by = "group",
  group.by = "CellType_l2",
  label = TRUE
) + coord_fixed(ratio = 1)

# Genes of interest
genes_of_interest <- c("MME","TLR2","FCGR2B","OLR1")

# Fetch pseudotime and gene expression data from the SCT assay data slot
data_df <- FetchData(seu_sub, vars = c("Lineage1", genes_of_interest))

# Perform min-max normalization for each gene independently
# For each gene: normalized_value = (value - min_gene) / (max_gene - min_gene)
data_norm <- data_df %>%
  mutate(across(all_of(genes_of_interest), 
                ~ ( . - min(.) ) / ( max(.) - min(.) )))

# Convert data to long format for ggplot
long_df <- data_norm %>%
  pivot_longer(
    cols = all_of(genes_of_interest), 
    names_to = "gene", 
    values_to = "expression"
  )

# Define colors for each gene
gene_colors <- c(
  "MME" = "#FFA500",    # corresponding to CD10
  "TLR2" = "#E68600",   # corresponding to CD282
  "FCGR2B" = "#1E90FF", # corresponding to CD32
  "OLR1" = "#87CEEB"    # corresponding to LOX-1
)

# Plot 2: Normalized gene expression over pseudotime with LOESS smoothing
# No points are shown, only smooth lines
p2 <- ggplot(long_df, aes(x = Lineage1, y = expression, color = gene)) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  scale_color_manual(values = gene_colors) +
  theme_bw() +
  labs(
    x = "Pseudotime",
    y = "Normalized Expression",
    title = "Normalized Gene Expression over Pseudotime"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Combine p1 and p2 vertically
final_plot <- p1 / p2 + 
  plot_layout(
    heights = c(1, 2)  # Adjust heights if needed
  )

# Display the combined plot
print(final_plot)

# Save the combined figure as PNG
ggsave(
  filename = "figures/fig13_gene_expression_fitted_lines.png",
  plot = final_plot,
  width = 150,
  height = 100,
  units = "mm",
  dpi = 300
)

# Save the combined figure as PDF with transparent background
ggsave(
  filename = "figures/fig13_gene_expression_fitted_lines.pdf",
  plot = final_plot,
  width = 150,
  height = 100,
  units = "mm",
  device = "pdf",
  bg = 'transparent'
)
