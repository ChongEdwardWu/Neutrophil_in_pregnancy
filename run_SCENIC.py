#!/usr/bin/env python
import os
import pandas as pd
import loompy as lp
from pyscenic.binarization import binarize

# Variables that need to be set manually
DATASET_ID = '02_pySCENIC_seu_filtered'
WORKDIR = "path_to_data"
AUXILLIARIES_FOLDERNAME = "path_to_auxiliaries"  # Change species if necessary!

# Set variables for file paths to read from and write to:
R_FOLDERNAME = os.path.join(WORKDIR, "03_R")
RESULTS_FOLDERNAME = os.path.join(WORKDIR, "04_SCENIC/results")
FIGURES_FOLDERNAME = os.path.join(WORKDIR, "04_SCENIC/figures")

# Path to loom file converted from the Seurat Object.
INPUTLOOM_FNAME = os.path.join(R_FOLDERNAME, f"{DATASET_ID}.loom")

# Ranking databases.
RANKING_DBS_FNAMES = " ".join(
    [os.path.join(AUXILLIARIES_FOLDERNAME, f) for f in os.listdir(AUXILLIARIES_FOLDERNAME) if f.endswith(".rankings.feather")]
)

# Motif annotations.
MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, "allTFs_hg38.txt")

# Output files
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.adjacencies.tsv")
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.motifs.csv")
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.regulons.dat")
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.auc.csv")
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.bin.csv")
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.thresholds.csv")
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.scenic.loom")
UMAP_FNAME = os.path.join(RESULTS_FOLDERNAME, f"{DATASET_ID}.umap.csv")

# STEP 1: Gene regulatory network inference
if os.system(
    f"pyscenic grn {INPUTLOOM_FNAME} {MM_TFS_FNAME} -o {ADJACENCIES_FNAME} --num_workers 4"
) == 0:
    print("SCENIC STEP 1 is done!")
else:
    print("Error in STEP 1: GRN inference")

# STEP 2: Regulon prediction (cisTarget)
if os.system(
    f"pyscenic ctx {ADJACENCIES_FNAME} {RANKING_DBS_FNAMES} --annotations_fname {MOTIF_ANNOTATIONS_FNAME} "
    f"--expression_mtx_fname {INPUTLOOM_FNAME} --output {MOTIFS_FNAME} --auc_threshold 0.05 --num_workers 4"
) == 0:
    print("SCENIC STEP 2 is done!")
else:
    print("Error in STEP 2: Regulon prediction")

# STEP 3: AUC matrix calculation and LOOM file creation
if os.system(
    f"pyscenic aucell {INPUTLOOM_FNAME} {MOTIFS_FNAME} --output {LOOM_FNAME} --num_workers 4"
) == 0:
    print("SCENIC STEP 3 is done!")
else:
    print("Error in STEP 3: AUC matrix calculation and LOOM creation")

# STEP 4: Regulon activity binarization
auc_mtx = pd.read_csv(AUCELL_MTX_FNAME, index_col=0)
bin_mtx, thresholds = binarize(auc_mtx, seed=123, num_workers=4)
bin_mtx.to_csv(BIN_MTX_FNAME)
thresholds.to_frame().rename(columns={0: 'threshold'}).to_csv(THR_FNAME)
print("SCENIC STEP 4 is done!")
