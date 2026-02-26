BiocManager::install("spacexr")

library(Seurat)
library(spacexr)
library(tidyverse)
# Use Seurat v5 to avoid error
# Load Visium data
CT_2_5 <- Load10X_Spatial(data.dir ="/spatial_transcriptomics/OSUMC/counts_and_images/2-5", 
                              filename = "filtered_feature_bc_matrix.h5")

CT_2_5 <- NormalizeData(CT_2_5)
CT_2_5 <- FindVariableFeatures(CT_2_5)
CT_2_5 <- ScaleData(CT_2_5)

# Sketch the cortical subset of the Visium HD dataset
CT_2_5 <- SketchData(object = CT_2_5, ncells = 50000, method = "LeverageScore", 
                     sketched.assay = "sketch")

DefaultAssay(CT_2_5) <- "sketch"
CT_2_5 <- FindVariableFeatures(CT_2_5)
CT_2_5 <- ScaleData(CT_2_5)
CT_2_5 <- RunPCA(CT_2_5, assay = "sketch", reduction.name = "pca.cortex.sketch", verbose = T)
CT_2_5 <- FindNeighbors(CT_2_5, reduction = "pca.cortex.sketch", dims = 1:50)
CT_2_5 <- FindClusters(CT_2_5, cluster.name = "seurat_cluster.sketched")
CT_2_5 <- RunUMAP(CT_2_5, reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", 
                  return.model = T, dims = 1:50, verbose = T)
DimPlot(AD2_3, label = T)

# Create the RCTD query object using 'SpatialRNA' function
counts_hd <- CT_2_5[["sketch"]]$counts
CT_2_5_cells_hd <- colnames(CT_2_5[["sketch"]])
coords <- GetTissueCoordinates(CT_2_5)[CT_2_5_cells_hd, 1:2]

query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

ref <- readRDS("/spatial_transcriptomics/Allen_middle-temporal-gyrus/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA

# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

# Creating RCTD Object
RCTD_CT_2_5 <- create.RCTD(query, reference, max_cores = 8) # max_cores: for parallel processing.
# run RCTD
RCTD_CT_2_5 <- run.RCTD(RCTD_CT_2_5, doublet_mode = "doublet")

# add results back to Seurat object
CT_2_5 <- AddMetaData(CT_2_5, metadata = RCTD_CT_2_5@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
CT_2_5 <- ProjectData(object = CT_2_5, assay = "Spatial",
                      full.reduction = "pca.cortex",   sketched.assay = "sketch",
                      sketched.reduction = "pca.cortex.sketch", umap.model = "umap.cortex.sketch",
                      dims = 1:50, refdata = list(full_first_type = "first_type"))
DefaultAssay(CT_2_5) <- "Spatial"

DimPlot(T4857, label = T, group.by = "full_first_type")
SpatialDimPlot(CT_2_5, group.by = "full_first_type", label = T, repel = T, label.size = 4)
