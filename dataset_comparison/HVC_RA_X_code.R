library(tidyverse)
library(Seurat)
library(cowplot)
library(egg)
library(qs)
library(Matrix)

#FIRST LOAD SEURAT OBJECTS BEFORE RUNNING CODE BELOW

#HVC_ZF

#CELLS
HVC_ZF <-
  WhichCells(
    object = ZFBF_seurat, 
    expression = (
      species == "zf" &
        position2 == "hvc"
    )
  )

#COUNTS
HVC_ZF_counts <-
  as(
    GetAssayData(ZFBF_seurat, slot = "counts", assay = "RNA")[, HVC_ZF], 
    "sparseMatrix"
  )

#SEURAT OBJECT
HVC_ZF_subset <-
  subset(
    ZFBF_seurat, 
    cells = HVC_ZF
  )

HVC_ZF_seurat <-
  CreateSeuratObject(
    counts = HVC_ZF_counts,
    project = "HVC_ZF"
  )

HVC_ZF_seurat$dataset <- 
  "HVC_ZF"

HVC_ZF_seurat$region <-
  "HVC"

HVC_ZF_seurat$species <-
  "ZF"

HVC_ZF_seurat$cluster_orig <-
  FetchData(HVC_ZF_subset, vars = "cluster_int_sub2")

#ADD MITOCHONDRIAL DATA
HVC_ZF_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

HVC_ZF_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    HVC_ZF_seurat,
    features = HVC_ZF_seurat_mito_genes
  )

#NORMALIZE DATA
HVC_ZF_seurat <-
  NormalizeData(
    HVC_ZF_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
HVC_ZF_seurat <-
  FindVariableFeatures(HVC_ZF_seurat)

#HVC_BF

#CELLS
HVC_BF <-
  WhichCells(
    object = ZFBF_seurat, 
    expression = (
      species == "bf" &
        position2 == "hvc"
    )
  )

#COUNTS
HVC_BF_counts <-
  as(
    GetAssayData(ZFBF_seurat, slot = "counts", assay = "RNA")[, HVC_BF], 
    "sparseMatrix"
  )

#SEURAT OBJECT
HVC_BF_subset <-
  subset(
    ZFBF_seurat, 
    cells = HVC_BF
  )

HVC_BF_seurat <-
  CreateSeuratObject(
    counts = HVC_BF_counts,
    project = "HVC_BF"
  )

HVC_BF_seurat$dataset <- 
  "HVC_BF"

HVC_BF_seurat$region <-
  "HVC"

HVC_BF_seurat$species <-
  "BF"

HVC_BF_seurat$cluster_orig <-
  FetchData(HVC_BF_subset, vars = "cluster_int_sub2")

#ADD MITOCHONDRIAL DATA
HVC_BF_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

HVC_BF_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    HVC_BF_seurat,
    features = HVC_BF_seurat_mito_genes
  )

#NORMALIZE DATA
HVC_BF_seurat <-
  NormalizeData(
    HVC_BF_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
HVC_BF_seurat <-
  FindVariableFeatures(HVC_BF_seurat)

#RA_ZF

#CELLS
RA_ZF <-
  WhichCells(
    object = ZFBF_seurat, 
    expression = (
      species == "zf" &
        position2 == "ra"
    )
  )

#COUNTS
RA_ZF_counts <-
  as(
    GetAssayData(ZFBF_seurat, slot = "counts", assay = "RNA")[, RA_ZF], 
    "sparseMatrix"
  )

#SEURAT OBJECT
RA_ZF_subset <-
  subset(
    ZFBF_seurat, 
    cells = RA_ZF
  )

RA_ZF_seurat <-
  CreateSeuratObject(
    counts = RA_ZF_counts,
    project = "RA_ZF"
  )

RA_ZF_seurat$dataset <- 
  "RA_ZF"

RA_ZF_seurat$region <-
  "RA"

RA_ZF_seurat$species <-
  "ZF"

RA_ZF_seurat$cluster_orig <-
  FetchData(RA_ZF_subset, vars = "cluster_int_sub2")

#ADD MITOCHONDRIAL DATA
RA_ZF_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

RA_ZF_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    RA_ZF_seurat,
    features = RA_ZF_seurat_mito_genes
  )

#NORMALIZE DATA
RA_ZF_seurat <-
  NormalizeData(
    RA_ZF_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
RA_ZF_seurat <-
  FindVariableFeatures(RA_ZF_seurat)

#RA_BF

#CELLS
RA_BF <-
  WhichCells(
    object = ZFBF_seurat, 
    expression = (
      species == "bf" &
        position2 == "ra"
    )
  )

#COUNTS
RA_BF_counts <-
  as(
    GetAssayData(ZFBF_seurat, slot = "counts", assay = "RNA")[, RA_BF], 
    "sparseMatrix"
  )

#SEURAT OBJECT
RA_BF_subset <-
  subset(
    ZFBF_seurat, 
    cells = RA_BF
  )

RA_BF_seurat <-
  CreateSeuratObject(
    counts = RA_BF_counts,
    project = "RA_BF"
  )

RA_BF_seurat$dataset <- 
  "RA_BF"

RA_BF_seurat$region <-
  "RA"

RA_BF_seurat$species <-
  "BF"

RA_BF_seurat$cluster_orig <-
  FetchData(RA_BF_subset, vars = "cluster_int_sub2")

#ADD MITOCHONDRIAL DATA
RA_BF_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

RA_BF_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    RA_BF_seurat,
    features = RA_BF_seurat_mito_genes
  )

#NORMALIZE DATA
RA_BF_seurat <-
  NormalizeData(
    RA_BF_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
RA_BF_seurat <-
  FindVariableFeatures(RA_BF_seurat)

#X
X_counts <-
  as(
    GetAssayData(MCX20_C_seurat, slot = "counts", assay = "RNA"), 
    "sparseMatrix"
  )

#SEURAT OBJECT
X_seurat <-
  CreateSeuratObject(
    counts = X_counts,
    project = "X"
  )

X_seurat$dataset <- 
  "X"

X_seurat$region <-
  "X"

X_seurat$species <-
  "ZF"

X_seurat$cluster_orig <-
  FetchData(MCX20_C_seurat, vars = "ident")

#ADD MITOCHONDRIAL DATA
X_seurat_mito_genes <- 
  c(
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND4L",
    "ND5",
    "ND6"
  )

X_seurat[["percent.mito"]] <- 
  PercentageFeatureSet(
    X_seurat,
    features = X_seurat_mito_genes
  )

#NORMALIZE DATA
X_seurat <-
  NormalizeData(
    X_seurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000)

#PLOT GENE VARIABILITY
X_seurat <-
  FindVariableFeatures(X_seurat)

rm(
  HVC_BF_counts,
  HVC_BF_subset,
  HVC_ZF_counts,
  HVC_ZF_subset,
  MCX20_C_seurat,
  RA_BF_counts,
  RA_BF_subset,
  RA_ZF_counts,
  RA_ZF_subset,
  X_counts,
  ZFBF_seurat,
  HVC_BF,
  HVC_ZF,
  RA_BF,
  RA_ZF
)

#INTEGRATE
HVCX_seurat_anchors <-
  FindIntegrationAnchors(
    object.list = list(
      HVC_ZF_seurat,
      HVC_BF_seurat,
      RA_ZF_seurat,
      RA_BF_seurat,
      X_seurat
    )
  )

rm(
  HVC_ZF_seurat,
  HVC_BF_seurat,
  RA_ZF_seurat,
  RA_BF_seurat,
  X_seurat
)

HVCX_seurat <-
  IntegrateData(
    anchorset = HVCX_seurat_anchors
  )

DefaultAssay(HVCX_seurat) <-
  "integrated"

#SCALE DATA
HVCX_seurat_REGION <-
  ScaleData(
    HVCX_seurat,
    vars.to.regress = c(
      "nCount_RNA", "percent.mito", "region"
    ),
    model.use ="linear"
  )

#RUN PCA
HVCX_seurat_REGION <-
  RunPCA(
    HVCX_seurat_REGION,
    features = VariableFeatures(
      object = HVCX_seurat_REGION
    )
  )

#PC ELBOW PLOT
ElbowPlot(HVCX_seurat_REGION)

#CLUSTER PLOT
HVCX_seurat_REGION_res1.1 <- 
  FindNeighbors(
    HVCX_seurat_REGION,
    dims = 1:16
  )

HVCX_seurat_REGION_res1.1 <-
  FindClusters(
    HVCX_seurat_REGION_res1.1,
    resolution = 1.1,
  )

HVCX_seurat_REGION_res1.1 <-
  RunUMAP(
    HVCX_seurat_REGION_res1.1,
    dims = 1:16,
    umap.method = "umap-learn",
    metric = "correlation"
  )

HVCX_seurat_REGION_res1.1 <-
  RenameIdents(
    HVCX_seurat_REGION_res1.1,
    "0" = "1",
    "1" = "2",
    "2" = "3",
    "3" = "4",
    "4" = "5",
    "5" = "6",
    "6" = "7",
    "7" = "8",
    "8" = "9",
    "9" = "10",
    "10" = "11",
    "11" = "12",
    "12" = "13",
    "13" = "14",
    "14" = "15",
    "15" = "16",
    "16" = "17",
    "17" = "18",
    "18" = "19",
    "19" = "20",
    "20" = "21",
    "21" = "22",
    "22" = "23",
    "23" = "24",
    "24" = "25",
    "25" = "26",
    "26" = "27",
    "27" = "28",
    "28" = "29",
    "29" = "30",
    "30" = "31",
    "31" = "32",
    "32" = "33",
    "33" = "34",
    "34" = "35",
    "35" = "36",
    "36" = "37",
    "37" = "38",
    "38" = "39",
    "39" = "40",
    "40" = "41"
  )