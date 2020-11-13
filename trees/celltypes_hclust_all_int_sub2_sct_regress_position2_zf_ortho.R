library(Seurat)
library(tidyverse)
library(qs)

library(pvclust)
library(future)
plan(sequential)
source("~/data2/rstudio/birds/utils/scRNA.R")

# Parameters ---------------------------------------------------------------

hclust_method = "complete"

# Directories -------------------------------------------------------------

## Data dir
dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
tree_dir = file.path(dev_out_dir, "trees")
script_name = "celltypes_hclust_all_int_sub2_sct_regress_position2_zf_ortho"
tree_dir = file.path(tree_dir, script_name)
tree_sub_dir = file.path(tree_dir,  sprintf("%s", hclust_method))
dir.create(tree_sub_dir, recursive = T)

data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered.qs")
data_out_avg_fname = file.path(tree_dir, "average_expr.qs")
data_out_pv_fname = file.path(tree_sub_dir, "pvclust.qs")


# Load data ---------------------------------------------------------------


res_to_use = "cluster_int_sub2"
redo = F
if (redo) {
  obj_int_filt = qread(data_fname)
  
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  idents_to_use = unique(Idents(obj_int_filt))
  
  idents_to_remove = c( "ARCO_Glut-1")
  idents_to_use = idents_to_use[!(idents_to_use %in% idents_to_remove)]
  obj_int_filt = subset(obj_int_filt, idents=idents_to_use)
  
  cells = Cells(obj_int_filt)[!is.na(obj_int_filt$position2)]
  obj_int_filt = subset(obj_int_filt, cells=cells)
  DefaultAssay(obj_int_filt) = "RNA"
  obj_int_filt = obj_int_filt %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("percent.mito", "position2", "nCount_RNA"))
  table(obj_int_filt$cluster_int_sub2)
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Marker id -----------------------------------------------------------------

assay_to_use = "RNA"
markers_fname = file.path(dev_data_dir, sprintf("marker_genes_cluster_int_sub2_%s.rds", assay_to_use))


redo_markers = T
if (redo_markers) {
  Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
  DefaultAssay(obj_int_filt) = assay_to_use
  plan(multiprocess(workers = 2))
  markers_int = FindAllMarkers(obj_int_filt, 
                               test.use = "wilcox", 
                               max.cells.per.ident = 200,
                               only.pos = F)
  saveRDS(markers_int, markers_fname)
} else {
  markers_int = readRDS(markers_fname)
}

# Average data ------------------------------------------------------------

Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
obj_int_avg = AverageExpression(obj_int_filt,assays = "RNA", slot="counts" )
obj_int_avg1 = log1p(obj_int_avg[["RNA"]])
qsave(obj_int_avg1, data_out_avg_fname)

# Cluster -----------------------------------------------------------------


nsig=50
markers_int_top = markers_int %>% 
  mutate(sign = avg_logFC>0) %>%
  group_by(cluster, sign) %>%
  top_n(nsig, abs(avg_logFC)) %>%
  ungroup() %>% 
  distinct(gene, .keep_all=T)

obj_int_avg_filt = obj_int_avg1[markers_int_top$gene,]

pv = pvclust(obj_int_avg_filt, method.dist = "correlation", method.hclust="complete", nboot = 100, parallel=4L)
plot(pv)


data_out_pv_fname2 = sub("pvclust", sprintf("pvclust_%s", nsig), data_out_pv_fname)
qsave(pv, data_out_pv_fname2)

