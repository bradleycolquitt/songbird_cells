library(Seurat)
library(tidyverse)
library(qs)
library(pvclust)
library(future)

source("~/data2/rstudio/birds/utils/scRNA.R")

# Parameters ---------------------------------------------------------------

hclust_method = "average"


# Directories -------------------------------------------------------------

dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

tree_dir = file.path(dev_dir, "trees")
script_name = "celltypes_hclust_all_int_sub2_sct_regress_position2"
tree_dir = file.path(tree_dir, script_name)
tree_sub_dir = file.path(tree_dir,  sprintf("%s", hclust_method))
dir.create(tree_sub_dir, recursive = T)

data_fname= file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered.qs")
data_out_avg_fname = file.path(tree_dir, "average_expr.qs")
data_out_pv_fname = file.path(tree_sub_dir, "pvclust.qs")


# Load data ---------------------------------------------------------------

res_to_use = "cluster_int_sub2"
redo = F
if (redo) {
  obj_int_filt = qread(data_fname)

  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  idents_to_use = Idents(obj_int_filt)
  
  idents_to_remove = c( "ARCO_Glut-1")
  idents_to_use = idents_to_use[!(idents_to_use %in% idents_to_remove)]
  obj_int_filt = subset(obj_int_filt, idents=idents_to_use)
  
  obj_int_filt = obj_int_filt %>%
    SCTransform(vars.to.regress=c("percent.mito","position2")) 
  
  plan(multiprocess(workers=10, gc=T))
  DefaultAssay(obj_int_filt) = "RNA"
  obj_int_filt = obj_int_filt %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("percent.mito", "position2", "nCount_RNA"))
  
  md = FetchData(obj_int_filt, c("cluster_int_sub2")) %>%
    rownames_to_column() %>%
    mutate(cluster_int_sub2 = if_else(cluster_int_sub2=="VLMC", "Mural", cluster_int_sub2)) %>%
    column_to_rownames()
  obj_int_filt = AddMetaData(obj_int_filt, md)
  table(obj_int_filt$cluster_int_sub2)
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Marker id -----------------------------------------------------------------

assay_to_use = "RNA"
markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_%s.rds", assay_to_use))

redo_markers = F
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

Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]


# Average data ------------------------------------------------------------


Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
obj_int_avg = AverageExpression(obj_int_filt,assays = "RNA", slot="counts" )
obj_int_avg1 = log1p(obj_int_avg[["RNA"]])
qsave(obj_int_avg1, data_out_avg_fname)

# Cluster -----------------------------------------------------------------

nsig=50
markers_int_top = markers_int %>% 
  group_by(cluster) %>% 
  top_n(nsig, avg_logFC) %>% 
  ungroup() %>% 
  distinct(gene, .keep_all=T)

obj_int_avg_filt = obj_int_avg1[markers_int_top$gene,]

pv = pvclust(obj_int_avg_filt, method.dist = "correlation", method.hclust=hclust_method, nboot = 100, parallel=4L)
plot(pv)


data_out_pv_fname2 = sub("pvclust", sprintf("pvclust_%s", nsig), data_out_pv_fname)
qsave(pv, data_out_pv_fname2)

