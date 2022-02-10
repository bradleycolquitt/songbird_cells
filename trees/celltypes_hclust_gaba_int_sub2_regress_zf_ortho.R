library(Seurat)
library(tidyverse)
library(qs)
library(pvclust)
library(future)
plan(sequential)

source("~/data2/rstudio/birds/utils/scRNA.R")


# Parameters ---------------------------------------------------------------


hclust_method = "average"


# Directories -------------------------------------------------------------


## Data dir
dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_GABA.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
tree_dir = file.path(dev_out_dir, "trees")
script_name = "celltypes_hclust_gaba_int_sub2_regress_zf_ortho"
tree_dir = file.path(tree_dir, script_name)
tree_sub_dir = file.path(tree_dir,  sprintf("%s", hclust_method))
dir.create(tree_sub_dir, recursive = T)

data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered.qs")
data_out_avg_fname = file.path(tree_dir, "average_expr.qs")
data_out_pv_fname = file.path(tree_sub_dir, "pvclust.qs")

# 


# dir_root = "~/sdd/data2/rstudio/birds/scRNA"
# dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
# dev_data_dir = file.path(dir_root, "devin_combined", "data")
# dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
# dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))
# 
# tree_dir = file.path(dev_dir, "trees")
# script_name = "celltypes_hclust_gaba_int_sub2_regress_zf_ortho"
# res_to_use = "cluster_int_sub2"
# tree_dir = file.path(tree_dir, script_name)
# 
# tree_sub_dir = file.path(tree_dir,  sprintf("%s_%s", res_to_use, hclust_method))
# dir.create(tree_sub_dir, recursive = T)
# #figures_dir = file.path(tree_dir, "figures", "all")
# #dir.create(figures_dir, recursive = T)
# 
# data_fname = file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")
# data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")
# data_out_avg_fname = file.path(tree_dir, "average_expr.qs")
# data_out_pv_fname = file.path(tree_sub_dir, "pvclust.qs")

# Load data ---------------------------------------------------------------


redo = T
if (redo) {
  obj_int_filt = qread(data_fname)
  #obj_int_filt$position2_cluster_int_sub2 = paste(obj_int_filt$position2, obj_int_filt$cluster_int_sub2, sep="_")
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  idents_to_use = Idents(obj_int_filt)
  
  ## Select GABA
  idents_to_use = idents_to_use[grepl("GABA", idents_to_use)]
#  idents_to_remove = c( "ra_GABA-Pre")
#  idents_to_use = idents_to_use[!(idents_to_use %in% idents_to_remove)]
  obj_int_filt = subset(obj_int_filt, idents=idents_to_use)
  
  DefaultAssay(obj_int_filt) = "RNA"
  obj_int_filt = obj_int_filt %>%
    NormalizeData() %>% 
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("percent.mito", "position2", "nCount_RNA"))
  
  obj_int_filt = obj_int_filt %>%
    SCTransform(vars.to.regress = c("percent.mito", "position2"))
    
  
  table(obj_int_filt$cluster_int_sub2)
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Marker id -----------------------------------------------------------------
#res_to_use = "position2_cluster_int_sub2"
assay_to_use = "SCT"

markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_%s_gaba_%s.rds", res_to_use, assay_to_use))
#markers_pos_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_gaba_%s_pos.rds", assay_to_use))

redo_markers = T
if (redo_markers) {
  Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
  DefaultAssay(obj_int_filt) = assay_to_use
  plan(multiprocess(workers = 2))
  markers_int = FindAllMarkers(obj_int_filt, 
                               test.use = "wilcox", 
                               #logfc.threshold = log(2), 
                               #min.pct = .5,
                               #min.diff.pct = .25,
                               max.cells.per.ident = 200,
                               only.pos = F)
  saveRDS(markers_int, markers_fname)
  
  # ## Within position
  # markers_int_pos = map(c("hvc", "ra"), function(pos) {
  #   cells = colnames(obj_int_filt)[obj_int_filt$position2==pos]
  #   tmp = subset(obj_int_filt, cells=cells)
  #   FindAllMarkers(tmp, 
  #                  test.use = "wilcox", 
  #                  #logfc.threshold = log(2), 
  #                  #min.pct = .5,
  #                  #min.diff.pct = .25,
  #                  max.cells.per.ident = 200,
  #                  only.pos = F)
  # })
  #saveRDS(markers_int_pos, markers_pos_fname)
} else {
  markers_int = readRDS(markers_fname)
  #markers_int_pos = readRDS(markers_pos_fname)
}

markers_export_fname = file.path(tree_dir, sprintf("marker_genes_%s_gaba_%s.txt", res_to_use, assay_to_use))
write.table(markers_int, markers_export_fname, sep="\t", row.names=F, col.names=T, quote=F)

Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]


# Average data ------------------------------------------------------------


obj_int_avg = AverageExpression(obj_int_filt, assays = c("RNA", 'SCT'), slot="counts")
obj_int_avg1 = log1p(obj_int_avg[["RNA"]])
qsave(obj_int_avg1, data_out_avg_fname)


# Cluster -----------------------------------------------------------------


nsig = 50
markers_int_top = markers_int %>% 
  mutate(sign = avg_logFC>0) %>%
  #group_by(cluster, sign) %>%
  group_by(cluster) %>% 
  #top_n(-1 * nsig, p_val_adj) %>%
  #top_n(nsig, abs(avg_logFC)) %>%
  top_n(nsig, avg_logFC) %>%
  distinct(gene, .keep_all=T)

obj_int_avg_filt = obj_int_avg1[markers_int_top$gene,]


pv = pvclust(obj_int_avg_filt, method.hclust=hclust_method, nboot = 100, parallel=4L)
plot(pv)
qsave(pv, data_out_pv_fname)

