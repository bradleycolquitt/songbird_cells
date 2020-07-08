library(Seurat)
library(tidyverse)
library(qs)
library(SCENIC)

source("~/data2/rstudio/birds/utils/go.R")
#np = import("np")

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "grn")
#script_data_dir = file.path(out_dir, "celltypes_hclust_glut_int_sub2")
#script_data_sub_dir = file.path(script_data_dir, sub_dir)
script_name = "export_to_numpy_glut"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")


# Load data ---------------------------------------------------------------

tfs = get_tf_genes_human()

obj_int_filt = qread(data_out_obj_fname)
obj_int_filt = subset(obj_int_filt, subset=species=="zf")
obj_int_filt = FindVariableFeatures(obj_int_filt, nfeatures=3000)

to_use = c(VariableFeatures(obj_int_filt))
exprMat = GetAssayData(obj_int_filt, assay="SCT", slot="counts")
exprMat = exprMat[to_use,]
exprMat_t <- t(exprMat)

inputTFs = intersect(tfs$external_gene_name, to_use)
write.table(inputTFs, file = file.path(out_dir, "1.1_inputTFs.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(exprMat_t, file = file.path(out_dir, "1.1_exprMatrix_filtered_t.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
