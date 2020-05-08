library(Seurat)
library(tidyverse)
library(qs)
library(scAlign)
library(readxl)
library(data.table)

# Env variables -----------------------------------------------------------
brad_dir = Sys.getenv("BRAD_HOME")
print(brad_dir)
# Directories -------------------------------------------------------------

dir_root = file.path(brad_dir, "data2/rstudio/birds/scRNA")
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "glut")
script_name = "comparative_tasic_zf_scalign_glut_prep"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")
data_out_obj_redo_fname = file.path(out_dir, "obj_integrated_subclustered_glut_SCT.qs")

## Markers
#assay_to_use = "SCT"
#markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_gaba_%s.rds", assay_to_use))

## Compare dir
tas_data_dir = file.path(brad_dir, "/data/scrna_datasets/tasic_2018")

# Load comparison ---------------------------------------------------------

cat("Load comparison...\n")
file_prefix = c("mouse_ALM_2018-06-14", "mouse_MOp_cells_2018-10-04", "mouse_VISp_2018-06-14")

fps = file_prefix[c(1,3)]
fps_sample = c("ALM", "VISp")
tas_samp = map(seq_along(fps), function(i) {
  print(fps[i])
  res = read_csv(file.path(tas_data_dir, sprintf("%s_samples-columns.csv", fps[i])))
  res = res %>% distinct(sample_name, .keep_all=T)
  rownames(res) = paste(fps_sample[i], res$sample_name, sep="_")
  res
}) 
#cnames = Reduce(intersect, map(tas_samp, ~colnames(.x)))

tas_samp = do.call(rbind, tas_samp)
#tas_samp = read_csv(file.path(tas_data_dir, sprintf("%s_samples-columns.csv", fp)))
#tas_samp = as.data.frame(tas_samp)


tas_genes = map(fps, function(fp) {
  read_csv(file.path(tas_data_dir, sprintf("%s_genes-rows.csv", fp)))
}) %>% bind_rows()

tas_inject_md = read_xlsx(file.path(tas_data_dir, "NIHMS1001532-supplement-Supplementary_Table_6.xlsx"))

tas_samp$injection_category = tas_inject_md$injection_target_category[match(tas_samp$injection_primary, tas_inject_md$injection_target)]

fp_name = "ALM_VISp"
fname_obj = file.path(tas_data_dir, sprintf("%s_seurat_glut.qs", fp_name))
fname_avg_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_list.rds", fp_name))
fname_avg_sc_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_subclass_list.rds", fp_name))
fname_avg_ip_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_injection_target_list.rds", fp_name))


redo_obj = F
load_obj = T
if (redo_obj) {

  
  compare_objs = map(fps, function(fp) {
    print(fp)
    tas_dat = fread(file.path(tas_data_dir, sprintf("%s_exon-matrix.csv", fp)))
    tas_dat = as.data.frame(tas_dat)
    rownames(tas_dat) = tas_dat$V1
    tas_dat = tas_dat[,-1]
    rownames(tas_dat) = toupper(tas_genes$gene_symbol[match(rownames(tas_dat), tas_genes$gene_entrez_id)])
    tas_dat = as.matrix(tas_dat)
    compare_obj = CreateSeuratObject(tas_dat)
    compare_obj
  })
  
  compare_obj = merge(compare_objs[[1]], compare_objs[[2]], add.cell.ids =  fps_sample)

  compare_obj = AddMetaData(compare_obj, metadata = tas_samp)
  compare_obj = subset(compare_obj, subset=class=="Glutamatergic")
  compare_obj = SCTransform(compare_obj, vars.to.regress = c( "percent_mt_exon_reads"), 
                        return.only.var.genes = T)
  
  qsave(compare_obj, fname_obj)
} else if (load_obj) {
  compare_obj = qread(fname_obj)
}

# Process finch data ------------------------------------------------------


cat("Load finch...\n")
redo = F
if (redo) {
  obj_int_filt = qread(data_out_obj_fname)
  obj_int_filt = subset(obj_int_filt, subset=species=="zf")
  obj_int_filt = SCTransform(obj_int_filt, 
                             assay="RNA",
                             vars.to.regress = "percent.mito",
                             return.only.var.genes = F
  )
  qsave(obj_int_filt, data_out_obj_redo_fname)
} else {
  obj_int_filt = qread(data_out_obj_redo_fname)
}



# Create object ------------------------------------------------------------
cat("Create object...\n")
fname_int = file.path(out_dir, "obj_scalign.qs")

res_to_use = "cluster_int_sub2"
compare_obj$dataset = "tasic_2018"
obj_int_filt$dataset = "finch"

compare_obj$species = "mouse"

compare_obj$celltype_comb = compare_obj$cluster
print(compare_obj)
obj_int_filt$celltype_comb = FetchData(obj_int_filt, res_to_use)

objs = c(obj_int_filt, compare_obj)  
objs = map(objs, function(x) {
  DefaultAssay(x) = "SCT"
  x
})

common_features = Reduce(intersect, map(objs, ~rownames(.x)))
common_cols = c("celltype_comb", "species", "dataset")
objs = map(objs, function(x) {
  x = x[common_features, ]
  md = FetchData(x, common_cols)
  x@meta.data = md
  x
})

assay_to_use = "SCT"
var_genes = SelectIntegrationFeatures(objs, assay=rep(assay_to_use, times=length(objs)))

objs_sce = map(objs, function(x) {
  md = FetchData(x, c("species", "dataset", "celltype_comb"))
  sce = SingleCellExperiment(
    assays = list(counts = as.matrix(GetAssayData(x, assay=assay_to_use, slot="counts")[var_genes,]),
                  logcounts = as.matrix(GetAssayData(x, assay=assay_to_use, slot="data")[var_genes,])#,
                  #scale.data = GetAssayData(x, assay=assay_to_use, slot="scale.data")[var_genes,]
		  ),
    colData = DataFrame(md)
)
  
})

  sca = scAlignCreateObject(sce.objects = objs_sce,
                                  genes.use = var_genes,
                                  data.use="logcounts",
                                  cca.reduce = TRUE,
                                  ccs.compute = 40,
                                  project.name = "joint_align")

qsave(sca, fname_int)

