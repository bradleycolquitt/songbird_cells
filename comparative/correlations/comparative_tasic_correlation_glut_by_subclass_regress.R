
source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
library(viridis)
library(pheatmap)
library(ggcorrplot)
library(Seurat)
library(qs)
library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)
library(rtracklayer)
library(ggsci)
library(proxy)


library(cowplot)
theme_set(theme_cowplot())
# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "anchor5_filter200_score_30")
script_name = "comparative_tasic_correlation_glut_by_subclass_regress"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")

## Markers
assay_to_use = "SCT"
markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))

# Load comparison data ----------------------------------------------------

tas_data_dir = "~/data/scrna_datasets/tasic_2018"
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

tas_samp = do.call(rbind, tas_samp)

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
load_obj = F
if (redo_obj) {

  
  tas_objs = map(fps, function(fp) {
    print(fp)
    fread(file.path(tas_data_dir, sprintf("%s_exon-matrix.csv", fp)))
    tas_dat = as.data.frame(tas_dat)
    rownames(tas_dat) = tas_dat$V1
    tas_dat = tas_dat[,-1]
    tas_dat = as.matrix(tas_dat)
    tas_obj = CreateSeuratObject(tas_dat)
    tas_obj
  })
  
  tas_obj = merge(tas_objs[[1]], tas_objs[[2]], add.cell.ids =  fps_sample)

  tas_obj = AddMetaData(tas_obj, metadata = tas_samp)
  tas_obj = subset(tas_obj, subset=class=="Glutamatergic")
  tas_obj = SCTransform(tas_obj, vars.to.regress = c( "percent_mt_exon_reads"), 
                        return.only.var.genes = T)
  
  qsave(tas_obj, fname_obj)
} else if (load_obj) {
  tas_obj = qread(fname_obj)
}


redo_avg = F
if (redo_avg) {
  Idents(tas_obj) = tas_obj@meta.data$cluster
  tas_obj_avg1 = AverageExpression(tas_obj, assays = c("SCT"), slot="data")
  #tas_obj_avg_alra = tas_obj_avg1[["alra"]]
  
  tas_obj_avg1 = map(tas_obj_avg1, function(x) {
    
    rownames(x) = tas_genes$gene_symbol[match(rownames(x), tas_genes$gene_entrez_id)]
    rownames(x) = toupper(rownames(x))
    x
  })
  
  saveRDS(tas_obj_avg1, fname_avg_list)
  Idents(tas_obj) = tas_obj@meta.data$subclass
  tas_obj_avg_sc = AverageExpression(tas_obj, assays = c("SCT"), slot="data")
  
  tas_obj_avg_sc = map(tas_obj_avg_sc, function(x) {
    
    rownames(x) = tas_genes$gene_symbol[match(rownames(x), tas_genes$gene_entrez_id)]
    rownames(x) = toupper(rownames(x))
    x
  })
  
  saveRDS(tas_obj_avg_sc, fname_avg_sc_list)
  
  ## Injection_primary
  Idents(tas_obj) = tas_obj@meta.data$injection_primary
  tas_obj_cur = subset(tas_obj, subset=injection_primary!="No Injection")
  
  tas_obj_avg_ip = AverageExpression(tas_obj_cur, assays = c("SCT"), slot="data")
  
  tas_obj_avg_ip = map(tas_obj_avg_ip, function(x) {
    rownames(x) = tas_genes$gene_symbol[match(rownames(x), tas_genes$gene_entrez_id)]
    rownames(x) = toupper(rownames(x))
    x
  })
  
  saveRDS(tas_obj_avg_ip, fname_avg_ip_list)
  
} else {
  tas_obj_avg1 = readRDS(fname_avg_list)
  tas_obj_avg_sc = readRDS(fname_avg_sc_list)
}
tas_obj_avg_sct = tas_obj_avg1[["SCT"]]
tas_obj_avg_sc_sct = tas_obj_avg_sc[["SCT"]]


# Load finch data ---------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)
obj_int_filt_markers = readRDS(markers_fname)

res_to_use = "cluster_int_sub2"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]



# By Cluster --------------------------------------------------------------

## Camparison markers
tas_exc_markers_fname = file.path(tas_data_dir, "tas_obj_exc_subclass_markers.rds")
redo_tas_markers = F
if (redo_tas_markers) {
  Idents(tas_obj) = tas_obj$subclass
  genes_to_test = na.omit(tas_genes$gene_entrez_id[match(rownames(obj_int_filt), toupper(tas_genes$gene_symbol))])
  genes_to_test = as.character(genes_to_test)
  genes_to_test = genes_to_test[genes_to_test %in% rownames(tas_obj)]
  tas_exc_markers = FindAllMarkers(tas_obj, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, only.pos=T,max.cells.per.ident = 200
                                   )
  tas_exc_markers = tas_exc_markers %>%
    rename(entrez_id = gene) %>%
    mutate(gene = tas_genes$gene_symbol[match(entrez_id, tas_genes$gene_entrez_id)]) 
  tas_exc_markers = tas_exc_markers %>%
    mutate(gene = toupper(gene))
  saveRDS(tas_exc_markers, tas_exc_markers_fname)
} else {
  tas_exc_markers = readRDS(tas_exc_markers_fname)
}


# Calc specificity --------------------------------------------------------

mm_mat_filt = tas_obj_avg_sc_sct[rownames(tas_obj_avg_sc_sct) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
mm_mat_filt = mm_mat_filt[,!grepl("CR|NP", colnames(mm_mat_filt))]

obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]

## Filter for markers or variable genes 
ngenes = 200
tas_inh_markers_sig = tas_exc_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  filter(cluster %in% colnames(obj_int_filt_avg_filt)) %>%
  group_by(cluster) %>%
  top_n(-1*ngenes, p_val_adj) %>%
  top_n(ngenes, avg_logFC) %>%
  filter(gene %in% tas_exc_markers$gene)
obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% obj_int_filt_markers_to_use$gene,]

mm_md_n = tas_samp %>% filter(grepl("Glutamatergic", class))
mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$subclass]
mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]

mat_a = log1p(obj_int_filt_avg_filt) + .1
mat_b = log1p(mm_mat_filt) + .1
obj_cor = specificity_correlate(mat_a, mat_b, method="spearman")

# Shuffle -----------------------------------------------------------------


obj_cor_shuf_fname = file.path(out_dir, "obj_cor_shuf.qs")
redo = T
if (redo) {
  nrep = 100
  plan(multiprocess(workers=12, gc=T))
  obj_cor_shuf = future_map(1:nrep, function(i) {
    mat_a_cur = mat_a
    mat_b_cur = mat_b
    rownames(mat_a_cur) = sample(rownames(mat_a_cur))
    obj_cor = specificity_correlate(mat_a_cur, mat_b_cur)
    
    
    obj_cor_df = melt(obj_cor)
    colnames(obj_cor_df) =  c("celltype", "compare", "value")
    obj_cor_df
  }) %>% bind_rows() 
  qsave(obj_cor_shuf, obj_cor_shuf_fname)
} else {
  obj_cor_shuf = qread(obj_cor_shuf_fname)
}

obj_cor_shuf_stat = obj_cor_shuf %>% group_by(celltype, compare) %>%
  summarize(value_mean = mean(value),
            value_sd = sd(value),
            value_q99 = quantile(value, .99),
            value_q95 = quantile(value, .95))

obj_cor_df = melt(obj_cor)
colnames(obj_cor_df) = c("celltype", "compare", "value")
obj_cor_df = obj_cor_df %>% left_join(obj_cor_shuf_stat)

sig_thresh = .05
obj_cor_df_filt = obj_cor_df %>% filter(value > value_q95)
obj_cor_filt = obj_cor[,colnames(obj_cor) %in% obj_cor_df_filt$compare]

obj_cor_df = obj_cor_df %>% mutate(sig_label = if_else(value>value_q95, "*", ""))
obj_cor_sig_adj_label = acast(obj_cor_df, celltype~compare, value.var="sig_label")




# Plot heatmap ------------------------------------------------------------

obj_cor1 = t(obj_cor_filt)
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks = seq(-0.25, .25, length.out=(palette_length+1))
hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

pheatmap(obj_cor1, 
         annotation_names_row = F,
         annotation_names_col = F,
         show_colnames=T, 
         cluster_cols=hc_c, 
         cluster_rows = hc_r, 
         treeheight_row = 8,
         treeheight_col = 8, 
         color=colors, 
         border_color = NA,
         breaks=breaks,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_vert.pdf"), height=3, width=3)

ct_order = c("L2/3 IT", "L4 IT","L5 IT", "L5 PT", "L6 IT", "L6 CT", "L6b")
obj_cor1 = obj_cor1[ct_order,]
obj_cor_sig_adj_label = obj_cor_sig_adj_label[,ct_order]
pheatmap(obj_cor1, 

         annotation_names_row = F,
         annotation_names_col = F,
         show_colnames=T, 
         cluster_cols=hc_c, 
         cluster_rows = F,
         treeheight_row = 8,
         treeheight_col = 8, 
         color=colors, 
         border_color = NA,
         breaks=breaks,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_clust.pdf"), height=3, width=3)

ct_order = c("L2/3 IT", "L4","L5 IT", "L5 PT", "L6 IT", "L6 CT", "L6b")
obj_cor1 = obj_cor1[ct_order,]
obj_cor_sig_adj_label = obj_cor_sig_adj_label[,ct_order]

ct_order2 = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-3", "HVC_Glut-5", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
obj_cor1 = obj_cor1[,ct_order2]
obj_cor_sig_adj_label = obj_cor_sig_adj_label[ct_order2,]

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0 

colors =  colorRampPalette(brewer_pal(palette="Blues")(9))(palette_length)

breaks =  seq(0, max(obj_cor1_floor), length.out=palette_length+1)
pheatmap(obj_cor1_floor, 
         annotation_names_row = F,
         annotation_names_col = F,
         show_colnames=T, 
         cluster_cols= F, 
         cluster_rows = F,
         treeheight_row = 8,
         treeheight_col = 8, 
         color=colors, 
         border_color = NA,
         breaks=breaks,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_clust2.pdf"), height=3, width=3)


pheatmap(obj_cor1, 
         annotation_names_row = F,
         annotation_names_col = F,
         annotation_legend=F,
         show_colnames=F, 
         show_rownames=F,
         cluster_cols=hc_c, 
         cluster_rows = hc_r, 
         treeheight_row = 8,
         treeheight_col = 8, 
         color=colors, 
         border_color = NA,
         breaks=breaks,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_names.pdf"), height=4, width=2)

