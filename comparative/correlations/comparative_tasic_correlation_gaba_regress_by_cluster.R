


library(Seurat)
library(qs)
library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)
library(rtracklayer)
library(proxy)

library(ggsci)
library(cowplot)
library(viridis)
library(pheatmap)
library(ggcorrplot)
theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "anchor5_filter200_score_30")
script_name = "comparative_tasic_correlation_gaba_regress_by_cluster"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_gaba_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")

## Markers
assay_to_use = "SCT"
markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_position2_cluster_int_sub2_gaba_%s.rds", assay_to_use))

# Load comparison data ----------------------------------------------------

tas_data_dir = "~/data/scrna_datasets/tasic_2018"
file_prefix = c("mouse_ALM_2018-06-14", "mouse_MOp_cells_2018-10-04", "mouse_VISp_2018-06-14")

fp = file_prefix[1]

tas_samp = read_csv(file.path(tas_data_dir, sprintf("%s_samples-columns.csv", fp)))
tas_samp = as.data.frame(tas_samp)
rownames(tas_samp) = tas_samp$sample_name

tas_genes = read_csv(file.path(tas_data_dir, sprintf("%s_genes-rows.csv", fp)))

tas_inject_md = read_xlsx(file.path(tas_data_dir, "NIHMS1001532-supplement-Supplementary_Table_6.xlsx"))

tas_samp$injection_category = tas_inject_md$injection_target_category[match(tas_samp$injection_primary, tas_inject_md$injection_target)]

fname_obj = file.path(tas_data_dir, sprintf("%s_seurat.rds", fp))
fname_avg_list = file.path(tas_data_dir, sprintf("%s_seurat_data_avg_list.rds", fp))
fname_avg_sc_list = file.path(tas_data_dir, sprintf("%s_seurat_data_avg_subclass_list.rds", fp))

redo_obj = F
load_obj = T
if (redo_obj) {

  
  tas_dat = fread(file.path(tas_data_dir, sprintf("%s_exon-matrix.csv", fp)))
  tas_dat = as.data.frame(tas_dat)
  rownames(tas_dat) = tas_dat$V1
  tas_dat = tas_dat[,-1]
  tas_dat = as.matrix(tas_dat)
  
  tas_obj = CreateSeuratObject(tas_dat)
  tas_obj = AddMetaData(tas_obj, metadata = tas_samp)
  tas_obj = SCTransform(tas_obj, vars.to.regress = c( "percent_mt_exon_reads"), return.only.var.genes = F)
  tas_obj = RunALRA(tas_obj)
  
  
  saveRDS(tas_obj, fname_obj)
} else if (load_obj) {
  tas_obj = readRDS(fname_obj)
}

redo_avg = F
if (redo_avg) {
  Idents(tas_obj) = tas_obj@meta.data$cluster
  tas_obj_avg1 = AverageExpression(tas_obj, assays = c("SCT"), slot="data")
  
  tas_obj_avg1 = map(tas_obj_avg1, function(x) {
    
    rownames(x) = tas_genes$gene_symbol[match(rownames(x), tas_genes$gene_entrez_id)]
    rownames(x) = toupper(rownames(x))
    
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
} else {
  tas_obj_avg = readRDS(fname_avg_list)
}
tas_obj_avg_sct = tas_obj_avg[["SCT"]]

# Load finch data ---------------------------------------------------------

obj_int_filt = qread(data_out_obj_fname)
obj_int_filt_markers = readRDS(markers_fname)

res_to_use = "cluster_int_sub2"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]

# By Cluster --------------------------------------------------------------

tas_inh_markers_fname = file.path(tas_data_dir, "tas_obj_inh_markers.rds")
redo_tas_markers = F
if (redo_tas_markers) {
  tas_obj_inh = subset(tas_obj, subset=class=="GABAergic")
  DefaultAssay(tas_obj_inh) = "SCT"
  tas_obj_inh_cur = subset(tas_obj_inh, subset=brain_region=="ALM")
  Idents(tas_obj_inh_cur) = tas_obj_inh_cur$cluster
  genes_to_test = na.omit(tas_genes$gene_entrez_id[match(rownames(obj_int_filt), toupper(tas_genes$gene_symbol))])
  genes_to_test = as.character(genes_to_test)
  genes_to_test = genes_to_test[genes_to_test %in% rownames(tas_obj_inh_cur)]
  tas_inh_markers = FindAllMarkers(tas_obj_inh_cur, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, only.pos=T)
  tas_inh_markers = tas_inh_markers %>%
    rename(entrez_id = gene) %>%
    mutate(gene = tas_genes$gene_symbol[match(entrez_id, tas_genes$gene_entrez_id)]) 
  tas_inh_markers = tas_inh_markers %>%
    mutate(gene = toupper(gene))
  saveRDS(tas_inh_markers, tas_inh_markers_fname)
} else {
  tas_inh_markers = readRDS(tas_inh_markers_fname)
}


# Calc specificity --------------------------------------------------------

mm_mat_filt = tas_obj_avg_sct[rownames(tas_obj_avg_sct) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]


obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]
obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("Pre", colnames(obj_int_filt_avg_filt))]


## Filter for markers or variable genes 
ngenes = 200
tas_inh_markers_sig = tas_inh_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  filter(p_val_adj < .05) %>% 
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC) %>%
  filter(gene %in% tas_inh_markers_sig$gene)

obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% obj_int_filt_markers_to_use$gene,]

mm_md_n = tas_samp %>% filter(grepl("GABAergic", class))
mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$cluster]
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

ann_col = mm_md_n %>% filter(cluster %in% colnames(mm_mat_filt)) %>%
  distinct(cluster, .keep_all=T) %>% select(cluster, subclass) %>% as.data.frame()
rownames(ann_col) = ann_col$cluster
ann_col = ann_col %>% select(-cluster)

classes = rev(c("Pvalb", "Sst", "Vip", "Sncg", "Serpinf1", "Lamp5", "Meis2"))
class_colors = pal_locuszoom()(length(classes))
names(class_colors) = classes
obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

obj_cor2 = obj_cor1
obj_cor2[obj_cor_sig_adj>.05] = 0

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")
pheatmap(obj_cor1, annotation_row = ann_col,show_colnames=T, 
         cluster_cols=hc_c, 
         cluster_rows = hc_r,
         annotation_names_col = T,
         treeheight_col = 0,
         treeheight_row = 8,
         color=colors, 
         breaks = breaks)
pheatmap(obj_cor1, annotation_row = ann_col, 
         show_colnames=T,  color=colors, 
         cluster_cols=hc_c, 
         cluster_rows = hc_r,
         treeheight_col = 0,
         treeheight_row = 8,
         border_color = NA,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap.pdf"), height=10, width=5)

pheatmap(t(obj_cor1), annotation_col = ann_col, show_rownames=T,  color=colors, 
         cluster_rows = hc_c, 
         cluster_cols = hc_r,
         treeheight_row = 0,
         treeheight_col = 10,
         border_color = NA,
         annotation_colors = list(subclass = class_colors),
         display_numbers = obj_cor_sig_adj_label,
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_horiz.pdf"), height=3.5, width=10)


# Network -----------------------------------------------------------------

n_top_cor = 4
cor_df = melt(obj_cor1)
colnames(cor_df) = c("target", "source", "value")
cor_df = cor_df %>% group_by(source) %>% 
  top_n(n_top_cor, value) %>%
  filter(value >= .2)

obj_cor_sig_adj_df1 = obj_cor_sig_adj_df %>% rename(p_val_adj = value,
                                                    target = Var2,
                                                    source = Var1)

cor_df = cor_df %>% inner_join(obj_cor_sig_adj_df1)
write.table(cor_df, file=file.path(out_dir, sprintf("cors_top%s.txt", n_top_cor)), sep="\t", row.names=F, quote=F)

md = data.frame(name = c(colnames(obj_cor1), rownames(obj_cor1)), 
                species = c(rep("mouse", length(colnames(obj_cor1))), 
                            rep("finch", length(rownames(obj_cor1))))
)
md = md %>%
  mutate(cluster = name) %>% 
  left_join(mm_md_n %>% distinct(cluster, cluster, subclass, brain_subregion, brain_region)) %>%
  rename(Region = subclass)
write.table(md, file=file.path(out_dir, "nodes.txt"), sep="\t", row.names=F, quote=F)
