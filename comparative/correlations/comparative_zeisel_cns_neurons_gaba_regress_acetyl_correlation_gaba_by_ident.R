
source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
library(viridis)
library(pheatmap)

library(ggcorrplot)
library(Seurat)

library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)
library(rtracklayer)
library(ggsci)
library(proxy)
library(qs)
library(scales)
library(colorspace)
library(future)


library(cowplot)
theme_set(theme_cowplot())
#here = here::here

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "anchor5_filter200_score_30")
script_name = "comparative_zeisel_cns_neurons_gaba_regress_acetyl_forebrain_correlation_gaba_by_ident"
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

## Comparison 
loom_data_dir = c("~/data/scrna_datasets/zeisel_2018/")
compare_prefix = "l6_r2_cns_neurons_seurat"

# Load finch data ---------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)
obj_int_filt_markers = readRDS(markers_fname)

res_to_use = "cluster_int_sub2"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]


# Load comparison data ----------------------------------------------------

compare_prefix1 = "l6_r2_cns_neurons_seurat_GABA_Acetyl_forebrain"


compare_sub_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix1))

redo = F
if (redo) {
  obj_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix))
  zei_obj_in = qread(obj_fname)
  zei_obj_in$idents = Idents(zei_obj_in)
  nt_select = "GABA|Acetyl"
  cells = Cells(zei_obj_in)[grepl(nt_select, zei_obj_in$Neurotransmitter)]
  

  zei_obj_in1 = subset(zei_obj_in, cells=cells)
  
  regions_to_use = c("Hippocampus,Cortex",
                     "Hippocampus",
                     "Olfactory bulb",
                     "Striatum dorsal",
                     "Striatum ventral",
                     "Striatum dorsal, Striatum ventral",
                     "Striatum dorsal, Striatum ventral,Amygdala",
                     "Pallidum")
  cells = Cells(zei_obj_in1)[zei_obj_in1$Region %in% regions_to_use]
  zei_obj_in1 = subset(zei_obj_in1, cells=cells)
  rm(zei_obj_in)
  
  options(future.globals.maxSize=10000 * 1024^2)
  plan(multiprocess(workers=3))
  zei_obj_in1 = SCTransform(zei_obj_in1, vars.to.regress = c("MitoRiboRatio"))
  qsave(zei_obj_in1, compare_sub_fname)
} else {
  zei_obj_in1 = qread(compare_sub_fname)
}

# Average -----------------------------------------------------------------


avg_list_fname = file.path(loom_data_dir, sprintf("%s_ident_sct_avg.rds", compare_prefix1))

redo_avg = F
if (redo_avg) {
  #Idents(zei_obj_in1) = zei_obj_in1$ident
  zei_obj_in_avg = AverageExpression(zei_obj_in1, assays = "SCT")
  qsave(zei_obj_in_avg, avg_list_fname)
} else {
  zei_obj_in_avg = qread(avg_list_fname)
}

zei_obj_in_avg1 = zei_obj_in_avg[["SCT"]]



# Compare marker ID --------------------------------------------------------------

compare_markers_fname = file.path(loom_data_dir, sprintf("%s_markers.rds", compare_prefix1))
redo_compare_markers = F
if (redo_compare_markers) {
  genes_to_test = intersect(rownames(obj_int_filt), rownames(zei_obj_in1))

  compare_markers = FindAllMarkers(zei_obj_in1, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, 
                                   max.cells.per.ident = 200,
                                   only.pos=T)

  saveRDS(compare_markers, compare_markers_fname)
} else {
  compare_markers = readRDS(compare_markers_fname)
}


# Calc specificity --------------------------------------------------------

mm_mat_filt = zei_obj_in_avg1[rownames(zei_obj_in_avg1) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
#mm_mat_filt = mm_mat_filt[,!grepl("CR", colnames(mm_mat_filt))]
#mm_mat_filt = mm_mat_filt[,grep("ALM", colnames(mm_mat_filt))]

obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]
obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("Pre", colnames(obj_int_filt_avg_filt))]
#obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("IEG", colnames(obj_int_filt_avg_filt))]

## Filter for markers or variable genes 
ngenes = 200
compare_markers_sig = compare_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  #top_n(-1*ngenes, p_val_adj) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  #filter(!(gene %in% tx_mt$gene_name)) %>%
  #filter(cluster %in% colnames(obj_int_filt_avg_filt)) #%>%
  filter(p_val_adj < .05) %>% 
  group_by(cluster) %>%
  #top_n(-200, p_val_adj) %>%
  top_n(ngenes, avg_logFC) %>%
  filter(gene %in% compare_markers_sig$gene)

obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% obj_int_filt_markers_to_use$gene,]

mm_md_n = zei_obj_in1@meta.data
mm_md_n$idents = Idents(zei_obj_in1)
mm_md_n = mm_md_n %>%
  filter(Region %in% c("Hippocampus,Cortex",
                       "Hippocampus",
                       "Olfactory bulb",
                       "Striatum dorsal",
                       "Striatum ventral",
                       "Striatum dorsal, Striatum ventral",
                       "Striatum dorsal, Striatum ventral,Amygdala",
                       "Pallidum"))
mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$ident]
#mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_avg_filt), rownames(mm_mat_filt))),]
mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]



mat_a = log1p(obj_int_filt_avg_filt) + .1
mat_b = log1p(mm_mat_filt) + .1
obj_cor = specificity_correlate(mat_a, mat_b, method="spearman")

obj_int_filt_avg_filt_log_spec = calc_specificity(mat_a)
mm_mat_filt_log_spec = na.omit(calc_specificity(mat_b))
obj_int_filt_avg_filt_log_spec = obj_int_filt_avg_filt_log_spec[match(rownames(mm_mat_filt_log_spec), rownames(obj_int_filt_avg_filt_log_spec)),]

obj_cor_sig_adj =cor_sig_matrix(obj_int_filt_avg_filt_log_spec, mm_mat_filt_log_spec, method="spearman", symmetric = F)

sig_thresh = .05
obj_cor_sig_adj_df = melt(obj_cor_sig_adj) %>% filter(value < sig_thresh)
obj_cor_filt = obj_cor[,colnames(obj_cor) %in% obj_cor_sig_adj_df$Var2]

obj_cor_sig_adj_label = matrix("", nrow=nrow(obj_cor_sig_adj), ncol=ncol(obj_cor_sig_adj), dimnames = dimnames(obj_cor))
obj_cor_sig_adj_label[obj_cor_sig_adj<sig_thresh] = "*"
obj_cor_sig_adj_label = obj_cor_sig_adj_label[, colnames(obj_cor_sig_adj_label) %in% colnames(obj_cor_filt)]


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
#obj_cor_sig_adj_label = matrix("", nrow=nrow(obj_cor_filt), ncol=ncol(obj_cor_filt), 
# dimnames = dimnames(obj_cor_filt))
#obj_cor_sig_adj_label[obj_cor_sig_adj<sig_thresh] = "*"
#obj_cor_sig_adj_label = obj_cor_sig_adj_label[, colnames(obj_cor_sig_adj_label) %in% colnames(obj_cor_filt)]



# Plot heatmap ------------------------------------------------------------


mm_md_n = mm_md_n %>%
  mutate(Region = as.character(Region)) %>% 
  mutate(Region = case_when(Region == "Striatum dorsal, Striatum ventral" ~ "Septum",
                            Region == "Striatum dorsal" ~ "Dorsal striatum",
                            Region == "Striatum ventral" ~ "Ventral striatum",
                            TRUE ~ Region))

ann_col = mm_md_n %>%
  filter(idents %in% colnames(mm_mat_filt)) %>%
  distinct(idents, .keep_all=T) %>% 
  select(idents, Taxonomy_group, Region) %>%
  column_to_rownames(var="idents") %>% 
  as.data.frame()


regions = na.omit(unique(mm_md_n$Region))
region_colors = qualitative_hcl(length(regions))
names(region_colors) = regions

obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label = obj_cor_sig_adj_label[,match(rownames(obj_cor1), colnames(obj_cor_sig_adj_label))]
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

pheatmap(obj_cor1, 
         cluster_rows = hc_r,
         cluster_cols = hc_c,
         annotation_row = ann_col, 
         annotation_colors = list(Region = region_colors),
         annotation_legend = F,
         show_colnames = T,  
         border_color = NA,
         annotation_names_row = F,
         annotation_names_col = F,
         color=colors, breaks=breaks,
         treeheight_row = 8, 
         treeheight_col = 5,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "idents_markers_specificity_heatmap_no_legend.pdf"),
         height=7, width=3)

pheatmap(obj_cor1, 
         cluster_rows = hc_r,
         cluster_cols = hc_c,
         annotation_row = ann_col, 
         annotation_colors = list(Region = region_colors),
         annotation_legend = T,
         show_colnames = T,  
         border_color = NA,
         annotation_names_row = F,
         annotation_names_col = F,
         color=colors, breaks=breaks,
         treeheight_row = 8, 
         treeheight_col = 5,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "idents_markers_specificity_heatmap_with_legend.pdf"),
         height=4, width=10)



# Plot heatmap, 2 ---------------------------------------------------------

mm_md_n = mm_md_n %>%
  mutate(Region = as.character(Region)) %>% 
  mutate(Region = case_when(Region == "Striatum dorsal, Striatum ventral" ~ "Septum",
                            Region == "Striatum dorsal" ~ "Dorsal striatum",
                            Region == "Striatum ventral" ~ "Ventral striatum",
                            TRUE ~ Region))

ann_col = mm_md_n %>%
  filter(idents %in% colnames(mm_mat_filt)) %>%
  distinct(idents, .keep_all=T) %>% 
  select(idents, Taxonomy_group, Region) %>%
  column_to_rownames(var="idents") %>% 
  as.data.frame()


regions = na.omit(unique(mm_md_n$Region))
region_colors = qualitative_hcl(length(regions))
names(region_colors) = regions

obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label = obj_cor_sig_adj_label[,match(rownames(obj_cor1), colnames(obj_cor_sig_adj_label))]
#palette_length = 50
#breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
#            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

breaks = c(min(obj_cor1), 0, max(obj_cor1))
colors =  colorRamp2(breaks =breaks,  colors= c("navy", "white", "firebrick3"))


hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

star_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", t(obj_cor_sig_adj_label)[i, j]), x, y, gp = gpar(fontsize = 6))}

target_width = 1
target_height = target_width * nrow(obj_cor1) / ncol(obj_cor1) 

ann_col = ann_col[rownames(obj_cor1), ]
regions = unique(ann_col$Region)
region_col = brewer_pal(type="qual")(length(regions))
names(region_col) = regions

tax = unique(ann_col$Taxonomy_group)
tax_col = brewer_pal(type="qual", palette = 2)(length(tax))
names(tax_col) = tax

ht_opt(legend_labels_gp = gpar(fontsize=6),
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(obj_cor1), "in"),
       legend_grid_height = unit(target_height / nrow(obj_cor1), "in"))
ha = rowAnnotation(df = ann_col, 
                   annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                   col = list(Region = region_col,
                              Taxonomy_group = tax_col),
                   annotation_width=unit(target_width*.2, "in"),
                   simple_anno_size_adjust=T,
                   annotation_legend_param = list(
                     Region = list(gp = gpar(fontsize=6)),
                     Taxonomy_group = list(gp = gpar(fontsize=6))
                   ))

hm = Heatmap(obj_cor1, col=colors, name = "rho",
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             clustering_method_rows = "ward.D",
             clustering_method_columns = "ward.D",
             show_row_dend = T,
             show_column_dend = T,
             cell_fun = star_fun,
             row_names_gp = gpar(fontsize=6),
             column_names_gp = gpar(fontsize=6),
             left_annotation = ha

             
)

pdf(file.path(out_dir, "idents_markers_specificity_heatmap.pdf"), 
height=target_height+4, 
width=target_width+6)
hm
dev.off()
hm

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
  mutate(idents = name) %>% 
  left_join(mm_md_n %>% distinct(idents, Description, Region))
write.table(md, file=file.path(out_dir, "nodes.txt"), sep="\t", row.names=F, quote=F)

# Marker genes ------------------------------------------------------------



markers_lake = data.frame(gene = c("LAMP5", "CUX2", "GLRA3", "CARTPT", "PRSS12", "RORB", "TPBG", "TOX", "FOXP2", "ETV1", "RPRM","RXFP1", "TLE4", "GRIK4", "NTNG2", "OPRK1", "NR4A2", "ADRA2A"),
                          layer = c("L2/3", "L2/3/4", "L2", "L2/3", "L3/4", "L4", "L3/5/6", "L5/6", "L6", "L5/6", "L5/6", "L5/6", "L6", "L4/6", "L6", "L6", "L6", "L6b"))
#markers_lake =



obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% markers_lake$gene,]
#mat_a = t(apply(obj_int_filt_avg_filt, 1, scale))
mat_a = t(apply(obj_int_filt_avg_filt, 1, scale_01))
#mat_a = obj_int_filt_avg_filt
colnames(mat_a) = colnames(obj_int_filt_avg_filt)
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(mat_a), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(mat_a)/palette_length, max(mat_a), length.out=floor(palette_length/2)))

ann_col = markers_lake %>% as.data.frame %>% column_to_rownames(var="gene")# %>% select(-gene)
pheatmap(mat_a, 
         annotation_row = ann_col, 
         show_colnames=T, 
         show_rownames=T,
         # cluster_cols=hc_r, 
         # cluster_rows = hc_c, 
         treeheight_row = 8,
         treeheight_col = 8, 
         #color=colors, 
         border_color = NA,
         # breaks=breaks,
         # display_numbers = obj_cor_sig_adj_label,
         filename=file.path(out_dir, "human_layer_marker_genes.pdf"), height=5, width=8)



# Using TF -------------------------------------------------------------------


# tf_genes = get_tf_genes_human()
# 
# 
# 
# 
# mm_mat_filt = tas_obj_avg_sct[rownames(tas_obj_avg_sct) %in% rownames(obj_int_filt_avg),]
# mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
# mm_mat_filt = mm_mat_filt[,!grepl("CR", colnames(mm_mat_filt))]
# mm_mat_filt = mm_mat_filt[,grep("ALM", colnames(mm_mat_filt))]
# 
# obj_int_filt_avg = as.matrix(obj_int_filt_avg)
# obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]
# obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("Pre", colnames(obj_int_filt_avg_filt))]
# obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("IEG", colnames(obj_int_filt_avg_filt))]
# #obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!(colnames(obj_int_filt_avg_filt) %in% c("Exc_NEUROD6"))]
# ## Filter for markers or variable genes 
# obj_int_filt_markers_to_use = obj_int_filt_markers %>%
#   #filter(!(gene %in% tx_mt$gene_name)) %>%
#   filter(gene %in% tf_genes$external_gene_name) %>% 
#   filter(p_val_adj < .05)
# #group_by(cluster) %>%
# #top_n(-100, p_val_adj) %>%
# #top_n(100, avg_logFC)
# obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% obj_int_filt_markers_to_use$gene,]
# 
# mm_md_n = tas_samp %>% filter(grepl("Glutamatergic", class))
# mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$cluster]
# mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]
# #obj_int_filt_avg_filt_log = log(obj_int_filt_avg_filt+1)
# #obj_int_filt_avg_filt_log = t(apply(obj_int_filt_avg_filt, 1, function(x) x - min(x) + 1))
# #mm_mat_filt_log = log(mm_mat_filt+1)
# #obj_int_filt_avg_filt_log = t(apply(obj_int_filt_avg_filt, 1, function(x) x - min(x)))
# #mm_mat_filt_log = t(apply(mm_mat_filt, 1, function(x) x - min(x)))
# obj_int_filt_avg_filt_log = obj_int_filt_avg_filt
# mm_mat_filt_log = mm_mat_filt
# #obj_int_filt_avg_filt_log_spec = t(apply(obj_int_filt_avg_filt_log, 1, calc_gene_specificity, scale=F))
# #mm_mat_filt_log_spec = na.omit(t(apply(mm_mat_filt_log, 1, calc_gene_specificity, scale=T)))
# obj_int_filt_avg_filt_log_spec = t(apply(obj_int_filt_avg_filt_log, 1, function(x) log(x / mean(x))))
# mm_mat_filt_log_spec = na.omit(t(apply(mm_mat_filt_log, 1, function(x) log(x / mean(x)))))
# obj_int_filt_avg_filt_log_spec = obj_int_filt_avg_filt_log_spec[match(rownames(mm_mat_filt_log_spec), rownames(obj_int_filt_avg_filt_log_spec)),]
# 
# obj_cor = cor(obj_int_filt_avg_filt_log_spec, mm_mat_filt_log_spec, method="spearman")
# 
# mat_a = obj_int_filt_avg_filt_log_spec
# mat_b = mm_mat_filt_log_spec
# obj_cor_sig = matrix(ncol = ncol(obj_cor), nrow=nrow(obj_cor), dimnames = list(colnames(mat_a), colnames(mat_b)))
# 
# for (i in 1:ncol(mat_a)) {
#   for (j in 1:ncol(mat_b)) {
#     obj_cor_sig[i,j] = cor.test(obj_int_filt_avg_filt_log_spec[,i], mm_mat_filt_log_spec[,j], method="spearman", exact=F)$p.value
#   }
# }
# 
# obj_cor_sig_adj = matrix(p.adjust(c(obj_cor_sig), method="BH"), ncol=ncol(obj_cor), nrow=nrow(obj_cor), dimnames = list(colnames(mat_a), colnames(mat_b)))
# 
# sig_thresh = .05
# obj_cor_sig_adj_df = melt(obj_cor_sig_adj) %>% filter(value < sig_thresh)
# obj_cor_filt = obj_cor[,colnames(obj_cor) %in% obj_cor_sig_adj_df$Var2]
# 
# obj_cor_sig_adj_label = matrix("", nrow=nrow(obj_cor_sig_adj), ncol=ncol(obj_cor_sig_adj), dimnames = dimnames(obj_cor))
# obj_cor_sig_adj_label[obj_cor_sig_adj<sig_thresh] = "*"
# obj_cor_sig_adj_label = obj_cor_sig_adj_label[, colnames(obj_cor_sig_adj_label) %in% colnames(obj_cor_filt)]
# 
# 
# ann_col = mm_md_n %>% filter(cluster %in% colnames(mm_mat_filt)) %>% distinct(cluster, .keep_all=T) %>% select(cluster, subclass) %>% as.data.frame()
# rownames(ann_col) = ann_col$cluster
# ann_col = ann_col %>% select(-cluster)
# 
# obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
# #obj_cor2 = obj_cor1
# #obj_cor2[obj_cor_sig_adj>.05] = 0
# palette_length = 50
# colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)
# 
# breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
#             seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))
# 
# hc_r = hclust(dist(obj_cor1), method="ward.D")
# hc_c = hclust(dist(t(obj_cor1)), method="ward.D")
# pheatmap(obj_cor1, annotation_row = ann_col,show_colnames=T, cluster_cols=T, cluster_rows = T, annotation_names_col = T, color=colors, breaks = breaks)
# pheatmap(obj_cor1, annotation_row = ann_col, show_colnames=T,  color=colors, filename=file.path(out_dir, "mouse_layers_bf_neurons_exc_tfs_specificity_heatmap_horiz.pdf"))
# 
# 
# pheatmap(t(obj_cor1), annotation_col = ann_col,show_colnames=T, cluster_cols=hc_r, cluster_rows = hc_c, annotation_names_col = T, color=colors, breaks = breaks, border_color = NA)
# pheatmap(t(obj_cor1), 
#          annotation_col = ann_col, 
#          show_colnames=T, 
#          cluster_cols=hc_r, 
#          cluster_rows = hc_c, 
#          treeheight_row = 8,
#          treeheight_col = 8, 
#          color=colors, 
#          border_color = NA,
#          breaks=breaks,
#          display_numbers = obj_cor_sig_adj_label,
#          filename=file.path(out_dir, "mouse_layers_bf_neurons_exc_tfs_specificity_heatmap.pdf"), height=5, width=8)
# 
# 
