
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
library(furrr)


library(cowplot)
theme_set(theme_cowplot())

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "anchor5_filter200_score_30")
script_name = "comparative_zeisel_cns_neurons_correlation_glut_by_ident"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")

## Markers
assay_to_use = "SCT"
markers_fname = file.path(dev_out_sub_dir, 
                          sprintf("marker_genes_position2_cluster_int_sub2_glut_%s.rds", assay_to_use))

## Comparison 
loom_data_dir = c("~/data/scrna_datasets/zeisel_2018/")
compare_prefix = "l5_Glutamate_seurat_telencephalon_projecting_excitatory"

# Load comparison ---------------------------------------------------------

compare_sub_fname = file.path(loom_data_dir, 
                              sprintf("%s.rds", compare_prefix))
redo = F
if (redo) {
  compare_obj1 = subset(compare_obj, subset=Taxonomy_group == "Telencephalon projecting excitatory neurons")
  compare_obj1 = SCTransform(compare_obj1, 
                             return.only.var.genes = F, 
                             vars.to.regress = c("MitoRiboRatio")
  )
  qsave(compare_obj1, compare_sub_fname)
} else {
  compare_obj1 = qread(compare_sub_fname)
}

compare_obj1$idents = Idents(compare_obj1)
md = FetchData(compare_obj1, c("idents", "Description", "Probable_location")) %>%
  distinct(idents, .keep_all=T)

# Average -----------------------------------------------------------------
avg_list_fname = file.path(loom_data_dir, sprintf("%s_ident_sct_avg.rds", compare_prefix))

redo_avg = F
if (redo_avg) {
  compare_avg = AverageExpression(compare_obj1, assays = "SCT")
  qsave(compare_avg, avg_list_fname)
} else {
  compare_avg = qread(avg_list_fname)
}

compare_avg1 = compare_avg[["SCT"]]

# Load finch data ---------------------------------------------------------

obj_int_filt = qread(data_out_obj_fname)
obj_int_filt_markers = readRDS(markers_fname)
res_to_use = "cluster_int_sub2"
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]


# Finch marker -----------------------------------------------------------


redo_markers = F
if (redo_markers) {
  Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
  DefaultAssay(obj_int_filt) = assay_to_use
  plan(multiprocess(workers = 2))
  obj_int_filt_markers = FindAllMarkers(obj_int_filt, 
                               test.use = "wilcox", 
                               max.cells.per.ident = 200,
                               only.pos = F)
  saveRDS(obj_int_filt_markers, markers_fname)
} else {
  obj_int_filter_markers = readRDS(markers_fname)
}

# Compare marker ID --------------------------------------------------------------

compare_markers_fname = file.path(loom_data_dir, sprintf("%s_markers.rds", compare_prefix))
redo_compare_markers = F
if (redo_compare_markers) {
  genes_to_test = intersect(rownames(obj_int_filt), rownames(compare_obj1))

  compare_markers = FindAllMarkers(compare_obj1, 
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

mm_mat_filt = compare_avg1[rownames(compare_avg1) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
mm_mat_filt = mm_mat_filt[,!(colnames(mm_mat_filt) %in% c("TEGLU11"))]

obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]
obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("Pre", colnames(obj_int_filt_avg_filt))]

## Filter for markers or variable genes 
ngenes = 200
compare_markers_sig = compare_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  filter(p_val_adj < .05) %>% 
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC) %>%
  filter(gene %in% compare_markers_sig$gene)

obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% obj_int_filt_markers_to_use$gene,]

mm_md_n = compare_obj1@meta.data
mm_md_n$idents = Idents(compare_obj1)

mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$ident]
mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]

mat_a = log1p(obj_int_filt_avg_filt) + .1
mat_b = log1p(mm_mat_filt) + .1
obj_cor = specificity_correlate(mat_a, mat_b, method="spearman")

# Shuffle -----------------------------------------------------------------


obj_cor_shuf_fname = file.path(out_dir, "obj_cor_shuf.qs")
redo = F
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

ann_col = mm_md_n %>%
  filter(idents %in% colnames(mm_mat_filt)) %>%
  distinct(idents, .keep_all=T) %>% 
  select(idents, Taxonomy_group, Region, Probable_location) %>%

  mutate(pl2 = case_when(grepl("Piriform", Probable_location)~"ventral pallium",
                         grepl("olfactory", Probable_location)~"ventral pallium",
                         grepl("deep", Probable_location) ~ "cortex / deep layer",
                         grepl("superficial", Probable_location)~"cortex / upper layer",
                      
                         grepl("amygdala", Probable_location)~"BLA",
                         grepl("Hippocampus", Region)~"HPC",
                         
                         grepl("layer 6", Probable_location)~"cortex / deep layer",
                         grepl("layer 5", Probable_location)~"cortex / deep layer",
                         grepl("Subiculum", Probable_location)~"HPC",
                         grepl("layer 4", Probable_location)~"cortex / layer 4",
                         grepl("layer 2", Probable_location)~"cortex / upper layer")
        
       ) %>%
  select(idents, pl2) %>% 
  mutate(pl2 = factor(pl2, levels=c("cortex / upper layer", "cortex / layer 4", "ventral pallium", "cortex / deep layer", "BLA", "HPC"))) %>%
  arrange(pl2) %>%
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

hc_r = hclust(dist(obj_cor1), method="average")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

pheatmap(obj_cor1, 
         cluster_rows = hc_r,
         cluster_cols = hc_c,
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
         height=4, width=3)

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
         height=4, width=4)



ct_order = rownames(ann_col)

ct_order2 = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5",  "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
obj_cor1 = obj_cor1[intersect(ct_order, rownames(obj_cor1)),ct_order2]
obj_cor_sig_adj_label_order = obj_cor_sig_adj_label[ct_order2,intersect(ct_order, colnames(obj_cor_sig_adj_label))]

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0 

colors =  colorRampPalette(brewer_pal(palette="Blues")(9))(palette_length)
breaks =  seq(0, max(obj_cor1_floor), length.out=palette_length+1)

pl2 = unique(ann_col$pl2)
pl2_colors = brewer_pal(type="qual")(length(pl2))
names(pl2_colors) = pl2
pheatmap(obj_cor1_floor, 
         annotation_row = ann_col, 
         annotation_colors = list(pl2 = pl2_colors),
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
         display_numbers = t(obj_cor_sig_adj_label_order),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_clust2.pdf"), height=3, width=4.5)


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
                species = c("mouse", "finch"))
md = md %>%
  mutate(idents = name) %>% 
  left_join(mm_md_n %>% distinct(idents, Description, Region))
write.table(md, file=file.path(out_dir, "nodes.txt"), sep="\t", row.names=F, quote=F)