
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


library(cowplot)
theme_set(theme_cowplot())

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "anchor5_filter200_score_30")
script_name = "comparative_turtle_correlation_gaba_regress_by_cluster"
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

## Compare dir
compare_data_dir = file.path("~/data/scrna_datasets/tosches_2018/")


# Load comparison data ----------------------------------------------------

compare_fname = file.path(compare_data_dir, "turtle.neurons.v3.rds")
print(compare_fname)
fp = "turtle.neurons.gaba.v3.SCT"
compare_process_fname = file.path(compare_data_dir, sprintf("%s.qs", fp))
redo = F
if (redo) {
  compare_obj = readRDS(compare_fname)
  cells = Cells(compare_obj)[grepl("i", compare_obj$clusters)]
  compare_obj = subset(compare_obj, cells=cells)
  compare_obj = SCTransform(compare_obj,
                            return.only.var.genes = F,
                            vars.to.regress = c("percent.mito")
  )
  qsave(compare_obj, compare_process_fname)
} else {
  compare_obj = qread(compare_process_fname)
}

compare_obj_avg = AverageExpression(compare_obj, assays="SCT")
compare_obj_avg_sct = compare_obj_avg[["SCT"]]
compare_md = compare_obj@meta.data


# Load finch data ---------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)
obj_int_filt_markers = readRDS(markers_fname)

res_to_use = "cluster_int_sub2"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]



# By Cluster --------------------------------------------------------------

compare_markers_fname = file.path(compare_data_dir, sprintf("%s_markers.rds", fp))
redo_compare_markers = F
if (redo_compare_markers) {
  DefaultAssay(compare_obj) = "SCT"
  Idents(compare_obj) = compare_obj$clusters
  genes_to_test = rownames(compare_obj)
  genes_to_test = genes_to_test[genes_to_test %in% rownames(compare_obj)]
  compare_markers = FindAllMarkers(compare_obj, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, only.pos=T)
  saveRDS(compare_markers, compare_markers_fname)
} else {
  compare_markers = readRDS(compare_markers_fname)
}


class_df = data.frame(clusters = c(paste0("i0",  1:9), paste0("i", 10:18)),
                      class = c("OB", "septum", "septum", "striatum", "amygdala", "amygdala", 
                                "SST", "SST", "SST", 'SST',
                                "PV-like", "PV-like", "PV-like",
                                "HTR3A VIP-like", "HTR3A VIP-like", "HTR3A VIP-like",
                                "HTR3A Reln", "HTR3A Reln"),
                      class2 = c("LGE-derived", "septum", "septum", "LGE-derived", "LGE-derived", "LGE-derived",
                                 "MGE-derived", "MGE-derived", "MGE-derived", "MGE-derived",
                                 "MGE-derived", "MGE-derived", "MGE-derived",
                                 "CGE-derived", "CGE-derived", "CGE-derived",
                                 "CGE-derived", "CGE-derived"))

# Calc specificity --------------------------------------------------------

mm_mat_filt = compare_obj_avg_sct[rownames(compare_obj_avg_sct) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
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

mm_md_n = compare_obj@meta.data
mm_md_n$idents = Idents(compare_obj)

mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$ident]
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

ann_col = mm_md_n %>% filter(clusters %in% colnames(mm_mat_filt)) %>%
  distinct(clusters, .keep_all=F) %>% 
  column_to_rownames(var="clusters")

classes = unique(class_df$class2)
class_colors = pal_locuszoom()(length(classes))
names(class_colors) = classes
obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")
pheatmap(obj_cor1,
         annotation_row = ann_col,show_colnames=T, 
         cluster_cols=hc_c, 
         cluster_rows = hc_r,
         annotation_names_col = T,
         treeheight_col = 0,
         treeheight_row = 8,
         color=colors, 
         breaks = breaks)

pheatmap(obj_cor1, annotation_row = ann_col, 
         show_colnames=T,  
         color=colors, 
         breaks=breaks,
         cluster_cols=hc_c, 
         cluster_rows = hc_r,
         treeheight_col = 0,
         treeheight_row = 8,
         border_color = NA,
         display_numbers = t(obj_cor_sig_adj_label),
         filename=file.path(out_dir, "cluster_markers_specificity_heatmap.pdf"), height=6, width=5)

pheatmap(t(obj_cor1), annotation_col = ann_col, show_rownames=T,  color=colors, 
         cluster_rows = hc_c, 
         cluster_cols = hc_r,
         breaks=breaks,
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
                species = c(rep("finch", length(colnames(obj_cor1))), 
                            rep("turtle", length(rownames(obj_cor1))))
)
md = md %>%
  mutate(clusters = name) %>% 
  left_join(mm_md_n %>% distinct(clusters, class, class2)) %>%
  rename(Region = class)
write.table(md, file=file.path(out_dir, "nodes.txt"), sep="\t", row.names=F, quote=F)
