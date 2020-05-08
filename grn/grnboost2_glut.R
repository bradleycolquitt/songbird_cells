library(Seurat)
library(tidyverse)
library(qs)
library(reshape2)
library(future)
library(furrr)

library(ComplexHeatmap)
library(scales)
library(cowplot)
library(circlize)

source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")

# Directories -------------------------------------------------------------

dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "grn", "export_to_numpy_glut")
script_name = "grnboost2_glut"

cor_method = "spearman"
clust_method = "ward.D"
out_dir = file.path(out_dir, script_name, cor_method, clust_method)
dir.create(out_dir, recursive=T)

## GRNBoost dir
grn_dir = "~/data2/grnboost2/glut/"
grn_fname = file.path(grn_dir, "output.tsv")

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")

## Dim vis dir
red_dir = file.path(dev_dir, "reduction_viz")
script_red_name = "combined_int_sub2_glut_umap_regress"
red_dir = file.path(red_dir, script_red_name)
red_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")

# Load expr data ----------------------------------------------------------

tfs = get_tf_genes_human()
obj_int_filt = qread(data_out_obj_fname)
obj_int_filt = subset(obj_int_filt, subset=species=="zf")
obj_int_filt = FindVariableFeatures(obj_int_filt, nfeatures=3000)

# Load GRNBoost2 results --------------------------------------------------

grn = read_delim(grn_fname, delim="\t", col_names = c("source", "target", "importance"))
grn_mat = acast(grn, source~target, fill=0)
grn_cor = cor(t(grn_mat), method=cor_method)

# Filter --------------------------------------------------


importance_thresh = 8
grn_source_sum = grn %>% group_by(source) %>% summarize(importance_mean = mean(importance),
                                                        importance_median = median(importance),
                                                        n_thresh = sum(importance>importance_thresh))
plot(density(grn_source_sum$n_thresh))
grn_source_sum_filt = grn_source_sum %>% filter(n_thresh>10)

grn_source_filt = grn %>% inner_join(grn_source_sum_filt)

importance_thresh = 10
grn_source_target_sum = grn_source_filt %>% 
  ungroup() %>% 
  group_by(target) %>% 
  summarize(importance_mean = mean(importance),
            importance_median = median(importance),
            n_thresh = sum(importance>importance_thresh))

plot(density(grn_source_target_sum$n_thresh))
grn_source_target_sum_filt = grn_source_target_sum %>% 
  filter(n_thresh>=2) %>%
  select(target)

grn_source_target_filt = grn_source_filt %>% inner_join(grn_source_target_sum_filt)

plot(density(log(grn_source_target_filt$importance)))

grn_source_target_filt = grn_source_target_filt %>% filter(importance>1)
write.table(grn_source_target_filt, file.path(grn_dir, "output_source_target_filtered.tsv"), sep="\t", row.names=F)
grn_mat = acast(grn_source_target_filt, source~target, fill=0, value.var="importance")
grn_cor = cor(t(grn_mat), method=cor_method)

grn_cor_hc = hclust(dist(grn_cor), method="ward")
pdf(file.path(out_dir, "source_filt_tf_cor.pdf"), height=15, width=14)

hm = Heatmap(grn_cor, cluster_rows = grn_cor_hc, cluster_columns = grn_cor_hc)
hm
dev.off()

grn_cor_df = melt(grn_cor)
colnames(grn_cor_df) = c("source", "target", "value")
grn_cor_df = grn_cor_df %>% filter(value<1) %>%
  filter(value>.1)
write.table(grn_cor_df, file.path(grn_dir, "output_source_filtered_cor.tsv"), sep="\t", row.names=F, quote=F)

# Shuffle correlation -----------------------------------------------------

nreps = 1000
grn_cor_shuf = map(1:nreps, function(i) {
  tmp = t(apply(grn_mat, 1, function(x)  sample(x)))
  tmp_cor = cor(t(tmp), method=cor_method)
  m = melt(tmp_cor)
  colnames(m) = c("source", "target", "value")
  m
}) %>% bind_rows()


# Filter cor --------------------------------------------------------------

grn_cor_shuf1 = grn_cor_shuf %>%
  filter(value<.98) %>% 
  group_by(source, target) %>%
  summarize(q99 = quantile(value, .99))

grn_cor_df = grn_cor_df %>% left_join(grn_cor_shuf1)
grn_cor_df = grn_cor_df %>% mutate(value_filt = if_else(value>q99, value, 0)) 

grn_cor_df_filt = grn_cor_df %>% filter(value>q99)

grn_cor_filt = acast(grn_cor_df, source~target, value.var="value_filt", fill=0)
grn_cor_filt_rs = t(apply(grn_cor_filt, 1, function(x) sum(x>0) / length(x)))[1,]
to_keep = names(grn_cor_filt_rs)[grn_cor_filt_rs>.1]
grn_cor_filt  = grn_cor_filt[to_keep,to_keep]


# Plot TF expression ----------------------------------------------------------

res_to_use = "cluster_int_sub2"
grn_cor_hc = hclust(dist(grn_cor), method=clust_method)


to_use_tf = rownames(grn_cor_filt)
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
avg = AverageExpression(obj_int_filt, assays = "SCT", features=to_use_tf, slot = "data")
avg1 = log1p(avg[[1]])
avg1 = avg1[match(grn_cor_hc$labels, rownames(avg1)),]
avg_scale1 = t(apply(avg1, 1, function(x) (x - min(x)) / (max(x) - min(x))))




cols = colorRamp2(breaks=c(0,1), colors=brewer_pal(palette="Blues")(9)[c(1,9)])

hc = hclust(dist(t(avg_scale1)), method="ward.D")
hm = Heatmap(avg_scale1,
             col = cols,
             cluster_columns = hc,
             column_names_side="top",
             cluster_rows = grn_cor_hc,
             show_row_dend=F,
             heatmap_width=unit(4, "in"),
             heatmap_legend_param = list(title="Scaled expression", direction="horizontal")
             )
hm
pdf(file.path(out_dir, "tf_expr_heatmap.pdf"), height=length(to_use_tf) / 4, width=8)
draw(hm, heatmap_legend_side="bottom")
dev.off()

grn_cor_filt_hc = hclust(dist(grn_cor_filt), method=clust_method)
to_use_tf = rownames(grn_cor_filt)

avg_filt = AverageExpression(obj_int_filt, assays = "SCT", features=to_use_tf, slot = "data")
avg_filt1 = log1p(avg_filt[[1]])
avg_filt1 = avg_filt1[match(grn_cor_filt_hc$labels, rownames(avg_filt1)),]

avg_filt_scale1 = t(apply(avg_filt1, 1, function(x) (x - min(x)) / (max(x) - min(x))))


cols = colorRamp2(breaks=c(0,1), colors=brewer_pal(palette="Blues")(9)[c(1,9)])

hc_filt = hclust(dist(t(avg_filt_scale1)), method="ward.D")
hm = Heatmap(avg_filt_scale1,
             col = cols,
             cluster_columns = hc_filt,
             column_names_side="top",
             cluster_rows = grn_cor_filt_hc,
             show_row_dend=F,
             heatmap_width=unit(4, "in"),
             heatmap_legend_param = list(title="Scaled expression", direction="horizontal")
)

pdf(file.path(out_dir, "tf_expr_filter_heatmap.pdf"), height=length(to_use_tf) / 4, width=8)
draw(hm, heatmap_legend_side="bottom")
dev.off()

# Plot correlation heatmap ----------------------------------------------------------

pdf(file.path(out_dir, "source_filt_tf_cor2.pdf"), height=20, width=20)
grn_cor[grn_cor>.98] = 0
cols = colorRamp2(c(0,max(c(grn_cor))), colors = brewer_pal(palette="Reds")(9)[c(1,9)])

hm = Heatmap(grn_cor, 
             col= cols, heatmap_height = unit(18,"in"), heatmap_width = unit(18, "in"),
             show_row_dend = T,
             show_column_dend = F,
             cluster_rows = grn_cor_hc, 
             cluster_columns = grn_cor_hc
)
hm
dev.off()


pdf(file.path(out_dir, "source_filt_tf_cor2_filt.pdf"), height=20, width=20)
grn_cor_filt[grn_cor_filt>.98] = 0
cols = colorRamp2(c(0,max(c(grn_cor_filt))), colors = brewer_pal(palette="Reds")(9)[c(1,9)])

hm = Heatmap(grn_cor_filt, 
             col= cols, heatmap_height = unit(18,"in"), heatmap_width = unit(18, "in"),
             show_row_dend = T,
             show_column_dend = F,
             cluster_rows = grn_cor_filt_hc, 
             cluster_columns = grn_cor_filt_hc
)
hm
dev.off()