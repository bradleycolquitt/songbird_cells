
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
library(future)
library(furrr)
library(qualpalr)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(scales)
theme_set(theme_cowplot())

# Directories -------------------------------------------------------------


## Data dir
dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_in_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut_zf_ortho.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
out_dir = file.path(dev_out_dir, "comparative", "correlations")
script_name = "comparative_turtle_correlation_glut_regress_by_area_tfs_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)


## Markers
assay_to_use = "SCT"
markers_fname = file.path(tree_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))

## Compare dir
compare_data_dir = file.path("~/data/scrna_datasets/tosches_2018/")

tfs = get_tf_genes_human()

# Load orthologs ----------------------------------------------------------


to_use = c("zf", "tse")


egg = read.delim("~/brad/assembly/eggnog/eggnog_orthologs.txt", sep="\t")
egg1 = egg %>% 
  distinct(predicted_name, species, .keep_all=T) 

egg_w = egg1 %>% select(-protein_id) %>%  
  pivot_wider(names_from = species, values_from = gene, values_fill=list(gene=NA)) %>%
  mutate_at(vars(one_of("zf", "bf", "mm", "tse")), as.character) %>%
  select_at(vars(one_of("predicted_name", to_use)))

egg_exclude = apply(as.matrix(egg_w[,to_use]), 1, function(x)any(is.na(x)))
names(egg_exclude) = egg_w$predicted_name
egg_include = egg_exclude[!egg_exclude]
egg_w_f = egg_w %>% filter(predicted_name %in% names(egg_include))
egg_l = egg_w_f %>% pivot_longer(names_to = "species", cols = one_of(to_use), values_to="gene" ) %>%
  filter(!(predicted_name == "" | is.na(predicted_name)))


# Load comparison data ----------------------------------------------------

compare_fname = file.path(compare_data_dir, "turtle.neurons.v3.rds")
print(compare_fname)
fp = "turtle.neurons.glut.v3.SCT"
compare_process_fname = file.path(compare_data_dir, sprintf("%s.qs", fp))
redo = F
if (redo) {
  compare_obj = readRDS(compare_fname)
  cells = Cells(compare_obj)[grepl("e", compare_obj$clusters)]
  compare_obj = subset(compare_obj, cells=cells)
  compare_obj = SCTransform(compare_obj,
                            return.only.var.genes = F,
                            vars.to.regress = c("percent.mito")
  )
  qsave(compare_obj, compare_process_fname)
} else {
  compare_obj = qread(compare_process_fname)
}


Idents(compare_obj) = compare_obj$pallial.area
compare_obj_avg = AverageExpression(compare_obj, assays="SCT")
compare_obj_avg_sct = compare_obj_avg[["SCT"]]
compare_md = compare_obj@meta.data %>% distinct(pallial.area, .keep_all=T)

compare_obj_avg = AverageExpression(compare_obj, assays="SCT")
compare_obj_avg_sct = compare_obj_avg[["SCT"]]

egg_l_tse = egg_l %>% filter(species=="tse")
compare_obj_avg_sct = compare_obj_avg_sct[intersect(rownames(compare_obj_avg_sct), egg_l_tse$gene),]
rownames(compare_obj_avg_sct) = egg_l_tse$predicted_name[match(rownames(compare_obj_avg_sct), egg_l_tse$gene)]

# Load finch data ---------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)

res_to_use = "cluster_int_sub2"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
obj_int_filt_avg = obj_int_filt_avg[[1]]



egg_l_zf = egg_l %>% filter(species=="zf")
obj_int_filt_avg = obj_int_filt_avg[intersect(rownames(obj_int_filt_avg), egg_l_zf$gene),]
rownames(obj_int_filt_avg) = egg_l_zf$predicted_name[match(rownames(obj_int_filt_avg), egg_l_zf$gene)]

# Finch markers -----------------------------------------------------------

redo_finch_markers = T
if (redo_finch_markers) {
  
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  genes_to_test = intersect(rownames(compare_obj_avg_sct), rownames(obj_int_filt))
  
  obj_int_filt_markers = FindAllMarkers(obj_int_filt, 
                                        assay="SCT", 
                                        features = genes_to_test, 
                                        test.use = "t",
                                        min.pct = .2, 
                                        only.pos=T,
                                        max.cells.per.ident = 200
  )
  saveRDS(obj_int_filt_markers, markers_fname)
} else {
  obj_int_filt_markers = readRDS(markers_fname)
}

obj_int_filt_markers = obj_int_filt_markers %>% 
  mutate(gene = egg_l_zf$predicted_name[match(gene, egg_l_zf$gene)]) 

obj_int_filt_markers = obj_int_filt_markers %>%
  filter(!is.na(gene))


# Compare markers --------------------------------------------------------------

compare_ident_to_use = "pallial.area"
compare_markers_fname = file.path(compare_data_dir, sprintf("%s_%s_markers.rds", fp, compare_ident_to_use))
redo_compare_markers = T
if (redo_compare_markers) {
  DefaultAssay(compare_obj) = "SCT"
  Idents(compare_obj) = FetchData(compare_obj, compare_ident_to_use)
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


# Calc specificity --------------------------------------------------------

mm_mat_filt = compare_obj_avg_sct[rownames(compare_obj_avg_sct) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]
obj_int_filt_avg_filt = obj_int_filt_avg_filt[,!grepl("Pre", colnames(obj_int_filt_avg_filt))]

## Filter for markers or variable genes 
ngenes = 500
compare_markers_sig = compare_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  filter(p_val_adj < .05) %>% 
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC) %>%
#  filter(gene %in% compare_markers_sig$gene) %>%
  filter(gene %in% tfs$external_gene_name)

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
  plan(multiprocess(workers=4, gc=T))
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
obj_cor_filt = obj_cor
obj_cor_df = obj_cor_df %>% mutate(sig_label = if_else(value>value_q95, "*", ""))
obj_cor_sig_adj_label = acast(obj_cor_df, celltype~compare, value.var="sig_label")


# Plot heatmap ------------------------------------------------------------


ann_col = compare_md %>% filter(pallial.area %in% colnames(obj_cor_filt)) %>%
  distinct(pallial.area, .keep_all=F) 

rownames(ann_col) = ann_col$pallial.area
obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label[,intersect(rownames(obj_cor1), colnames(obj_cor_sig_adj_label))])


classes = unique(ann_col$pallial.area)
class_colors = qualpalr::qualpal(length(classes))$hex
names(class_colors) = classes

palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

ct_order2 = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5",  "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
obj_cor1 = obj_cor1[,ct_order2]
obj_cor_sig_adj_label1 = obj_cor_sig_adj_label1[,ct_order2]

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0 

colors =  colorRampPalette(brewer_pal(palette="Blues")(9))(palette_length)


star_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", obj_cor_sig_adj_label1[i, j]), x, y, gp = gpar(fontsize = 6))}

target_width = 1
target_height = target_width * nrow(obj_cor1) / ncol(obj_cor1) 

ht_opt(legend_labels_gp = gpar(fontsize=6),
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(obj_cor1), "in"),
       legend_grid_height = unit(target_height / nrow(obj_cor1), "in"))


region_col = list(pallial.area = class_colors)
ha = rowAnnotation(pallial.area = ann_col$pallial.area, col = region_col)
hm = Heatmap(obj_cor1_floor, col=colors, name = "rho",
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             cluster_rows = T,
             cluster_columns=F,
             cell_fun = star_fun,
             row_names_gp = gpar(fontsize=6),
             column_names_gp = gpar(fontsize=6),
             left_annotation = ha
)
hm

pdf(file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_clust2.pdf"), 
    height=target_height+4, 
    width=target_width+6)
print(hm)
dev.off()
