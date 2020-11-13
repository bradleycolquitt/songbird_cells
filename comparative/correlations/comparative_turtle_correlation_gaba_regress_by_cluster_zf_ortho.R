
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

library(cowplot)
theme_set(theme_cowplot())
library(ComplexHeatmap)
library(scales)
library(circlize)
# Directories -------------------------------------------------------------

## Data dir
dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_gaba_int_sub2_regress_zf_ortho"
tree_dir = file.path(tree_dir, script_data_name)
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
out_dir = file.path(dev_out_dir, "comparative", "correlations")
script_name = "comparative_turtle_correlation_gaba_regress_by_cluster_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Markers
assay_to_use = "SCT"
markers_fname = file.path(tree_dir, sprintf("marker_genes_position2_cluster_int_sub2_gaba_%s.rds", assay_to_use))

## Compare dir
compare_data_dir = file.path("~/data/scrna_datasets/tosches_2018/")

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
fp = "turtle.neurons.gaba.v3.SCT"
compare_process_fname = file.path(compare_data_dir, sprintf("%s.qs", fp))
redo = F
if (redo) {
  compare_obj = readRDS(compare_fname)
  cells = Cells(compare_obj)[grepl("i", compare_obj$clusters)]
  compare_obj = subset(compare_obj, cells=cells)
  plan(multiprocess(workers=6))
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
compare_md = compare_obj@meta.data %>% distinct(clusters, .keep_all=T)

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


# Comparison markers --------------------------------------------------------------


compare_markers_fname = file.path(compare_data_dir, sprintf("%s_markers.rds", fp))
redo_compare_markers = F
if (redo_compare_markers) {
  DefaultAssay(compare_obj) = "SCT"
  Idents(compare_obj) = compare_obj$clusters
  genes_to_test = rownames(compare_obj)
  genes_to_test = genes_to_test[genes_to_test %in% rownames(compare_obj_avg_sct)]
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
obj_cor_filt = obj_cor[,colnames(obj_cor) %in% obj_cor_df_filt$compare]

obj_cor_df = obj_cor_df %>% mutate(sig_label = if_else(value>value_q95, "*", ""))
obj_cor_sig_adj_label = acast(obj_cor_df, celltype~compare, value.var="sig_label")


# Plot heatmap ------------------------------------------------------------

ann_col = mm_md_n %>% filter(clusters %in% colnames(obj_cor_filt)) %>%
  distinct(clusters, .keep_all=F) %>% 
  left_join(class_df) %>% 
  column_to_rownames(var="clusters")

classes = unique(class_df$class)
class_colors = qualpalr::qualpal(length(classes))$hex
names(class_colors) = classes

classes2 = unique(class_df$class2)
class2_colors = pal_aaas()(length(classes2))
names(class2_colors) = classes2

obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label)

obj_cor_sig_adj_label1 = obj_cor_sig_adj_label1[rownames(obj_cor1), colnames(obj_cor1)]
palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

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

ha = rowAnnotation(df = ann_col, 
                   annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                   col = list(class = class_colors, class2 = class2_colors),
                   annotation_width=unit(target_width*.2, "in"),
                   simple_anno_size_adjust=T,
                   annotation_legend_param = list(
                     class = list(gp = gpar(fontsize=6)),
                     class2 = list(gp = gpar(fontsize=6)))
)

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

# Network -----------------------------------------------------------------

md = data.frame(name = c(colnames(obj_cor1), rownames(obj_cor1)), 
                species = c(rep("finch", length(colnames(obj_cor1))), 
                            rep("turtle", length(rownames(obj_cor1))))
)

md = md %>%
  mutate(clusters = name) %>% 
  left_join(class_df) %>%

  mutate(source_label = if_else(species!="finch", paste( class, clusters, sep=" - "), clusters)) %>%
  select(-name) %>% 
  dplyr::rename(Region = class,
                source = clusters,
                name = source_label) 
n_top_cor = 4

cor_df = obj_cor_df
colnames(cor_df)[1:2] = c("target", "source")
cor_df = cor_df %>% group_by(target) %>% 
  top_n(n_top_cor, value) %>%
  filter(value >= .15) %>% 
  filter(value > value_q95) %>% 
  left_join(md %>% filter(species=="turtle")) %>%
  mutate(source = if_else(species=="turtle", name, source)) %>%
  select(target, source, value)


write.table(cor_df, file=file.path(out_dir, sprintf("cors_top%s.txt", n_top_cor)), sep="\t", row.names=F, quote=F)


write.table(md, file=file.path(out_dir, "nodes.txt"), sep="\t", row.names=F, quote=F)
