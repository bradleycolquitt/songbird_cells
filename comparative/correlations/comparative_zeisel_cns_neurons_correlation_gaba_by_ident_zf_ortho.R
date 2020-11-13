
source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
library(viridis)
library(pheatmap)
#library(here)
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
library(furrr)
library(future)
library(ComplexHeatmap)
library(scales)
library(circlize)

library(cowplot)
theme_set(theme_cowplot())
#here = here::here

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
script_name = "comparative_zeisel_cns_neurons_correlation_gaba_by_ident_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Markers
assay_to_use = "SCT"
markers_fname = file.path(tree_dir, sprintf("marker_genes_position2_cluster_int_sub2_gaba_%s.rds", assay_to_use))

## Comparison 
loom_data_dir = c("~/data/scrna_datasets/zeisel_2018/")
compare_prefix = "l6_r2_cns_neurons_seurat"

# Load orthologs ----------------------------------------------------------


to_use = c("zf", "mm")


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

compare_prefix1 = "l6_r2_cns_neurons_seurat_GABA_Acetyl"
compare_sub_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix1))

redo = F
if (redo) {
  obj_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix))
  zei_obj_in = qread(obj_fname)
  zei_obj_in$idents = Idents(zei_obj_in)
  nt_select = "GABA|Acetyl"
  cells = Cells(zei_obj_in)[grepl(nt_select, zei_obj_in$Neurotransmitter)]
  zei_obj_in1 = subset(zei_obj_in, cells=cells)
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

egg_l_mm = egg_l %>% filter(species=="mm") %>% 
  distinct(predicted_name, gene)
zei_obj_in_avg1 = zei_obj_in_avg1[intersect(rownames(zei_obj_in_avg1), toupper(egg_l_mm$gene)),]
rownames(zei_obj_in_avg1) = egg_l_mm$predicted_name[match(rownames(zei_obj_in_avg1), toupper(egg_l_mm$gene))]

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

redo_finch_markers = F
if (redo_finch_markers) {
  
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  genes_to_test = intersect(rownames(tas_obj_avg_sc_sct), rownames(obj_int_filt))
  
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

# Compare marker ID --------------------------------------------------------------

compare_markers_fname = file.path(loom_data_dir, sprintf("%s_markers.rds", compare_prefix1))
redo_compare_markers = F

plan(multiprocess(workers = 6))
if (redo_compare_markers) {
  
  compare_markers = FindAllMarkers(zei_obj_in1, 
                                   assay="SCT", 
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

mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]
obj_int_filt_avg_filt_log = obj_int_filt_avg_filt
mm_mat_filt_log = mm_mat_filt

obj_int_filt_avg_filt_log_spec = t(apply(obj_int_filt_avg_filt_log, 1, function(x) log(x / mean(x))))
mm_mat_filt_log_spec = na.omit(t(apply(mm_mat_filt_log, 1, function(x) log(x / mean(x)))))
obj_int_filt_avg_filt_log_spec = obj_int_filt_avg_filt_log_spec[match(rownames(mm_mat_filt_log_spec), rownames(obj_int_filt_avg_filt_log_spec)),]

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

# Plot heatmap  ---------------------------------------------------------



mm_md_n = mm_md_n %>%
  mutate(Region = as.character(Region)) %>% 
  mutate(Region = case_when(Region == "Striatum dorsal, Striatum ventral" ~ "Septum",
                            Region == "Striatum dorsal" ~ "Dorsal striatum",
                            Region == "Striatum ventral" ~ "Ventral striatum",
                            TRUE ~ Region))

ann_col = mm_md_n %>%
  filter(idents %in% colnames(obj_cor_filt)) %>%
  distinct(idents, .keep_all=T) %>% 
  select(idents, Taxonomy_group, Region) %>%
  column_to_rownames(var="idents") %>% 
  as.data.frame()


regions = na.omit(unique(ann_col$Region))
region_colors = pal_aaas()(length(regions))
names(region_colors) = regions

taxs = unique(ann_col$Taxonomy_group)
tax_colors = qualpalr::qualpal(length(taxs))$hex
names(tax_colors) = taxs

obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label)

palette_length = 50
breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0


colors =  colorRampPalette(brewer_pal(palette="Blues")(9))(palette_length)

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

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
                   col = list(Region = region_colors,
                              Taxonomy_group = tax_colors),
                   annotation_width=unit(target_width*.2, "in"),
                   simple_anno_size_adjust=T,
                   annotation_legend_param = list(
                     Region = list(gp = gpar(fontsize=6)),
                     Taxonomy_group = list(gp = gpar(fontsize=6)))
)

hm = Heatmap(obj_cor1_floor, col=colors, name = "rho",
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
print(hm)
dev.off()
hm

