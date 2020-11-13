
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
library(furrr)
library(future)
library(ComplexHeatmap)
library(scales)
library(circlize)

library(cowplot)
theme_set(theme_cowplot())


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
script_name = "comparative_tasic_correlation_gaba_regress_by_cluster_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Markers
assay_to_use = "SCT"
markers_fname = file.path(out_dir, sprintf("marker_genes_position2_cluster_int_sub2_gaba_%s.rds", assay_to_use))

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

egg_l_mm = egg_l %>% filter(species=="mm") %>% rename(gene_symbol = gene) %>% select(-species)
tas_genes = tas_genes %>% left_join(egg_l_mm)
tas_genes_filt = tas_genes %>% filter(!is.na(predicted_name)) %>%
  distinct(predicted_name, .keep_all=T)

tas_inject_md = read_xlsx(file.path(tas_data_dir, "NIHMS1001532-supplement-Supplementary_Table_6.xlsx"))
tas_samp$injection_category = tas_inject_md$injection_target_category[match(tas_samp$injection_primary, tas_inject_md$injection_target)]

fp_name = "ALM_VISp"
fname_obj = file.path(tas_data_dir, sprintf("%s_seurat_gaba_zf_ortho.qs", fp_name))
fname_avg_list = file.path(tas_data_dir, sprintf("%s_seurat_gaba_data_zf_ortho_avg_list.rds", fp_name))
fname_avg_sc_list = file.path(tas_data_dir, sprintf("%s_seurat_gaba_data_zf_ortho_avg_subclass_list.rds", fp_name))

redo_obj = F
load_obj = T
if (redo_obj) {
  tas_objs = map(fps, function(fp) {
    print(fp)
    tas_dat = fread(file.path(tas_data_dir, sprintf("%s_exon-matrix.csv", fp)))
    tas_dat = as.data.frame(tas_dat)
    rownames(tas_dat) = tas_dat$V1
    tas_dat_filt = tas_dat[rownames(tas_dat) %in% tas_genes_filt$gene_entrez_id,]
    rownames(tas_dat_filt) = tas_genes_filt$predicted_name[match(rownames(tas_dat_filt), tas_genes_filt$gene_entrez_id)]
    
    tas_dat_filt = tas_dat_filt[,-1]
    tas_dat_filt = as.matrix(tas_dat_filt)
    tas_obj = CreateSeuratObject(tas_dat_filt)
    tas_obj
  })
  
  tas_obj = merge(tas_objs[[1]], tas_objs[[2]], add.cell.ids =  fps_sample)
  tas_obj = AddMetaData(tas_obj, metadata = tas_samp)
  tas_obj = subset(tas_obj, subset=class=="GABAergic")
  options(future.globals.maxSize= 15000 * 1024^2)
  tas_obj = SCTransform(tas_obj, vars.to.regress = c( "percent_mt_exon_reads"), 
                        return.only.var.genes = T)
  
  qsave(tas_obj, fname_obj)
} else if (load_obj) {
  tas_obj = qread(fname_obj)
}


redo_avg = T
if (redo_avg) {
  Idents(tas_obj) = tas_obj@meta.data$cluster
  tas_obj_avg1 = AverageExpression(tas_obj, assays = c("SCT"), slot="data")
  saveRDS(tas_obj_avg1, fname_avg_list)
  
  Idents(tas_obj) = tas_obj@meta.data$subclass
  tas_obj_avg_sc = AverageExpression(tas_obj, assays = c("SCT"), slot="data")
  saveRDS(tas_obj_avg_sc, fname_avg_sc_list)

} else {
  tas_obj_avg1 = readRDS(fname_avg_list)
  tas_obj_avg_sc = readRDS(fname_avg_sc_list)
}
tas_obj_avg_sct = tas_obj_avg1[["SCT"]]
tas_obj_avg_sc_sct = tas_obj_avg_sc[["SCT"]]


# Load finch data ---------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)
#obj_int_filt_markers = readRDS(markers_fname)

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

# By Cluster --------------------------------------------------------------

tas_inh_markers_fname = file.path(tas_data_dir, "tas_obj_inh_markers.rds")
redo_tas_markers = F
if (redo_tas_markers) {
  tas_obj_inh = subset(tas_obj, subset=class=="GABAergic")
  DefaultAssay(tas_obj_inh) = "SCT"
  Idents(tas_obj_inh) = tas_obj_inh$cluster
  genes_to_test = intersect(rownames(tas_obj_avg_sc_sct), rownames(obj_int_filt_avg))
  tas_inh_markers = FindAllMarkers(tas_obj_inh, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, only.pos=T)
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
  plan(multiprocess(workers=4, gc=T))
  obj_cor_shuf = future_map(1:nrep, function(i) {
    mat_a_cur = mat_a
    mat_b_cur = mat_b
    rownames(mat_a_cur) = sample(rownames(mat_a_cur))
    obj_cor = specificity_correlate(mat_a_cur, mat_b_cur)
    
    
    obj_cor_df = melt(obj_cor)
    colnames(obj_cor_df) =  c("celltype", "compare", "value")
    obj_cor_df
  }, .progress=T) %>% bind_rows() 
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

# Plot heatmap ---------------------------------------------------------


ann_col = mm_md_n %>% filter(cluster %in% colnames(mm_mat_filt)) %>%
  distinct(cluster, .keep_all=T) %>%
  select(cluster, subclass) %>% 
  as.data.frame() %>%
  column_to_rownames("cluster")

classes = rev(c("Pvalb", "Sst", "Vip", "Sncg", "Serpinf1", "Lamp5", "Meis2"))
class_colors = pal_locuszoom()(length(classes))
names(class_colors) = classes

obj_cor1 = t(obj_cor_filt[,na.omit(match(rownames(ann_col), colnames(obj_cor_filt)))])
obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label)

palette_length = 50
colors =  colorRampPalette(c("navy", "white", "firebrick3"))(palette_length)

breaks =  c(seq(min(obj_cor1), 0, length.out=ceiling(palette_length/2) + 1), 
            seq(max(obj_cor1)/palette_length, max(obj_cor1), length.out=floor(palette_length/2)))

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0

palette_length = 50
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
                   col = list(subclass = class_colors),
                   annotation_width=unit(target_width*.2, "in"),
                   simple_anno_size_adjust=T,
                   annotation_legend_param = list(
                     subclass = list(gp = gpar(fontsize=6)))
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

# Network -----------------------------------------------------------------

n_top_cor = 4
cor_df = obj_cor_df
colnames(cor_df)[1:2] = c("target", "source")
cor_df = cor_df %>% group_by(target) %>% 
  top_n(n_top_cor, value) %>%
  filter(value >= .15) %>% 
  filter(value > value_q95)


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
