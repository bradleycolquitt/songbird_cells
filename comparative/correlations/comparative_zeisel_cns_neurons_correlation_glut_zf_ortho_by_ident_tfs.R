
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
library(ComplexHeatmap)
library(circlize)
library(scales)
#here = here::here

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
script_name = "comparative_zeisel_cns_neurons_correlation_glut_zf_ortho_by_ident_tfs"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Markers
assay_to_use = "SCT"
markers_fname = file.path(tree_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))

## Comparison 
loom_data_dir = c("~/data/scrna_datasets/zeisel_2018/")
compare_prefix = "l5_Glutamate_seurat_telencephalon_projecting_excitatory"

tfs = get_tf_genes_human()

# Load orthologs ----------------------------------------------------------


to_use = c("zf", "mm")


egg = read.delim("~/brad/assembly/eggnog/eggnog_orthologs.txt", sep="\t")
egg1 = egg %>% 
  mutate(predicted_name = toupper(predicted_name)) %>% 
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


egg_l_mm = egg_l %>% filter(species=="mm") %>%
  filter(!is.na(predicted_name)) %>%
  filter(predicted_name!="") %>%
  mutate(external_gene_name = toupper(gene))
tfs = tfs %>% left_join(egg_l_mm) %>%
  select(-external_gene_name) %>%
  dplyr::rename(external_gene_name = predicted_name)

# Load comparison ---------------------------------------------------------

compare_prefix1 = "l5_Glutamate_seurat_telencephalon_projecting_excitatory_zf_ortho"
compare_sub_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix1))

redo = F
if (redo) {
  na = "mm"
  egg_l_cur = egg_l %>% filter(species==na) %>%
    filter(!is.na(predicted_name)) %>%
    filter(predicted_name!="")
  
  obj_fname = file.path(loom_data_dir, sprintf("%s.qs", compare_prefix))
  compare_obj = qread(obj_fname)
  compare_obj$idents = Idents(compare_obj)

  dat = GetAssayData(compare_obj, assay="RNA", slot="counts")
  dat_new = dat[intersect(rownames(dat), toupper(egg_l_cur$gene)),]
  new_rnames = egg_l_cur$predicted_name[match(rownames(dat_new), toupper(egg_l_cur$gene))]
  dimnames(dat_new) = list(new_rnames, colnames(dat_new))
  assay_new =  CreateAssayObject(counts=dat_new)
  compare_obj[["ortho"]] = assay_new
  DefaultAssay(compare_obj) = "ortho"
  
  options(future.globals.maxSize=10000 * 1024^2)
  plan(multiprocess(workers=3))
  compare_obj = SCTransform(compare_obj, assay = "ortho",
                            return.only.var.genes = F,
                            vars.to.regress = c("MitoRiboRatio"))
  
  qsave(compare_obj, compare_sub_fname)
} else {
  compare_obj = qread(compare_sub_fname)
}

compare_obj$idents = Idents(compare_obj)
md = FetchData(compare_obj, c("idents", "Description", "Probable_location")) %>%
  distinct(idents, .keep_all=T)

compare_obj = FindVariableFeatures(compare_obj, nfeatures=6000)
compare_obj_var = VariableFeatures(compare_obj)

avg_list_fname = file.path(loom_data_dir, sprintf("%s_ident_sct_avg.rds", compare_prefix1))

redo_avg = T
if (redo_avg) {
  compare_avg = AverageExpression(compare_obj, assays = "SCT")
  qsave(compare_avg, avg_list_fname)
} else {
  compare_avg = qread(avg_list_fname)
}

compare_avg1 = compare_avg[["SCT"]]

# Load finch data ---------------------------------------------------------

redo = F
if (redo) {
  na = "zf"
  egg_l_cur = egg_l %>% filter(species==na) %>%
    filter(!is.na(predicted_name)) %>%
    filter(predicted_name!="")
  
  x = qread(data_in_obj_fname)
  x = subset(x, subset=species=="zf")
  dat = GetAssayData(x, assay="RNA", slot="counts")
  dat_new = dat[intersect(rownames(dat), egg_l_cur$gene),]
  new_rnames = egg_l_cur$predicted_name[match(rownames(dat_new), egg_l_cur$gene)]
  dimnames(dat_new) = list(new_rnames, colnames(dat_new))
  assay_new =  CreateAssayObject(counts=dat_new)
  x[["ortho"]] = assay_new
  DefaultAssay(x) = "ortho"
  
  options(future.globals.maxSize=10000 * 1024^2)
  plan(multiprocess(workers=3))
  obj_int_filt = SCTransform(x, assay="ortho",
                             return.only.var.genes = F,
                             vars.to.regress = c("percent.mito"))
  rm(x)
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

res_to_use = "cluster_int_sub2"
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt_avg1 = AverageExpression(obj_int_filt, assays = assay_to_use)
obj_int_filt_avg = obj_int_filt_avg1[[assay_to_use]]

egg_l_zf = egg_l %>% filter(species=="zf")
obj_int_filt_avg = obj_int_filt_avg[intersect(rownames(obj_int_filt_avg), egg_l_zf$gene),]
rownames(obj_int_filt_avg) = egg_l_zf$predicted_name[match(rownames(obj_int_filt_avg), egg_l_zf$gene)]

obj_int_filt = FindVariableFeatures(obj_int_filt, nfeatures=6000)
obj_int_filt_var = VariableFeatures(obj_int_filt)

# Calc specificity --------------------------------------------------------

mm_mat_filt = compare_avg1[rownames(compare_avg1) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]

obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]

## Filter for markers or variable genes 
to_use = Reduce(intersect, list(obj_int_filt_var, compare_obj_var, tfs$external_gene_name))
obj_int_filt_markers_to_use = data.frame(gene = to_use)
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

qsave(obj_cor_df, file.path(out_dir, "obj_cor_df.qs"))


# Heatmap ordered ---------------------------------------------------------
  
  pl_levels = c("Cortical pyramidal layer 2/3",
                "Cortical pyramidal layer 4",
                "Cortical pyramidal layer 5",
                "Cortex pyramidal layer 6",
                "Lateral cortex layer 6: gustatory, barrel field, auditory",
                "Cortex pyramidal layer 6b",
                "Cortical pyramidal (poor markers)",
                "Cingulate/Retrosplenial area, layer 2",
                "Cingulate/Retrosplenial area, layer 5",
                "Cingulate/Retrosplenial area, layer 6",
                "Cortex pyramidal layer 5, Cingulate/Retrosplenial area (superficial and deep)",
                "CA1",
                "CA1 (posterior),subiculum,entorhinal",
                "CA3",
                "Subiculum",
                "Entorhinal superficial layers",
                "Piriform pyramidal",
                "Anterior olfactory nucleus and ventral striatum",
                "Anterior olfactory nucleus, deep layer",
                "Basolateral amygdala"
  )
  
  ct_levels = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
  ann_col = mm_md_n %>%
    filter(idents %in% colnames(mm_mat_filt)) %>%
    distinct(idents, .keep_all=T) %>% 
    select(idents, Taxonomy_group, Region, Probable_location) %>%
    mutate(Probable_location = factor(Probable_location, levels=pl_levels)) %>% 
    arrange(Probable_location) %>%
    column_to_rownames(var="idents") %>% 
    as.data.frame()
  
  
  obj_cor1 = t(obj_cor[ct_levels,rownames(ann_col)])
  obj_cor_sig_adj_label = acast(obj_cor_df, celltype~compare, value.var="sig_label")
  obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label[ct_levels,rownames(ann_col)])

  obj_cor1_floor = obj_cor1
  obj_cor1_floor[obj_cor1_floor<0] = 0 

  colors = colorRamp2(breaks=c(0, max(obj_cor1_floor)), colors = brewer_pal(palette="Blues")(9)[c(1,9)])
  

  
  star_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%s", obj_cor_sig_adj_label1[i, j]), x, y, gp = gpar(fontsize = 6))}
  
  target_width = 1
  target_height = target_width * nrow(obj_cor1) / ncol(obj_cor1) 
  
  hm = Heatmap(obj_cor1_floor,
               col=colors,
               cluster_rows = F,
               cluster_columns = F,
               name = "rho",
               width = unit(target_width, "in"),
               height = unit(target_height, "in"),
               clustering_method_rows = "ward.D",
               clustering_method_columns = "ward.D",
               show_row_dend = F,
               show_column_dend = F,
               cell_fun = star_fun,
               row_names_gp = gpar(fontsize=6),
               column_names_gp = gpar(fontsize=6),
               #left_annotation = ha
               
               
  )
  hm
  pdf(file.path(out_dir, "heatmap_no_order.pdf"), 10, 10)
  print(hm)
  dev.off()
  