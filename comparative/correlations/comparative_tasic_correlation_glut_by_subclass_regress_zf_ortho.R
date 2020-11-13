
source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
library(viridis)
library(pheatmap)
library(Seurat)
library(qs)
library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)
library(rtracklayer)
library(ggsci)
library(proxy)
library(scales)
library(future)
library(furrr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
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
script_name = "comparative_tasic_correlation_glut_by_subclass_regress_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Markers
assay_to_use = "SCT"
markers_fname = file.path(out_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))

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
fname_obj = file.path(tas_data_dir, sprintf("%s_seurat_glut_ortho_zf.qs", fp_name))
fname_avg_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_list.rds", fp_name))
fname_avg_sc_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_subclass_list.rds", fp_name))
fname_avg_ip_list = file.path(tas_data_dir, sprintf("%s_seurat_glut_data_avg_injection_target_list.rds", fp_name))


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
  tas_obj = subset(tas_obj, subset=class=="Glutamatergic")
  tas_obj = SCTransform(tas_obj, vars.to.regress = c( "percent_mt_exon_reads"), 
                        return.only.var.genes = T)
  
  qsave(tas_obj, fname_obj)
} else if (load_obj) {
  tas_obj = qread(fname_obj)
}


redo_avg = F
if (redo_avg) {
  Idents(tas_obj) = tas_obj@meta.data$cluster
  tas_obj_avg1 = AverageExpression(tas_obj, assays = c("SCT"), slot="data")

  
  saveRDS(tas_obj_avg1, fname_avg_list)
  Idents(tas_obj) = tas_obj@meta.data$subclass
  tas_obj_avg_sc = AverageExpression(tas_obj, assays = c("SCT"), slot="data")

  saveRDS(tas_obj_avg_sc, fname_avg_sc_list)
  
  ## Injection_primary
  Idents(tas_obj) = tas_obj@meta.data$injection_primary
  tas_obj_cur = subset(tas_obj, subset=injection_primary!="No Injection")
  
  tas_obj_avg_ip = AverageExpression(tas_obj_cur, assays = c("SCT"), slot="data")

  
  saveRDS(tas_obj_avg_ip, fname_avg_ip_list)
  
} else {
  tas_obj_avg1 = readRDS(fname_avg_list)
  tas_obj_avg_sc = readRDS(fname_avg_sc_list)
}
tas_obj_avg_sct = tas_obj_avg1[["SCT"]]
tas_obj_avg_sc_sct = tas_obj_avg_sc[["SCT"]]

tab_obj = FindVariableFeatures(tas_obj, assay = "SCT", nfeatures = 2000)


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


obj_int_filt = FindVariableFeatures(obj_int_filt, assay = "SCT", nfeatures = 2000)
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
                                        only.pos = T,
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

# Comparison markers ------------------------------------------------------

tas_exc_markers_fname = file.path(tas_data_dir, "tas_obj_glut_ortho_zf_subclass_markers.rds")
redo_tas_markers = F
if (redo_tas_markers) {
  
  Idents(tas_obj) = tas_obj$subclass
  genes_to_test = intersect(rownames(tas_obj), rownames(obj_int_filt))
  
  cells = Cells(tas_obj)[!(Idents(tas_obj)%in% c("NP", "CR"))]
  tas_obj_sub = subset(tas_obj, cells=cells)
  tas_exc_markers = FindAllMarkers(tas_obj_sub, 
                                   assay="SCT", 
                                   features = genes_to_test, 
                                   test.use = "t",
                                   min.pct = .2, 
                                   only.pos=T,
                                   max.cells.per.ident = 200
  )
  
  saveRDS(tas_exc_markers, tas_exc_markers_fname)
} else {
  tas_exc_markers = readRDS(tas_exc_markers_fname)
}


# Calc specificity --------------------------------------------------------

mm_mat_filt = tas_obj_avg_sc_sct[rownames(tas_obj_avg_sc_sct) %in% rownames(obj_int_filt_avg),]
mm_mat_filt = mm_mat_filt[!duplicated(rownames(mm_mat_filt)),]
mm_mat_filt = mm_mat_filt[,!grepl("CR|NP", colnames(mm_mat_filt))]

obj_int_filt_avg = as.matrix(obj_int_filt_avg)
obj_int_filt_avg_filt = obj_int_filt_avg[rownames(obj_int_filt_avg) %in% rownames(mm_mat_filt),]

## Filter for markers or variable genes 
ngenes = 400
tas_exc_markers_sig = tas_exc_markers %>% 
  filter(p_val_adj<.05) %>%
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC)

obj_int_filt_markers_to_use = obj_int_filt_markers %>%
  filter(cluster %in% colnames(obj_int_filt_avg_filt)) %>%
  filter(p_val_adj <.05) %>% 
  group_by(cluster) %>%
  top_n(ngenes, avg_logFC) %>%
  filter(gene %in% tas_exc_markers_sig$gene)
  #filter(gene %in% tas_exc_markers$gene)

#to_use = intersect(VariableFeatures(obj_int_filt,assay = "SCT" ), VariableFeatures(tas_obj, assay="SCT"))
to_use = obj_int_filt_markers_to_use$gene
obj_int_filt_avg_filt = obj_int_filt_avg_filt[rownames(obj_int_filt_avg_filt) %in% to_use,]

mm_md_n = tas_samp %>% filter(grepl("Glutamatergic", class))
mm_mat_filt = mm_mat_filt[,colnames(mm_mat_filt) %in% mm_md_n$subclass]
mm_mat_filt = mm_mat_filt[na.omit(match(rownames(obj_int_filt_avg_filt), rownames(mm_mat_filt))),]

colnames(mm_mat_filt) = if_else(colnames(mm_mat_filt)=="L4", "L4 IT", colnames(mm_mat_filt))

mat_a = log1p(obj_int_filt_avg_filt) + .1
mat_b = log1p(mm_mat_filt)  + .1
obj_cor = specificity_correlate(mat_a, mat_b, method="spearman")

# Shuffle -----------------------------------------------------------------


obj_cor_shuf_fname = file.path(out_dir, "obj_cor_shuf.qs")
redo = F
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


# Plot heatmap ------------------------------------------------------------

obj_cor1 = t(obj_cor)
obj_cor_sig_adj_label1 = t(obj_cor_sig_adj_label)
obj_cor_sig_adj_label1 = obj_cor_sig_adj_label1[rownames(obj_cor1),colnames(obj_cor1)]
palette_length = 50

hc_r = hclust(dist(obj_cor1), method="ward.D")
hc_c = hclust(dist(t(obj_cor1)), method="ward.D")

ct_order = c("L2/3 IT", "L4 IT","L5 IT", "L5 PT", "L6 IT", "L6 CT", "L6b")
obj_cor1 = obj_cor1[ct_order,]
obj_cor_sig_adj_label1 = obj_cor_sig_adj_label1[ct_order,]

ct_order2 = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5",  "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
obj_cor1 = obj_cor1[,ct_order2]
obj_cor_sig_adj_label1 = obj_cor_sig_adj_label1[,ct_order2]

obj_cor1_floor = obj_cor1
obj_cor1_floor[obj_cor1_floor<0] = 0 

colors =  colorRampPalette(brewer_pal(palette="Blues")(9))(palette_length)
breaks =  seq(0, max(obj_cor1_floor), length.out=palette_length+1)

star_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%s", obj_cor_sig_adj_label1[i, j]), x, y, gp = gpar(fontsize = 6))}

target_width = 1
target_height = target_width * nrow(obj_cor1) / ncol(obj_cor1) 

ht_opt(legend_labels_gp = gpar(fontsize=6),
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(obj_cor1), "in"),
       legend_grid_height = unit(target_height / nrow(obj_cor1), "in"))

hm = Heatmap(obj_cor1_floor, col=colors, name = "rho",
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             cluster_rows = F,
             cluster_columns=F,
             cell_fun = star_fun,
             row_names_gp = gpar(fontsize=6),
             column_names_gp = gpar(fontsize=6)
             
)

pdf(file.path(out_dir, "cluster_markers_specificity_heatmap_vert_no_clust2.pdf"), 
    height=target_height+4, 
    width=target_width+6)
print(hm)
dev.off()
