
source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
source("~/data2/rstudio/birds/utils/scRNA.R")
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


dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integration")
script_name = "comparative_tasic_seurat_glut_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name =  "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_in_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut_zf_ortho.qs")

# Params ------------------------------------------------------------------

params = list(k.anchor = 5,
              k.filter = 200,
              k.score = 30,
              max.features = 200,
              var.features = 2000)
dev_out_sub_dir= file.path(out_dir, 
                           sprintf("anchor%s_filter%s_score%s_maxfeatures%s_varfeatures%s", 
                                   params[["k.anchor"]], 
                                   params[["k.filter"]], 
                                   params[["k.score"]],
                                   params[["max.features"]],
                                   params[["var.features"]]))
dir.create(dev_out_sub_dir, recursive=T)

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

compare_obj = tas_obj

cells = Cells(compare_obj)[!(compare_obj$subclass %in% c("CR", "NP"))]
compare_obj = subset(compare_obj, cells=cells)
compare_obj$species = "mm"
compare_obj$cluster_int_sub2 = compare_obj$subclass

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

# Integrate ---------------------------------------------------------------

fname_anchors = file.path(dev_out_sub_dir, "obj_anchors.qs")
fname_int = file.path(dev_out_sub_dir, "obj_integrated_clustered.qs")

redo = T
redo_anchors = T

obj_list = list(zf = obj_int_filt, mm = compare_obj)
if (redo) {
  options(future.globals.maxSize= 15000 * 1024^2)
  
  obj_features <- SelectIntegrationFeatures(object.list = obj_list,
                                            assay=rep("SCT", times=length(obj_list), nfeatures = params[["var.features"]])
  )
  
  obj_list <- PrepSCTIntegration(object.list = obj_list, 
                                 assay = "SCT", 
                                 anchor.features = obj_features, 
                                 verbose = TRUE)
  
  if (redo_anchors) {
    obj_anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                          assay = rep("SCT", times=length(obj_list)),
                                          #reference = 1,
                                          normalization.method = "SCT",
                                          anchor.features = obj_features, 
                                          reduction = "cca",
                                          k.anchor = params[["k.anchor"]],
                                          k.filter = params[["k.filter"]],
                                          k.score = params[["k.score"]],
                                          max.features = params[["max.features"]],
                                          verbose = TRUE)
    qsave(obj_anchors, fname_anchors)
  } else {
    obj_anchors = qread(fname_anchors)
  }
 
  common_features = Reduce(intersect, map(obj_list, function(x) rownames(x@assays$SCT@data)))
  plan(multiprocess(workers=1))
  obj_int <- IntegrateData(anchorset = obj_anchors,
                           normalization.method = "SCT", 
                           features.to.integrate = common_features,
                           verbose = T)

  rm(obj_anchors)
  obj_int <- RunPCA(obj_int,  verbose = FALSE)
  obj_int <- RunUMAP(obj_int, dims = 1:30) 
  
  ress = seq(.2, 1, .2)
  obj_int = seurat_find_clusters(obj_int, ress, cluster_prefix = "cluster_int", dims.use=1:30, algorithm = 2)
  
  qsave(obj_int, fname_int)
} else {
  obj_int = qread(fname_int)
}

# Dimplot -----------------------------------------------------------------

redo_plot = T
if (redo_plot) {
  reduction = "umap"
  cats = c("orig.ident", "position", "species", "cluster_int_sub2")
  
  em = Embeddings(obj_int, "umap") %>% as.data.frame() %>% rownames_to_column()
  md = FetchData(obj_int, cats) %>% rownames_to_column()
  em = em %>% left_join(md)
  walk(cats, function(ca) {
    print(ca)
    em_ca = em %>% group_by(!!as.symbol(ca)) %>% 
      summarize(UMAP_1 = mean(UMAP_1),
                UMAP_2 = mean(UMAP_2))
    cas = unique(FetchData(obj_int, ca)[[1]])
    col = qualpalr::qualpal(length(cas))$hex
    names(col) = cas
    gg = ggplot(em, aes_string("UMAP_1", "UMAP_2", color=ca)) + 
      geom_point(alpha=.5, size=.5, shape=20) + 
      theme_void() + 
      scale_color_manual(values=col) + 
      guides(color = guide_legend(override.aes = list(size = 3) ) ) + 
      coord_equal()
    gg_leg = get_legend(gg)
    gg = gg + theme(legend.position="none")
    png(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s.png", ca)), width=4, height=4, units = "in", res = 300)
    print(gg)
    dev.off()
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_legend.pdf", ca)), gg_leg, base_width=4)
    
    gg = ggplot(em_ca, aes(UMAP_1, UMAP_2)) + 
      geom_text(aes_string(label=ca, color=ca)) + 
      theme_void() + 
      scale_color_manual(values=col)  + 
      coord_equal()
    gg = gg + theme(legend.position="none")
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_label.pdf", ca)), gg, base_width=4)
  })
  
  ca = "cluster_int_sub2"
  em_ca = em %>% group_by(!!as.symbol(ca), species) %>% 
    summarize(UMAP_1 = mean(UMAP_1),
              UMAP_2 = mean(UMAP_2))
  cas = unique(FetchData(obj_int, ca)[[1]])
  col = qualpalr::qualpal(length(cas))$hex
  names(col) = cas
  gg = ggplot(em, aes_string("UMAP_1", "UMAP_2", color=ca)) + 
    geom_point(alpha=.5, size=.5, shape=20) + 
    theme_void() + 
    scale_color_manual(values=col) + 
    guides(color = guide_legend(override.aes = list(size = 3) ) ) + 
    coord_equal() + 
    facet_grid(.~species)
  gg_leg = get_legend(gg)
  gg = gg + theme(legend.position="none")
  png(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_split-by-species.png", ca)), width=8, height=4, units = "in", res = 300)
  print(gg)
  dev.off()
  save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_split-by-species_legend.pdf", ca)), gg_leg, base_width=4, ncol=1)
  
  gg = ggplot(em_ca, aes(UMAP_1, UMAP_2)) + 
    geom_text(aes_string(label=ca, color=ca)) + 
    theme_void() + 
    scale_color_manual(values=col)  + 
    coord_equal() + 
    facet_grid(.~species)
  gg = gg + theme(legend.position="none")
  save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_split-by-species_label.pdf", ca)), gg, base_width=4, ncol=2)

}

# Transfer data -----------------------------------------------------------


obj_ref = obj_list[[2]]
  obj_query = obj_list[[1]]

transfer_anchors_fname = file.path(dev_out_sub_dir, "transfer_anchors.qs")
redo_transfer = F
if (redo_transfer) {
  
  
  obj_features <- SelectIntegrationFeatures(object.list = list(obj_ref, obj_query),
                                            assay=rep("SCT", times=2), nfeatures = params[["var.features"]]
  )
  transfer_anchors <- FindTransferAnchors(reference = obj_ref,
                                        query = obj_query, 
                                        features = obj_features,
                                        reference.assay="SCT",
                                        query.assay="SCT",
                                        reduction = "cca",
                                        k.anchor = params[["k.anchor"]],
                                        k.filter = params[["k.filter"]],
                                        k.score = params[["k.score"]],
                                        max.features = params[["max.features"]],
                                       
                                        dims = 1:30)
  qsave(transfer_anchors, transfer_anchors_fname)
  
  
} else {
  transfer_anchors = qread(transfer_anchors_fname)  
}


# Transfer labels -----------------------------------------------------------

ref_group = "cluster_int_sub2"
ref_group_sym = as.symbol(ref_group)
k.weight = 100
predictions <- TransferData(anchorset = transfer_anchors, 
                            weight.reduction = "cca", 
                            k.weight = k.weight, sd.weight = 1,
                            refdata = FetchData(obj_ref, ref_group)[[1]],
                            dims = 1:30 )
obj_query <- AddMetaData(obj_query, metadata = predictions)

# Plot transfer, predicted.id -----------------------------------------------------------


obj_query_md = obj_query@meta.data %>%
  as.data.frame()

obj_ref_md = obj_ref@meta.data %>%
  as.data.frame()

query_group = "cluster_int_sub2"
query_group_sym = as.symbol(query_group)
obj_query_md_stat = obj_query_md %>% group_by(!!query_group_sym, predicted.id) %>%
  summarize(n=n()) %>%
  ungroup() %>% 
  group_by(!!query_group_sym) %>%
  mutate(n_frac = n / sum(n))
obj_query_md_stat_mat = obj_query_md_stat %>% select_at(vars(one_of(query_group, "predicted.id", "n_frac"))) %>%
  pivot_wider(names_from = query_group, values_from = "n_frac", values_fill = list(n_frac = 0)) %>% 
  column_to_rownames("predicted.id") %>% 
  as.matrix()


mat_cur = obj_query_md_stat_mat
cols = colorRamp2(breaks = c(0, 1), colors=c("white", "black"))

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 


ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))

ct_levels = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
sc_levels = c("L2/3 IT", "L4", "L5 IT", "L5 PT", "L6 IT", "L6 CT", "L6b")
mat_cur = mat_cur[sc_levels,ct_levels]
hm = Heatmap(mat_cur, 
             cluster_rows = F, 
             clustering_method_rows = "ward.D",
             cluster_columns=F,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in")
             )

hm
pdf(file.path(dev_out_sub_dir, sprintf("transfer_predicted_id_kweight%s.pdf", k.weight)), width= target_width+4, height = target_height + 4)
print(hm)
dev.off()

# Plot transfer, prob -----------------------------------------------------

obj_query_md = obj_query@meta.data %>%
  as.data.frame() 

obj_ref_md = obj_ref@meta.data %>%
  as.data.frame() %>% 
  rownames_to_column()

query_group = "cluster_int_sub2"
query_group_sym = as.symbol(query_group)

obj_query_md_stat = obj_query_md %>% select_at(vars(contains("prediction.score"), one_of(c(query_group)), one_of("rowname"))) %>%
  select(-contains("GABA"), -rowname) %>% 
  pivot_longer(contains("prediction.score"), names_to = "class", values_to = "value") %>% 
  mutate(class = sub("prediction.score.", "", class)) %>%
  filter(class != "max") %>%
  group_by(!!query_group_sym, class) %>%
  summarize(value_mean = mean(value))

obj_query_md_stat_mat = obj_query_md_stat %>% select_at(vars(one_of(query_group, "class", "value_mean"))) %>%
  pivot_wider(names_from = query_group, values_from = "value_mean", values_fill = list(value_mean = 0)) %>% 
  column_to_rownames("class") %>% 
  as.matrix()

mat_cur = obj_query_md_stat_mat

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 

ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))

ct_levels = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")
sc_levels = c("L2.3.IT", "L4", "L5.IT", "L5.PT", "L6.IT", "L6.CT", "L6b")
mat_cur = mat_cur[sc_levels,ct_levels]

cols = colorRamp2(breaks = c(0, max(mat_cur)), colors=c("white", "black"))

hm = Heatmap(mat_cur, 
             cluster_rows = F, 
             clustering_method_rows = "complete",
             clustering_method_columns = "complete",
             cluster_columns=F,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in")
             )
hm
pdf(file.path(dev_out_sub_dir, sprintf("transfer_mean_prediction_kweight%s.pdf", k.weight)), width= target_width+4, height = target_height + 4)
print(hm)
dev.off()
