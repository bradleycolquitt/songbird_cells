
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
script_name = "comparative_saunders_seurat_glut_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name =  "celltypes_hclust_glut_int_sub2_regress"
tree_dir = file.path(tree_dir, script_data_name)
data_in_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut.qs")
data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_glut_zf_ortho.qs")

## Comparison 
compare_dir = "~/data/scrna_datasets/dropviz/"
compare_prefix = "saunders_pallial_obj"

compare_fname_prefix = list(fc = "F_GRCm38.81.P60Cortex_noRep5_FRONTALonly",
                            pc = "F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly",
                            hc = "F_GRCm38.81.P60Hippocampus",
                            str = "F_GRCm38.81.P60Striatum")
compare_dat_fnames = map(compare_fname_prefix, ~file.path(compare_dir, sprintf("%s.raw.dge.txt.gz", .x)))
compare_subcluster_fnames = map(compare_fname_prefix, ~file.path(compare_dir, sprintf("%s.subcluster.assign.RDS", .x)))

tfs = get_tf_genes_human()

# Params ------------------------------------------------------------------

params = list(k.anchor = 5,
              k.filter = 200,
              k.score = 30,
              max.features=200)
dev_out_sub_dir= file.path(out_dir, 
                           sprintf("anchor%s_filter%s_score%s_maxfeatures%s", params[["k.anchor"]], params[["k.filter"]], params[["k.score"]], params[["max.features"]]))
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

# Load comparison ---------------------------------------------------------

compare_prefix1 = "saunders_pallial_obj_zf_ortho"
compare_sub_fname = file.path(compare_dir, sprintf("%s.qs", compare_prefix1))

redo = F
if (redo) {
  
  md = read.csv(file.path(compare_dir, "annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv"))
  md_glut = md %>% filter(grepl("Slc17", class_marker)) %>%
    filter(!grepl("Gad", class_marker))
  na = "mm"
  egg_l_cur = egg_l %>% filter(species==na) %>%
    filter(!is.na(predicted_name)) %>%
    filter(predicted_name!="")
  
  dats = map(compare_dat_fnames, loadSparseDge )
  names(dats) = names(compare_dat_fnames)
  clusters = map(compare_subcluster_fnames, readRDS)
  clusters_df = map_dfr(names(clusters), ~data.frame(cell=names(clusters[[.x]]), subcluster = clusters[[.x]], tissue=.x))
  clusters_df = clusters_df %>% 
    mutate(tissue = toupper(tissue)) %>% 
    mutate(tissue_subcluster = paste(tissue, subcluster, sep="_"))
  
  clusters_df_glut = clusters_df %>% inner_join(md_glut %>% select(tissue_subcluster))
  
  subsample = .1
  
  clusters_df_glut_samp = clusters_df_glut %>% 
    group_by(tissue_subcluster) %>% 
    sample_frac(size=subsample)
  clusters_df_glut_samp = clusters_df_glut_samp %>% left_join(md_glut)
  clusters_df_glut_samp_md = clusters_df_glut_samp %>% column_to_rownames("cell")
  dats_glut = map(names(dats), function(x) {
    clusters_cur = clusters_df_glut_samp %>% filter(tissue == toupper(x))
    dat_cur = dats[[x]]
    dat_cur[,clusters_cur$cell]
  })
  
  common_genes = Reduce(intersect, map(dats_glut, rownames))
  dats_glut = map(dats_glut, function(dat) {
    dat = dat[common_genes,]
    dat
  })
  
  names(dats_glut) = names(dats)
  compare_objs = map(names(dats_glut), function(na) {
    obj = CreateSeuratObject(counts =dats_glut[[na]],
                             project = sprintf("saunders_%s", na), 
                             assay="RNA")
    obj$percent.mito = PercentageFeatureSet(obj, pattern = "^mt-")
    obj$tissue = toupper(na)
    obj
  })
  compare_obj = merge(compare_objs[[1]], compare_objs[2:4])
  compare_obj = AddMetaData(compare_obj, clusters_df_glut_samp_md)
  
  dat = GetAssayData(compare_obj, assay="RNA", slot="counts")
  dat_new = dat[intersect(rownames(dat), egg_l_cur$gene),]
  new_rnames = egg_l_cur$predicted_name[match(rownames(dat_new), egg_l_cur$gene)]
  dimnames(dat_new) = list(new_rnames, colnames(dat_new))
  assay_new =  CreateAssayObject(counts=dat_new)
  compare_obj[["ortho"]] = assay_new
  DefaultAssay(compare_obj) = "ortho"
  
  options(future.globals.maxSize=10000 * 1024^2)
  plan(multiprocess(workers=3))
  compare_obj = SCTransform(compare_obj, assay = "ortho",
                            return.only.var.genes = F,
                            vars.to.regress = c("percent.mito"))
  
  qsave(compare_obj, compare_sub_fname)
} else {
  compare_obj = qread(compare_sub_fname)
}

Idents(compare_obj) = FetchData(compare_obj, "tissue_subcluster")
compare_obj$idents = FetchData(compare_obj, "tissue_subcluster")
md = FetchData(compare_obj, c( "class_marker", "tissue", "idents", "common_name", "full_name", "type_marker")) %>%
  distinct(idents, .keep_all=T)

compare_obj = FindVariableFeatures(compare_obj, nfeatures=6000)
compare_obj_var = VariableFeatures(compare_obj)

compare_obj$cluster_int_sub2 = compare_obj$tissue_subcluster

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

# Filter ------------------------------------------------------------------

common_genes = Reduce(intersect, list(rownames(obj_int_filt), rownames(compare_obj)))

obj_list = list(zf = obj_int_filt, mm = compare_obj)
obj_list = map(obj_list, function(x) {
  x[common_genes,]
})


# Integrate ---------------------------------------------------------------

fname_anchors = file.path(dev_out_sub_dir, "obj_anchors.qs")
fname_int = file.path(dev_out_sub_dir, "obj_integrated_clustered.qs")

redo = T
redo_anchors = T


if (redo) {
  options(future.globals.maxSize= 15000 * 1024^2)
  
  obj_features <- SelectIntegrationFeatures(object.list = obj_list,
                                            assay=rep("SCT", times=length(obj_list)), fvf.nfeatures = 3000
  )
  
  obj_list <- PrepSCTIntegration(object.list = obj_list, 
                                 assay = "SCT", 
                                 anchor.features = obj_features, 
                                 verbose = TRUE)
  
  if (redo_anchors) {
    obj_anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                          assay = rep("SCT", times=length(obj_list)),
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

obj_int$cluster_int_sub3 =  if_else(grepl("_[0-9]+-[0-9]", obj_int$cluster_int_sub2), 
                                    str_extract(obj_int$cluster_int_sub2, "[A-Z]+_[0-9]"), 
                                    obj_int$cluster_int_sub2)

# Dimplot -----------------------------------------------------------------


redo_plot = T
if (redo_plot) {
  reduction = "umap"
  cats = c("orig.ident", "position", "species", "cluster_int_sub3")
  
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
  
  walk(cats, function(ca) {
    print(ca)
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
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_split-by-species_legend.pdf", ca)), gg_leg, base_width=4, ncol=2)
    
    gg = ggplot(em_ca, aes(UMAP_1, UMAP_2)) + 
      geom_text(aes_string(label=ca, color=ca)) + 
      theme_void() + 
      scale_color_manual(values=col)  + 
      coord_equal()+ 
      facet_grid(.~species)
    gg = gg + theme(legend.position="none")
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_by_%s_split-by-species_label.pdf", ca)), gg, base_width=4, ncol=2)
  })
}
  

# Transfer data -----------------------------------------------------------


obj_ref = obj_list[[2]]
obj_query = obj_list[[1]]

transfer_anchors_fname = file.path(dev_out_sub_dir, "transfer_anchors.qs")
redo_transfer = T
if (redo_transfer) {
  
  
  obj_features <- SelectIntegrationFeatures(object.list = list(obj_ref, obj_query),
                                            assay=rep("SCT", times=2)
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


# Annotation --------------------------------------------------------------

md_ann = FetchData(compare_obj, c("tissue_subcluster", "type_marker", "common_name", "full_name", "tissue" )) %>% distinct_all()

md_ann = md_ann %>%
  mutate(class = case_when(grepl("Deep|Layer\\s*[6/5]", common_name) ~ "Deep cortical",
                           grepl("Superficial|2/3", common_name) ~ "Upper cortical",
                           grepl("Claust", common_name) ~ "Claustrum",
                           grepl("Cajal", common_name) ~ "CR",
                           grepl("Subplate", common_name, ignore.case=T) ~ "Subplate",
                           grepl("Subiculum", common_name, ignore.case = T) ~ "Subiculum",
                           grepl("Retrosplen", common_name) ~ "RSP",
                           grepl("Entorhi", common_name, ignore.case = T) ~ "Entorhinal",
                           grepl("Dentate", common_name, ignore.case = T) ~ "DG",
                           grepl("CA1", common_name) ~ "CA1",
                           grepl("CA3", common_name) ~ "CA3",
                           grepl("CA2", common_name) ~ "CA2"),
  ) %>%
  filter(!is.na(class)) %>%
  filter(!grepl("IEG", common_name))

# Plot transfer, predicted.id -----------------------------------------------------------
ct_levels = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3", "RA_Glut-1", "RA_Glut-2", "RA_Glut-3")

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
mat_cur = mat_cur[,ct_levels]

rsums = rowSums(mat_cur)

cols = colorRamp2(breaks = c(0, 1), colors=c("white", "black"))

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 

ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))


hm = Heatmap(mat_cur, 
             cluster_rows = T, 
             clustering_method_rows = "ward.D",
             cluster_columns=F,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
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
mat_cur = mat_cur[,ct_levels]

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 

cols = colorRamp2(breaks = c(0, max(mat_cur)), colors=c("white", "black"))
ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))

hm = Heatmap(mat_cur, 
             cluster_rows = T, 
             clustering_method_rows = "ward.D",
             cluster_columns=F,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in")
             )
hm
pdf(file.path(dev_out_sub_dir, sprintf("transfer_mean_prediction_kweight%s.pdf", k.weight)), width= target_width+4, height = target_height + 4)
print(hm)
dev.off()
