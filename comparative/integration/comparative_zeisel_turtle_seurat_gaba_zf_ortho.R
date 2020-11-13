
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
script_name = "comparative_zeisel_turtle_seurat_gaba_zf_ortho"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

## Data dir
tree_dir = file.path(dev_dir, "trees")
script_data_name = "celltypes_hclust_gaba_int_sub2_regress_zf_ortho"
tree_dir = file.path(tree_dir, script_data_name)
data_in_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")
data_out_obj_fname = file.path(tree_dir, "obj_intergrated_subclustered_gaba_zf_ortho.qs")

## Compare dir
compare_data_dirs = list(
  mm = "~/data/scrna_datasets/zeisel_2018/",
  tse = "~/data/scrna_datasets/tosches_2018/")


# Params ------------------------------------------------------------------

params = list(k.anchor = 5,
              k.filter = 100,
              k.score = 30,
              max.features=100)
dev_out_sub_dir= file.path(out_dir, 
                           sprintf("anchor%s_filter%s_score%s_maxfeatures%s", params[["k.anchor"]], params[["k.filter"]], params[["k.score"]], params[["max.features"]]))
dir.create(dev_out_sub_dir, recursive=T)

# Load orthologs ----------------------------------------------------------


to_use = c("zf", "mm", "tse")

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

# Load comparison data, mm ----------------------------------------------------

compare_prefix_1 = "l6_r2_cns_neurons_seurat"
compare_prefix_1_out = "l6_r2_cns_neurons_seurat_GABA_Acetyl_zf_ortho"
compare_1_out_fname = file.path(compare_data_dirs[["mm"]], sprintf("%s.qs", compare_prefix_1_out))

redo = F
if (redo) {
  na = "mm"
  egg_l_cur = egg_l %>% filter(species==na) %>%
    filter(!is.na(predicted_name)) %>%
    filter(predicted_name!="")
  
  obj_fname = file.path(compare_data_dirs[["mm"]], sprintf("%s.qs", compare_prefix_1))
  compare_obj = qread(obj_fname)
  compare_obj$idents = Idents(compare_obj)
  nt_select = "GABA|Acetyl"
  cells = Cells(compare_obj)[grepl(nt_select, compare_obj$Neurotransmitter)]
  x = subset(compare_obj, cells=cells)
  rm(compare_obj)
  dat = GetAssayData(x, assay="RNA", slot="counts")
  dat_new = dat[intersect(rownames(dat), toupper(egg_l_cur$gene)),]
  new_rnames = egg_l_cur$predicted_name[match(rownames(dat_new), toupper(egg_l_cur$gene))]
  dimnames(dat_new) = list(new_rnames, colnames(dat_new))
  assay_new =  CreateAssayObject(counts=dat_new)
  x[["ortho"]] = assay_new
  DefaultAssay(x) = "ortho"
  
  options(future.globals.maxSize=10000 * 1024^2)
  plan(multiprocess(workers=3))
  compare_obj = SCTransform(x, assay = "ortho",
                  return.only.var.genes = F,
                  vars.to.regress = c("MitoRiboRatio"))

  rm(x)

  qsave(compare_obj, compare_1_out_fname)
} else {
  compare_obj = qread(compare_1_out_fname)
}

compare_obj_1 = compare_obj
compare_obj_1$species = "mm"
compare_obj_1$cluster_int_sub2 = Idents(compare_obj_1)

regions_to_use = c("Hippocampus,Cortex",
                   "Hippocampus",
                   "Olfactory bulb",
                   "Striatum dorsal",
                   "Striatum ventral",
                   "Striatum dorsal, Striatum ventral",
                   "Striatum dorsal, Striatum ventral,Amygdala",
                   "Pallidum")

mm_anno = FetchData(compare_obj_1, c("Description", "idents"))
mm_anno = mm_anno %>% distinct(Description, idents)

idents_levels = c("MSN1", "MSN2", "MSN3", "MSN4", "MSN5", "MSN6", 
                  "OBINH1", "OBINH2", "OBINH3", "OBINH4", "OBINH5",
                  "OBDOP1", "OBDOP2", "OBNBL1", "OBNBL2", "OBNBL3", "OBNBL4", "OBNBL5",
                  "TEINH1", "TEINH2", "TEINH3", "TECHO", 
                  "TEINH21", "TEINH19", "TEINH20", "TEINH13", "TEINH17", "TEINH18",
                  "TEINH14", "TEINH15", "TEINH16", "TEINH7", "TEINH8", "TEINH4", "TEINH5", "TEINH6",
                  "TEINH12", "TEINH9", "TEINH10", "TEINH11")

cells = Cells(compare_obj_1)[compare_obj_1$Region %in% regions_to_use]
compare_obj_1 = subset(compare_obj_1, cells=cells)

# Load comparison data, tse ----------------------------------------------------

compare_2_fname = file.path(compare_data_dirs[["tse"]], "turtle.neurons.v3.rds")
fp = "turtle.neurons.gaba.v3.SCT_zf_ortho"
compare_2_out_fname = file.path(compare_data_dirs[["tse"]], sprintf("%s.qs", fp))

redo = F
if (redo) {
  
  na = "tse"
  egg_l_cur = egg_l %>% filter(species==na) %>%
    filter(!is.na(predicted_name)) %>%
    filter(predicted_name!="")
  
  compare_obj = readRDS(compare_2_fname)
  cells = Cells(compare_obj)[grepl("i", compare_obj$clusters)]
  compare_obj = subset(compare_obj, cells=cells)
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
                            vars.to.regress = c("percent.mito"))
  qsave(compare_obj, compare_2_out_fname)
} else {
  compare_obj = qread(compare_2_out_fname)
}

compare_obj_2 = compare_obj
compare_obj_2$species = "tse"
compare_obj_2$cluster_int_sub2 = compare_obj_2$clusters

class_df = data.frame(cluster_int_sub2 = c(paste0("i0",  1:9), paste0("i", 10:18)),
                      class = c("OB", "septum", "septum", "striatum", "amygdala", "amygdala", 
                                "SST", "SST", "SST", 'SST',
                                "PV-like", "PV-like", "PV-like",
                                "HTR3A VIP-like", "HTR3A VIP-like", "HTR3A VIP-like",
                                "HTR3A Reln", "HTR3A Reln"),
                      class2 = factor(c("LGE-derived", "septum", "septum", "LGE-derived", "LGE-derived", "LGE-derived",
                                 "MGE-derived", "MGE-derived", "MGE-derived", "MGE-derived",
                                 "MGE-derived", "MGE-derived", "MGE-derived",
                                 "CGE-derived", "CGE-derived", "CGE-derived",
                                 "CGE-derived", "CGE-derived"),
                                 levels = c("LGE-derived", "septum", "MGE-derived", "CGE-derived")))


rm(compare_obj)

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

# Integrate ---------------------------------------------------------------

fname_anchors = file.path(dev_out_sub_dir, "obj_anchors.qs")
fname_int = file.path(dev_out_sub_dir, "obj_integrated_clustered.qs")

redo = T
redo_anchors = T

obj_list = list(zf = obj_int_filt, mm = compare_obj_1, tse = compare_obj_2)
if (redo) {
  options(future.globals.maxSize= 15000 * 1024^2)
  
  obj_features <- SelectIntegrationFeatures(object.list = obj_list,
                                            assay=rep("SCT", times=length(obj_list))
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
  qsave(obj_int, fname_int)
} else {
  obj_int = qread(fname_int)
}


redo_reduction = T

dims = 30
n.neighbors = 30L
min.dist = .3
reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

if (redo_reduction) {
  obj_int <- RunPCA(obj_int,  verbose = FALSE)
  obj_int = RunUMAP(obj_int,
                    dims = 1:dims,
                    min.dist = min.dist,
                    n.neighbors = n.neighbors,
                    reduction.name=reduction.name
  )
  ress = seq(.2, 1, .2)

}

# Dimplot -----------------------------------------------------------------

redo_plot = T
if (redo_plot) {
  reduction = reduction.name
  cats = c("orig.ident", "position", "species", "cluster_int_sub2")
  
  em = Embeddings(obj_int, reduction) %>% as.data.frame() %>% rownames_to_column()
  colnames(em) = c("rowname", "UMAP_1", "UMAP_2")
  md = FetchData(obj_int, cats) %>% rownames_to_column()
  em = em %>% left_join(md)
  walk(cats, function(ca) {
    
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
    png(file.path(dev_out_sub_dir, sprintf("umap_int_%s_by_%s.png", reduction, ca)), width=4, height=4, units = "in", res = 300)
    print(gg)
    dev.off()
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_%s_by_%s_legend.pdf", reduction, ca)), gg_leg, base_width=4)
    
    gg = ggplot(em_ca, aes(UMAP_1, UMAP_2)) + 
      geom_text(aes_string(label=ca, color=ca)) + 
      theme_void() + 
      scale_color_manual(values=col)  + 
      coord_equal()
    gg = gg + theme(legend.position="none")
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_%s_by_%s_label.pdf", reduction, ca)), gg, base_width=4)
  })
  
    cats = c("orig.ident", "position", "species")
  walk(cats, function(ca) {
    cas = unique(FetchData(obj_int, ca)[[1]])
    col = qualpalr::qualpal(length(cas))$hex
    names(col) = cas
    gg = ggplot(em, aes_string("UMAP_1", "UMAP_2", color=ca)) +
      geom_point(alpha=.1) + 
      theme_void() + 
      scale_color_manual(values=col)  + 
      facet_wrap(ca)
    print(gg)
    save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_%s_split-by_%s.png", reduction, ca)), gg, base_width = 4)
  })
  
  ## Split by species, color by cluster
  ca = "cluster_int_sub2"
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
    facet_grid(.~species) +
    coord_equal()
  gg_leg = get_legend(gg)
  gg = gg + theme(legend.position="none")
  png(file.path(dev_out_sub_dir, sprintf("umap_int_%s_split-by-species_%s.png", reduction, ca)), width=12, height=4, units = "in", res = 300)
  print(gg)
  dev.off()
  save_plot(file.path(dev_out_sub_dir, sprintf("umap_int_%s_split-by-species_by_%s_legend.pdf", reduction, ca)), gg_leg, base_width=4)
}


# Transfer data -----------------------------------------------------------


obj_refs = obj_list[2:3]
obj_query = obj_list[[1]]

## causes error in FindTransferAnchors
obj_refs = map(obj_refs, function(obj_ref) {
  obj_ref@reductions$ica = NULL
  obj_ref@reductions$tsne = NULL
  obj_ref
})

transfer_anchors_fname = file.path(dev_out_sub_dir, "transfer_anchors.qs")
redo_transfer = T
if (redo_transfer) {
  
  transfer_anchors_list = map(obj_refs, function(obj_ref) {
    obj_features <- SelectIntegrationFeatures(object.list = list(obj_ref, obj_query),
                                              assay=rep("SCT", times=2)
    )
    transfer_anchors <- FindTransferAnchors(reference = obj_ref,
                                            query = obj_query, 
                                            features = obj_features,
                                            reference.assay = "SCT",
                                            query.assay = "SCT",
                                            reduction = "cca",
                                            k.anchor = params[["k.anchor"]],
                                            k.filter = params[["k.filter"]],
                                            k.score = params[["k.score"]],
                                            max.features = params[["max.features"]],
                                            dims = 1:30)
    transfer_anchors
  }) %>% set_names(names(obj_refs))
  qsave(transfer_anchors_list, transfer_anchors_fname)
} else {
  transfer_anchors_list = qread(transfer_anchors_fname)  
}

ref_group = "cluster_int_sub2"
ref_group_sym = as.symbol(ref_group)
predictions =  map_dfc(names(transfer_anchors_list), function(na) {
  transfer_anchors = transfer_anchors_list[[na]]
  obj_ref = obj_refs[[na]]
  preds = TransferData(anchorset = transfer_anchors, 
               weight.reduction = "cca", 
               k.weight = 100, 
               sd.weight = 1,
               refdata = FetchData(obj_ref, ref_group)[[1]],
               dims = 1:30)
  colnames(preds) = paste(na, colnames(preds), sep=".")
  preds
})
rownames(predictions) = Cells(obj_query)

obj_query <- AddMetaData(obj_query, metadata = predictions)


# Plot transfer, predicted.id -----------------------------------------------------------

walk(names(obj_refs), function(species) {
  print(species)
species_not = names(obj_refs)[names(obj_refs) != species]
obj_query_md = obj_query@meta.data %>%
  as.data.frame() %>%
  select(-matches(sprintf("^%s", species_not)))

colnames(obj_query_md) = sub(sprintf("^%s\\.", species), "", colnames(obj_query_md))

obj_ref = obj_refs[[species]]
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

mat_cur_rmax = apply(mat_cur, 1, max)
thresh = .1
mat_cur = mat_cur[mat_cur_rmax>thresh,]

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 

if (species=="mm") {
  ann_row = obj_ref_md %>% distinct(!!ref_group_sym, .keep_all=T) %>% 
    select(!!ref_group_sym, Taxonomy_group, Region) %>% 
    mutate(!!ref_group_sym := factor(!!ref_group_sym, levels=idents_levels)) %>%
    arrange(!!ref_group_sym) %>% 
    column_to_rownames(ref_group)
  
  ann_row = ann_row[rownames(ann_row) %in% rownames(mat_cur),]
  mat_cur = mat_cur[match(rownames(ann_row), rownames(mat_cur)),]
  cols = colorRamp2(breaks = c(0, 1), colors=c("white", "black"))
  
  regions = na.omit(unique(ann_row$Region))
  region_colors = pal_aaas()(length(regions))
  names(region_colors) = regions
  
  taxs = unique(ann_row$Taxonomy_group)
  tax_colors = qualpalr::qualpal(length(taxs))$hex
  names(tax_colors) = taxs
  ha = rowAnnotation(df = ann_row, 
                     annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                     col = list(Region = region_colors, Taxonomy_group = tax_colors),
                     annotation_width=unit(target_width*.2, "in"),
                     simple_anno_size_adjust=T,
                     annotation_legend_param = list(
                       Region = list(gp = gpar(fontsize=6)),
                       Taxonomy_group = list(gp = gpar(fontsize=6)))
  )
  
} else if (species=="tse") {
  ann_row = obj_ref_md %>% distinct(!!ref_group_sym, .keep_all=T) %>% 
    select(!!ref_group_sym) %>% left_join(class_df) %>% 
    arrange(class2, class) %>% 
    column_to_rownames(ref_group)
  
  ann_row = ann_row[rownames(ann_row) %in% rownames(mat_cur),]
  mat_cur = mat_cur[match(rownames(ann_row), rownames(mat_cur)),]
  cols = colorRamp2(breaks = c(0, 1), colors=c("white", "black"))
  classes = unique(class_df$class)
  class_colors = qualpalr::qualpal(length(classes))$hex
  names(class_colors) = classes
  
  classes2 = unique(class_df$class2)
  class2_colors = pal_aaas()(length(classes2))
  names(class2_colors) = classes2
  
  ha = rowAnnotation(df = ann_row, 
                     annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                     col = list(class = class_colors, class2 = class2_colors),
                     annotation_width=unit(target_width*.2, "in"),
                     simple_anno_size_adjust=T,
                     annotation_legend_param = list(
                       class = list(gp = gpar(fontsize=6)),
                       class2 = list(gp = gpar(fontsize=6)))
  )
  
  
  
}

ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))

hc_r = hclust(dist(mat_cur, method="euclidean"), method="ward.D")
hm = Heatmap(mat_cur, 
             cluster_rows = F,
             cluster_columns=F,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             left_annotation = ha)
hm
pdf(file.path(dev_out_sub_dir, sprintf("transfer_predicted_id_%s.pdf", species)), width= target_width+4, height = target_height + 4)
print(hm)
dev.off()
})

# Plot transfer, prob -----------------------------------------------------


walk(names(obj_refs), function(species) {
  print(species)
  species_not = names(obj_refs)[names(obj_refs) != species]
  obj_query_md = obj_query@meta.data %>%
    as.data.frame() %>%
    select(-matches(sprintf("^%s", species_not)))
  
  colnames(obj_query_md) = sub(sprintf("^%s\\.", species), "", colnames(obj_query_md))
  
  obj_ref = obj_refs[[species]]
  obj_ref_md = obj_ref@meta.data %>%
    as.data.frame()
  query_group = "cluster_int_sub2"
  query_group_sym = as.symbol(query_group)
  
obj_query_md_stat = obj_query_md %>% select_at(vars(contains("prediction.score"), one_of(c(query_group)), one_of("rowname"))) %>%
  select(-contains("GABA"), -rowname) %>% 
  pivot_longer(contains("prediction.score"), names_to = "class", values_to = "value") %>% 
  mutate(class = sub("prediction.score.", "", class)) %>%
  filter(class != "max") %>%
  group_by(!!query_group_sym, class) %>%
  summarize(value_mean = mean(value))# %>%

obj_query_md_stat_mat = obj_query_md_stat %>% select_at(vars(one_of(query_group, "class", "value_mean"))) %>%
  pivot_wider(names_from = query_group, values_from = "value_mean", values_fill = list(value_mean = 0)) %>% 
  column_to_rownames("class") %>% 
  as.matrix()

mat_cur = obj_query_md_stat_mat


mat_cur_rsums = rowSums(mat_cur)
thresh = .25
mat_cur = mat_cur[mat_cur_rsums>thresh,]

target_width = 1
target_height = target_width * nrow(mat_cur) / ncol(mat_cur) 


if (species=="mm") {
  ann_row = obj_ref_md %>% distinct(!!ref_group_sym, .keep_all=T) %>% 
    select(!!ref_group_sym, Taxonomy_group, Region) %>% 
    column_to_rownames(ref_group)
  
  ann_row = ann_row[rownames(mat_cur),]

  regions = na.omit(unique(ann_row$Region))
  region_colors = pal_aaas()(length(regions))
  names(region_colors) = regions
  
  taxs = unique(ann_row$Taxonomy_group)
  tax_colors = qualpalr::qualpal(length(taxs))$hex
  names(tax_colors) = taxs
  ha = rowAnnotation(df = ann_row, 
                     annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                     col = list(Region = region_colors, Taxonomy_group = tax_colors),
                     annotation_width=unit(target_width*.2, "in"),
                     simple_anno_size_adjust=T,
                     annotation_legend_param = list(
                       Region = list(gp = gpar(fontsize=6)),
                       Taxonomy_group = list(gp = gpar(fontsize=6)))
  )
  
} else if (species=="tse") {
  ann_row = obj_ref_md %>% distinct(!!ref_group_sym, .keep_all=T) %>% 
    select(!!ref_group_sym) %>% left_join(class_df) %>%
    column_to_rownames(ref_group)
  
  ann_row = ann_row[rownames(mat_cur),]
  classes = unique(class_df$class)
  class_colors = qualpalr::qualpal(length(classes))$hex
  names(class_colors) = classes
  
  classes2 = unique(class_df$class2)
  class2_colors = pal_aaas()(length(classes2))
  names(class2_colors) = classes2
  
  ha = rowAnnotation(df = ann_row, 
                     annotation_name_gp =  gpar(fontsize=6, fontface="bold"),
                     col = list(class = class_colors, class2 = class2_colors),
                     annotation_width=unit(target_width*.2, "in"),
                     simple_anno_size_adjust=T,
                     annotation_legend_param = list(
                       class = list(gp = gpar(fontsize=6)),
                       class2 = list(gp = gpar(fontsize=6)))
  )
  
}


cols = colorRamp2(breaks = c(0, max(mat_cur)), colors=c("white", "black"))

ht_opt(legend_labels_gp = gpar(fontsize=6), 
       legend_title_gp = gpar(fontsize=6, fontface="bold"),
       legend_grid_width = unit(target_width / ncol(mat_cur), "in"),
       legend_grid_height = unit(target_height / nrow(mat_cur), "in"),
       heatmap_row_names_gp = gpar(fontsize=6),  
       heatmap_column_names_gp = gpar(fontsize=6))


hc_r = hclust(dist(mat_cur, method="correlation"), method="average")
hm = Heatmap(mat_cur, 
             cluster_rows = hc_r, 
             cluster_columns=T,
             col=cols, 
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             left_annotation = ha)

hm
pdf(file.path(dev_out_sub_dir, sprintf("transfer_predicted_mean_%s.pdf", species)), width= target_width+4, height = target_height + 4)
print(hm)
dev.off()
})
