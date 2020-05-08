library(Seurat)
library(tidyverse)
library(qs)
library(cowplot)
theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

red_dir = file.path(dev_dir, "reduction_viz")
script_name = "combined_int_sub2_glut_pre_umap_regress"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_fname= file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")
data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")

# Load data ---------------------------------------------------------------

res_to_use = "cluster_int_sub2"

dims = 12
n.neighbors = 30L
min.dist = .3
reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

redo = F
if (redo) {
  

  obj_int_filt = qread(data_fname)
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)

  cells = colnames(obj_int_filt)[grepl("Glut|^Pre", Idents(obj_int_filt))]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  
  cells = colnames(obj_int_filt)[!(is.na(obj_int_filt$position2) | obj_int_filt$position2=="arco")]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  
  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>% 
    FindVariableFeatures() %>%
    ScaleData(
    )
  obj_int_filt = RunPCA(obj_int_filt)
  obj_int_filt = RunUMAP(obj_int_filt,
                         dims = 1:dims,
                         min.dist = min.dist,
                         n.neighbors = n.neighbors,
                         reduction.name=reduction.name)
  qsave(obj_int_filt, data_out_obj_fname)
  
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Plot --------------------------------------------------------------------


cats = c("species", "position2", res_to_use)
for ( ca in cats ) {
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none", 
       ) +
  labs(x="", y="UMAP2") 
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", ca, reduction.name)), gg, base_height=7, base_asp =1)
}

ca = "position2"
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none", 
        ) +
  scale_color_manual(values=position_colors2[c("hvc", "ra")])
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1 )


