library(Seurat)
library(tidyverse)
library(qs)
library(cowplot)
theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

red_dir = file.path(dev_dir, "reduction_viz")
script_name = "combined_int_sub2_gaba_pre_regress_umap"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_fname= file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")
data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered_gaba.qs")


# Load data ---------------------------------------------------------------

dims = 15
n.neighbors = 30
min.dist = .5
spread = 2.5
reduction.name = sprintf("dims%snn%smindist%sspread%s", dims, n.neighbors, min.dist, spread)

redo = T
if (redo) {
  obj_int_filt = qread(data_fname)
  
  res_to_use = "cluster_int_sub2"
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  idents_to_use = Idents(obj_int_filt)
  
  ## Select GABA
  idents_to_remove = c("ra_GABA-Pre")
  idents_to_use = idents_to_use[grepl("GABA", idents_to_use)]
  obj_int_filt = subset(obj_int_filt, idents=idents_to_use)
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>%
    FindVariableFeatures(nfeatures=3000) %>%
    RunPCA()
  
  obj_int_filt = obj_int_filt %>%
    RunUMAP(
      dims = 1:dims,
      min.dist = min.dist,
      n.neighbors = n.neighbors, 
      spread=spread,
      reduction.name=reduction.name
    ) 
  
  table(obj_int_filt$cluster_int_sub)
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

DimPlot(obj_int_filt, reduction=reduction.name)
DimPlot(obj_int_filt, reduction=reduction.name, group.by="species")

# Cluster UMAP ------------------------------------------------------------

res_to_use1 = "cluster_int_sub2"
cats = c("species", "position2", res_to_use1)
for ( ca in cats ) {
  gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position="none"
    ) +
    labs(x="", y="")
  gg
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1.1 )
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1.1 )
}

celltype_colors = c('GABA-1' = "#2CA05A",
                    'GABA-2' = "#E85F5F",
                    'GABA-3' = "#AE5FA2",
                    'GABA-4' = "#D3135F",
                    'GABA-5' = "#37ABC8",
                    'GABA-6' = "#D3C230",
                    'GABA-7' = "#005FD2",
                    'GABA-8' = "#C87137",
                    'GABA-Pre' = "#1AD7AE")

ncells = length(Cells(obj_int_filt))
ca = res_to_use1
em = Embeddings(obj_int_filt, reduction.name)
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none"
  ) +
  scale_color_manual(values=celltype_colors)

gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=5, base_asp =1)
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=5, base_asp =1)

gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none"
  ) +
  scale_color_manual(values=celltype_colors)

gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_no_label.pdf", reduction.name, ca)), gg, base_height=5, base_asp =1)
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_no_label.png", reduction.name, ca)), gg, base_height=5, base_asp =1)




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
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=5, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=5, base_asp =1 )
