library(Seurat)
library(tidyverse)
library(qs)

library(cowplot)

theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")

# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_GLUT.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_glut_ra_umap_zf_ortho"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")

res_to_use = "cluster_int_sub2"
dims = 5
n.neighbors = 30L
min.dist = 0.3
reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

redo = F
if (redo) {

  # Load data ---------------------------------------------------------------
  
  obj_int_filt = qread(data_fname)
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)

  cells = colnames(obj_int_filt)[grepl("GLUT", Idents(obj_int_filt))]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  obj_int_filt = subset(obj_int_filt, subset=position2=="ra")
  

  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>% 
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() 
  obj_int_filt = RunUMAP(obj_int_filt,
                         dims = 1:dims,
                         min.dist = min.dist,
                         n.neighbors = n.neighbors,
                         reduction.name=reduction.name
  )

  obj_int_filt@reductions$umap = NULL
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}

# Plot clusters -----------------------------------------------------------


res_to_use1 = "cluster_int_sub2"
cats = c("species", "position2", res_to_use1)
for ( ca in cats ) {
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none", panel.border = element_rect(fill=NA, color=1)) +
  labs(x="UMAP1", y="UMAP2")  + 
  coord_equal()
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", ca, reduction.name)), gg, base_height=5, base_width=5 )
}

# Plot gene expression ----------------------------------------------------

n_to_plot = 5

genes_to_plot = unique(c( "CRHR2", "PCP4", "NR4A2", "SORCS1", "LRP1B", "VIPR1", "FOSL2", "HOMER1", "EGR1", "NDNF", "PVALB",
                          "MTNR1B", "ADAMTS2", "ADAMTS18", "NFATC1", "COL6A3",
                          "percent.mito", "nCount_RNA"))

assays_to_use = c("SCT")
for (assay_to_use in assays_to_use) {
DefaultAssay(obj_int_filt) = assay_to_use
em = Embeddings(obj_int_filt, reduction.name) 
xmin = min(em[,1])
for (gene in genes_to_plot) {
  gg = FeaturePlot(obj_int_filt, reduction=reduction.name, feature=gene) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(), title = element_blank(),
      legend.position="none", 
      panel.border = element_rect(fill=NA, color=1)
    ) + 
    annotate(
      x=min(em[,1]) + .1*(max(em[,1]) - min(em[,1]) ),
      y=max(em[,2]),
      label=gene,
      geom="text",
      fontface="italic",
      size=6
    )
  gg
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s_%s.pdf", assay_to_use, gene, reduction.name)), gg, base_height=4, base_asp =1.1 )
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s_%s.png", assay_to_use, gene, reduction.name)), gg, base_height=4, base_asp =1.1 )
}
}


# Dotplot -----------------------------------------------------------------
s
genes_to_plot = c("SLC17A6",  "COL6A3", "ADAMTS18", "NFATC1")
genes_to_plot = c("SLC17A6", "NFATC1", "ADAMTS18", "COL6A3")
Idents(obj_int_filt) = obj_int_filt$cluster_int_sub2
#Idents(obj_int_filt) = factor(Idents(obj_int_filt), levels=rev(c("RA_GLUT-3", "RA_GLUT-2", "RA_GLUT-1")), labels=c("RA_Glut-3", "RA_Glut-2", "RA_Glut-1"))
gg = DotPlot(obj_int_filt, features=genes_to_plot, assay="SCT") + 
  scale_size_area() + 
  labs(x="", y = "")+ coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_manual.pdf"), gg, base_height=2.5)

