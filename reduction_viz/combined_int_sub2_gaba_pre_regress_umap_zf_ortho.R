library(Seurat)
library(tidyverse)
library(qs)
library(reticulate)
use_condaenv("anaconda3")
library(cowplot)
theme_set(theme_cowplot())

library(ComplexHeatmap)
library(circlize)
source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")


# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_GABA.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_int_sub2_gaba_pre_regress_umap_zf_ortho"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")


# Load data ---------------------------------------------------------------

dims = 20
n.neighbors = 30
min.dist = 1
spread = 2.5
reduction.name = sprintf("dims%snn%smindist%sspread%s", dims, n.neighbors, min.dist, spread)

redo = F
open_orig = F
if (redo) {
  if (open_orig) {
    obj_int_filt = qread(data_fname)
  } else {
    obj_int_filt = qread(data_out_obj_fname)
  }
  res_to_use = "cluster_int_sub2"
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
  idents_to_use = Idents(obj_int_filt)
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>%
    FindVariableFeatures(nfeatures=3000) %>%
    RunPCA()
  
  obj_int_filt = obj_int_filt %>%
    RunUMAP(
      seed.use = 42L,
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
DimPlot(obj_int_filt, reduction="umap")
DimPlot(obj_int_filt, reduction=reduction.name, group.by="species")
  
# Cluster UMAP ------------------------------------------------------------


res_to_use1 = "cluster_int_sub"
cats = c("species", "position2", res_to_use1)
for ( ca in cats ) {
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none"#, 
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
      legend.position="none"
      ) +
  scale_color_manual(values=position_colors2[c("hvc", "ra")])
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=5, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=5, base_asp =1 )



# Cluster UMAP, cluster_int_sub2 ------------------------------------------------------------


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

celltype_colors = c('GABA-1-1' = "#2CA05A",
                    'GABA-1-2' = "#4bce7f",
                    'GABA-2' = "#E85F5F",
                    'GABA-3' = "#AE5FA2",
                    'GABA-4' = "#D3135F",
                    'GABA-5-1' = "#27788c",
                    'GABA-5-2' = "#37ABC8",
                    'GABA-5-3' = "#73c4d9",
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


# Heatmap -----------------------------------------------------------------

genes_to_plot = c("GAD1", "FOXP2", "NPY", "NXPH1", "TTLL5", "PENK", "LAMP5", "LHX8", "CALB2")
ct_order = c("GABA-1", "GABA-2", "GABA-3", "GABA-4", "GABA-5", "GABA-6", "GABA-8", "GABA-7")
assay_to_use = "SCT"
DefaultAssay(obj_int_filt) = assay_to_use

idents_to_use = "cluster_int_sub"
Idents(obj_int_filt) = FetchData(obj_int_filt, idents_to_use)
mat_avg = AverageExpression(obj_int_filt, assays=assay_to_use)
mat_avg = log1p(mat_avg[[assay_to_use]])
mat_avg_cur = t(mat_avg[genes_to_plot,ct_order])
#mat_avg_cur = t(apply(mat_avg_cur, 1, scale_01))

height = 1
width = height * ncol(mat_avg_cur)  / nrow(mat_avg_cur)
cols = colorRamp2(breaks=c(0,max(mat_avg_cur)), colors=c("white", "black"))
hm = Heatmap(mat_avg_cur, cluster_rows=F, cluster_columns = F, col=cols,
             height = unit(height, "in"),
             width = unit(width, "in"))
hm

pdf(file.path(figures_dir, c("marker_heatmap.pdf")), height = height+2, width = width +2)
print(hm)
dev.off()

# Plot gene expression ----------------------------------------------------

cat("UMAP\n")
genes_to_plot = c(
  "LAMP5", "PVALB", "SST", "NPY", 
  "CRH", "PENK", 
  "FOXP2", 
  "SOX4", "DCX",
  "LHX8",
  "CEMIP", "NXPH1",
  "VIP", "NOS1", "CCK", "RELN", "CALB2",
  "CNTN6", "PTN", "CXCL12",
  "SLC17A6", "CDH9", "LSAMP"
)


em = Embeddings(obj_int_filt, reduction.name)
xmin = min(em[,1])
assays_to_use = c( "SCT")
for (assay_to_use in assays_to_use) {
  DefaultAssay(obj_int_filt) = assay_to_use
for (gene in genes_to_plot) {
  gg = FeaturePlot(obj_int_filt, reduction=reduction.name, feature=gene) +
    labs(title="") + 
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position="none", 
      panel.border = element_rect(fill=NA, color=1)
    ) + 
    annotate(
      x=max(em[,1]),
      y=max(em[,2]),
      label=gene,
      geom="text",
      fontface="italic",
      size=8,
      hjust=1
    ) + 
    coord_equal()
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", assay_to_use, gene)), gg, base_height=7, base_asp =1.1 )
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", assay_to_use, gene)), gg, base_height=7, base_asp =1.1 )
}
}
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use1)
gg = DotPlot(obj_int_filt, assay = "SCT", features =  genes_to_plot) + 
  coord_flip()
gg
