library(Seurat)
library(tidyverse)
library(qs)
library(seriation)
library(proxy)
library(cowplot)

theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")


# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_glut_pre_umap_zf_ortho"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

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
  
  cells = colnames(obj_int_filt)[!is.na(obj_int_filt$position2)]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  
  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>% 
    FindVariableFeatures() %>%
    ScaleData()
  obj_int_filt = RunPCA(obj_int_filt)
  
  obj_int_filt = RunUMAP(obj_int_filt,
                         dims = 1:dims,
                         min.dist = min.dist,
                         n.neighbors = n.neighbors,
                         reduction.name=reduction.name
  )
  

  
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
        legend.position="none"
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
        legend.position="none"
        ) +
  scale_color_manual(values=position_colors2[c("hvc", "ra")])

gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1 )

# Plot gene expression, neurogenesis ----------------------------------------------------

cat("UMAP, neurogenesis\n")
source("~/data2/rstudio/birds/utils/gene_lists.R")

figures_ng_dir = file.path(figures_dir, "neurogenesis")
dir.create(figures_ng_dir)

genes_to_plot = c(genes_b1_cells, genes_c_cells, genes_qNSC, genes_aNSC, c("ASCL1", "NR2E1", "SOX4", "SOX11", "DCX",
                                                                           "NECTIN3",
                                                                           "SLC17A6", 
                                                                           "NEUROG2",
                                                                           "NEUROD1", "PPP1R1B"
                                                                           ))
genes_to_plot = unique(genes_to_plot)
#"TBR2", "DLL3", "EGFR1",
genes_to_plot = genes_to_plot[!(genes_to_plot %in% c( "SRRT", "HDAC1", "NEUROD1"))]
#assays_to_use = c( "SCT", "RNA")
obj_int_filt_sub = subset(obj_int_filt, idents=ct_order)
dat = AverageExpression(obj_int_filt, assays = "SCT")
dat = log1p(dat[["SCT"]])
dat = dat[,!grepl("RA", colnames(dat))]
ct_order = c("Pre-1", "Pre-2", "Pre-3", "Pre-4", "HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3")
dat = dat[intersect(genes_to_plot, rownames(dat)), ct_order]
dat_sums = rowSums(dat)
dat_maxes = apply(dat, 1, max)
plot(density(dat_maxes))
dat = dat[dat_maxes>.25,]
gene_dist = dist(dat, method="euclidean")
hc_gene = reorder(hclust(gene_dist, method="average"), gene_dist)
plot(hc_gene)
genes_to_plot_order = hc_gene$labels[hc_gene$order]


Idents(obj_int_filt_sub) = factor(Idents(obj_int_filt_sub), levels=ct_order)
gg = DotPlot(obj_int_filt_sub, assay = "SCT", features =  rev(genes_to_plot_order), col.min = 0, dot.min=.1, scale = T) +
  scale_size_area() + 
  theme(
    axis.text.x = element_text(angle=90, face="italic", hjust=1, vjust=.5),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") 
# coord_flip()
gg
save_plot(file.path(figures_ng_dir, "dotplot_vert.pdf"), gg,
          base_height=length(unique(Idents(obj_int_filt_sub))) / 2.5, 
          base_width = length(genes_to_plot_order) /3.5)

save_plot(file.path(figures_ng_dir, "dotplot_vert_legend.pdf"), get_legend(gg),
          base_height=length(unique(Idents(obj_int_filt_sub))), 
          base_width = length(genes_to_plot_order))
