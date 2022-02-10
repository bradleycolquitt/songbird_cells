library(Seurat)
library(tidyverse)
library(qs)

library(cowplot)
theme_set(theme_cowplot())

library(ComplexHeatmap)
library(circlize)
source("~/data2/rstudio/birds/utils/scRNA.R")

source("~/data2/rstudio/birds/utils/stats.R")

#source("~/data2/rstudio/birds/utils/common_aesthetics.R")


# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_int_sub2_glut_hvc_regress_umap"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")
# 
# dir_root = "~/sdd/data2/rstudio/birds/scRNA"
# dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
# dev_data_dir = file.path(dir_root, "devin_combined", "data")
# dev_out = file.path(dev_dir, "preprocessing", "integrate", "SCT_regress", "hvc")
# dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s", 5, 200, 30))
# 
# red_dir = file.path(dev_dir, "reduction_viz")
# script_name = "combined_int_sub2_glut_hvc_regress_umap"
# red_dir = file.path(red_dir, script_name)
# dir.create(red_dir, recursive = T)
# figures_dir = file.path(red_dir)
# dir.create(figures_dir)
# 
# data_fname= file.path(dev_out_sub_dir, "obj_integrated_hvc_subclustered.qs")
# data_out_obj_fname = file.path(red_dir, "obj_integrated_hvc_subclustered.qs")


res_to_use = "cluster_int_sub2"

dims = 8
n.neighbors = 30L
min.dist = .3
reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

redo = F
if (redo) {
  
  # Load data ---------------------------------------------------------------
  
  obj_int_filt = qread(data_fname)
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)

  cells = colnames(obj_int_filt)[grepl("Glut", Idents(obj_int_filt))]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  obj_int_filt = subset(obj_int_filt, subset=position2=="hvc")
 
  # md = FetchData(obj_int_filt, c(res_to_use)) %>%
  #   rownames_to_column() 
  # 
  # # md_un = md %>%
  # #   distinct(cluster_sub_hvc, .keep_all=T) %>% 
  # #   filter(grepl("Glut", cluster_sub_hvc)) %>%
  # #   mutate(tmp1 = map(cluster_sub_hvc, ~unlist(str_split(.x, "-")))) %>%
  # #   mutate(tmp1a = map_chr(tmp1, 1),
  # #          tmp1b = as.integer(map_chr(tmp1, 2))) %>% 
  # #   mutate(new_index = tmp1b + 1L) %>% 
  # #   mutate(cluster_sub_hvc1 = paste(tmp1a, new_index, sep="-")) %>%
  # #   select(cluster_sub_hvc, cluster_sub_hvc1)
  # # # 
  # # md = md %>% left_join(md_un) %>%
  # #   column_to_rownames()
  # obj_int_filt = AddMetaData(obj_int_filt, md)
  
  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>% 
    FindVariableFeatures() %>%
    ScaleData(
      #vars.to.regress = c("orig.ident")
              )
  obj_int_filt = RunPCA(obj_int_filt)
  
  obj_int_filt = RunUMAP(obj_int_filt,
                         dims = 1:dims,
                         min.dist = min.dist,
                         n.neighbors = n.neighbors,
                         reduction.name=reduction.name
  )

#  obj_int_filt = seurat_run_alra(obj_int_filt, assay="SCT")
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
}


# Plot --------------------------------------------------------------------


res_to_use1 = "cluster_int_sub2"
cats = c("species", "position2", "orig.ident", res_to_use1)
for ( ca in cats ) {
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none", panel.border = element_rect(fill=NA, color=1)) +
  labs(x="UMAP1", y="UMAP2") + 
  coord_equal() 
#gg = AugmentPlot(gg)
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", ca, reduction.name)), gg, base_height=7, base_asp =1.1 )
}

cats = c("species", "position2", "orig.ident", res_to_use1)
for ( ca in cats ) {
  gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          #legend.position="none", 
          panel.border = element_rect(fill=NA, color=1)) +
    labs(x="UMAP1", y="UMAP2") + 
    coord_equal() 
  #gg = AugmentPlot(gg)
  gg
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s_legend.pdf", ca, reduction.name)), gg, base_height=7, base_asp =1.2 )
}

# Heatmap -----------------------------------------------------------------

genes_to_plot = c("SLC17A6", "GFRA1", "GRIA4", "SCUBE1", "CACNA1G", "BDNF")
ct_order = c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-5", "HVC_Glut-2", "HVC_Glut-3")
assay_to_use = "SCT"
DefaultAssay(obj_int_filt) = assay_to_use

mat_avg = AverageExpression(obj_int_filt, assays=assay_to_use)
mat_avg = log1p(mat_avg[[assay_to_use]])
mat_avg_cur = mat_avg[genes_to_plot,ct_order]
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


genes_to_plot = c("FOSL2", "BDNF", "HOMER1", "UCHL1", 
                  "GRM8", "FOXP1", "UTS2B", "NTS",
                  "GRIA4", "SLIT3", "MYO5B", "CHRM2", "CACNA1G",
                  "TESC", "CDH9", "ADARB2", "GPC5", "DPP10", "GFRA1", "FOXP2", "SCUBE1", "GFRA1", "DUSP1", "FOS", "EGR1",
                  "GAP43")
assay_to_use = "SCT"
DefaultAssay(obj_int_filt) = assay_to_use
em = Embeddings(obj_int_filt, reduction.name)
xmin = min(em[,1])
for (gene in genes_to_plot) {
  gg = FeaturePlot(obj_int_filt, reduction=reduction.name, feature=gene) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position="none", 
          #panel.border = element_rect(fill=NA, color=1)
          ) +
    labs(x="UMAP1", y="UMAP2", title = "") + 
    coord_equal() +
    annotate(
      x=min(em[,1]),
      y=max(em[,2]),
      label=gene,
      geom="text",
      fontface="italic",
      size=4
    )
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s_%s.pdf", assay_to_use, gene, reduction.name)), gg, base_height=4, base_asp =1.1 )
  save_plot(file.path(figures_dir, sprintf("umap_%s_%s_%s.png", assay_to_use, gene, reduction.name)), gg, base_height=4, base_asp =1.1 )
}





# Dotplot -----------------------------------------------------------------

#Idents(obj_int_filt) = paste(toupper(obj_int_filt$position2), obj_int_filt$cluster_int_sub2, sep="_")
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
Idents(obj_int_filt) = factor(Idents(obj_int_filt), levels=rev(c("HVC_Glut-1", "HVC_Glut-4",  "HVC_Glut-2", 
                                                                "HVC_Glut-3", "HVC_Glut-5")))

genes_to_plot = c("GFRA1", "GRIA4", "SCUBE1", "SLIT3", "CACNA1G")
gg = DotPlot(obj_int_filt, assay = "SCT", features =  rev(genes_to_plot), col.min = 0, dot.min=.1) +
  theme(
    axis.text.x = element_text(angle=30, face="italic", hjust=1),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") 
# coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_vert.pdf"), gg, base_height=length(unique(Idents(obj_int_filt)))/2, 
          base_width = length(genes_to_plot))


genes_to_plot = c("FOSL2", "BDNF", "HOMER1", "UCHL1", 
                  "GRM8", "FOXP1", "UTS2B", "NTS",
                  "GRIA4", "SLIT3", "MYO5B", "CHRM2", "CACNA1G",
                  "TESC", "CDH9", "ADARB2", "GPC5", "DPP10", "GFRA1", "FOXP2", "SCUBE1", 
                  "DUSP1", "FOS", "EGR1",
                  "GAP43")
gg = DotPlot(obj_int_filt, assay = "SCT", features =  rev(genes_to_plot), col.min = 0, dot.min=.1) +
  theme(
    axis.text.x = element_text(angle=30, face="italic", hjust=1),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") 
# coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_vert_full.pdf"), gg, base_height=length(unique(Idents(obj_int_filt)))/2, 
          base_width = length(genes_to_plot))

# Dotplot, cholinergic recepors -----------------------------------------------------------------

#Idents(obj_int_filt) = paste(toupper(obj_int_filt$position2), obj_int_filt$cluster_int_sub2, sep="_")
Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
Idents(obj_int_filt) = factor(Idents(obj_int_filt), levels=rev(c("HVC_Glut-1", "HVC_Glut-4",  "HVC_Glut-2", 
                                                                 "HVC_Glut-3", "HVC_Glut-5")))

DefaultAssay(obj_int_filt) = "SCT"
genes_to_plot = sort(grep("^CHR[M|N]", rownames(obj_int_filt), value=T))
gg = DotPlot(obj_int_filt, assay = "SCT", features =  rev(genes_to_plot), col.min = 0, dot.min=.1) +
  scale_size_area() + 
  theme(
    axis.text.x = element_text(angle=30, face="italic", hjust=1),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") 
# coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_CHR_vert.pdf"), gg, base_height=length(unique(Idents(obj_int_filt)))/2, 
          base_width = length(genes_to_plot)/2)

genes_to_plot = rev(sort(grep("^CHRM", rownames(obj_int_filt), value=T)))
gg = DotPlot(obj_int_filt, assay = "SCT", features =  rev(genes_to_plot), col.min = 0, dot.min=.1) +
  scale_size_area() + 
  theme(
    axis.text.x = element_text(angle=30, face="italic", hjust=1),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") 
# coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_CHRM_vert.pdf"), gg, base_height=length(unique(Idents(obj_int_filt)))/2, 
          base_width = length(genes_to_plot)/2)

Idents(obj_int_filt) = factor(Idents(obj_int_filt), levels=rev(levels(Idents(obj_int_filt))))
gg = DotPlot(obj_int_filt, assay = "SCT", features =  genes_to_plot, col.min = 0, dot.min=.1) +
  scale_size_area() + 
  theme(
    axis.text.y = element_text(angle=0, face="italic", hjust=1),
    axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
    legend.position="bottom",
    legend.direction="horizontal"
  ) + 
  labs(x="", y="") +
  coord_flip()
gg
save_plot(file.path(figures_dir, "dotplot_CHRM_horiz.pdf"), gg, base_width=length(unique(Idents(obj_int_filt)))/2, 
          base_height = length(genes_to_plot) * .75)


save_plot(file.path(figures_dir, "dotplot_CHRM_horiz_legend.pdf"), get_legend(gg), base_width=20, 
          base_height = 6)
# genes_to_plot = c("FOSL2", "BDNF", "HOMER1", "UCHL1", 
#                   "GRM8", "FOXP1", "UTS2B", "NTS",
#                   "GRIA4", "SLIT3", "MYO5B", "CHRM2", "CACNA1G",
#                   "TESC", "CDH9", "ADARB2", "GPC5", "DPP10", "GFRA1", "FOXP2", "SCUBE1", 
#                   "DUSP1", "FOS", "EGR1",
#                   "GAP43")
# gg = DotPlot(obj_int_filt, assay = "SCT", features =  rev(genes_to_plot), col.min = 0, dot.min=.1) +
#   theme(
#     axis.text.x = element_text(angle=30, face="italic", hjust=1),
#     legend.position="bottom",
#     legend.direction="horizontal"
#   ) + 
#   labs(x="", y="") 
# # coord_flip()
# gg
# save_plot(file.path(figures_dir, "dotplot_vert_full.pdf"), gg, base_height=length(unique(Idents(obj_int_filt)))/2, 
#           base_width = length(genes_to_plot))
# 


