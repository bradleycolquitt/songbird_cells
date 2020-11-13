library(Seurat)
library(tidyverse)
library(ggdendro)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(scales)
library(ggsci)
library(dendextend)
library(scales)
library(colorspace)
library(ComplexHeatmap)
library(pvclust)
library(future)
library(seriation)
plan(sequential)

theme_set(theme_cowplot())
source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
source("~/data2/rstudio/birds/utils/stats.R")

# Parameters --------------------------------------------------------------

hclust_method = "ward.D"
sub_dir = sprintf("%s", hclust_method)

# Directories -------------------------------------------------------------


## Data dir
dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
tree_dir = file.path(dev_out_dir, "trees")
script_data_dir = file.path(tree_dir, "celltypes_hclust_glut_int_sub2_regress_zf_ortho")
script_data_sub_dir = file.path(script_data_dir, sub_dir)

script_name =  "celltypes_dendro_markers_violin_glut_regress_split"
tree_dir = file.path(tree_dir, script_name)
tree_sub_dir = file.path(tree_dir,  sprintf("%s", hclust_method))
dir.create(tree_sub_dir, recursive = T)
figures_dir = tree_sub_dir

data_out_obj_fname = file.path(script_data_dir, "obj_integrated_subclustered.qs")
data_out_avg_fname = file.path(script_data_dir, "average_expr.qs")
data_out_pv_fname = file.path(script_data_sub_dir, "pvclust.qs")


# dir_root = "~/data2/rstudio/birds/scRNA"
# dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
# dev_data_dir = file.path(dir_root, "devin_combined", "data")
# dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
# dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))
# 
# tree_dir = file.path(dev_dir, "trees")
# script_data_dir = file.path(tree_dir, "celltypes_hclust_glut_int_sub2_regress")
# script_data_sub_dir = file.path(script_data_dir, sub_dir)
# script_name = "celltypes_dendro_markers_violin_glut_regress_split"
# tree_dir = file.path(script_data_dir, script_name, sub_dir)
# figures_dir = file.path(tree_dir, "figures")
# dir.create(figures_dir, recursive = T)
# 
# data_fname= file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")
# data_out_obj_fname = file.path(script_data_dir, "obj_integrated_subclustered_glut.qs")
# data_out_avg_fname = file.path(script_data_dir, "average_expr.qs")
# data_out_pv_fname = file.path(script_data_sub_dir, "pvclust.qs")

# Load data ---------------------------------------------------------------


obj_int_filt = qread(data_out_obj_fname)
obj_int_avg1 = qread(data_out_avg_fname)
pv = qread(data_out_pv_fname)

res_to_use = "cluster_int_sub2"
Idents(obj_int_filt) = obj_int_filt$cluster_int_sub2


# Marker id -----------------------------------------------------------------

assay_to_use = "SCT"
markers_fname = file.path(dev_data_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))
markers_pos_fname = file.path(dev_data_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s_pos.rds", assay_to_use))

redo_markers = F
if (redo_markers) {
  Idents(obj_int_filt) = obj_int_filt@meta.data[,res_to_use]
  DefaultAssay(obj_int_filt) = assay_to_use
  plan(multiprocess(workers = 2))
  markers_int = FindAllMarkers(obj_int_filt, 
                               test.use = "wilcox", 
                               #logfc.threshold = log(2), 
                               #min.pct = .5,
                               #min.diff.pct = .25,
                               max.cells.per.ident = 200,
                               only.pos = F)
  saveRDS(markers_int, markers_fname)
  
  ## Within position
  markers_int_pos = map(c("hvc", "ra"), function(pos) {
    cells = colnames(obj_int_filt)[obj_int_filt$position2==pos]
    tmp = subset(obj_int_filt, cells=cells)
    FindAllMarkers(tmp, 
                   test.use = "wilcox", 
                   #logfc.threshold = log(2), 
                   #min.pct = .5,
                   #min.diff.pct = .25,
                   max.cells.per.ident = 200,
                   only.pos = F)
  })
  saveRDS(markers_int_pos, markers_pos_fname)
} else {
  markers_int = readRDS(markers_fname)
  markers_int_pos = readRDS(markers_pos_fname)
}


# Idents(obj_int_filt) = obj_int_filt$cluster_int_sub2
# 
# #DefaultAssay(obj_int) = "SCT"
# assay_to_use = "SCT"
# markers_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s.rds", assay_to_use))
# 
# markers_int = readRDS(markers_fname)
# 
# markers_pos_fname = file.path(dev_out_sub_dir, sprintf("marker_genes_cluster_int_sub2_glut_%s_pos.rds", assay_to_use))
# 
# markers_int_pos = readRDS(markers_pos_fname)
 markers_int_pos = markers_int_pos %>% bind_rows()


# Calculate gene specificity ----------------------------------------------

 obj_int_avg1_spec = na.omit(apply(obj_int_avg1, 1, function(x) sum(calc_gene_specificity(x)))) 
 
 obj_int_avg1_m = melt(as.matrix(obj_int_avg1))
 colnames(obj_int_avg1_m) = c("gene", "cluster", "value")
 
 obj_int_avg1_m = obj_int_avg1_m %>% mutate(position2 = str_extract(cluster, "^[A-Z]+")) %>%
   group_by(position2, gene) %>%
   mutate(value_spec = calc_gene_specificity(value+.1))
 

 
 regions = c("HVC", "RA")
 walk(regions, function(reg) {
   
   nspec = 40
   obj_int_avg1_m_top = obj_int_avg1_m %>% 
     filter(position2==reg) %>% 
     group_by(cluster) %>%
     top_n(nspec, value_spec) %>%
     ungroup() %>%
     distinct(gene, .keep_all=T)
   
   obj_int_avg1_m1 = obj_int_avg1_m %>%
     filter(position2==reg) %>% 
     group_by(gene) %>%
     summarize(value_spec_sum = sum(value_spec))

   obj_int_avg1_m1_top = obj_int_avg1_m1 %>%
     top_n(200, value_spec_sum)

   
   nsig = 25
   markers_int_pos_top = markers_int_pos %>%
     filter(grepl(reg, cluster)) %>% 
     filter(avg_logFC>0) %>%
     #mutate(sign = avg_logFC>0) %>%
     group_by(cluster) %>%
     top_n(-1 * nsig, p_val_adj) %>%
     top_n(nsig, avg_logFC) %>%
     ungroup() %>%
     distinct(gene, .keep_all=T)
   
   genes = Reduce(intersect, list(markers_int_pos_top$gene, obj_int_avg1_m_top$gene, obj_int_avg1_m1_top$gene))
   genes = unique(c(genes, switch(reg, "HVC" = c("UTS2B", "NTS", "ALDH1A2", "UCHL1", "NFATC1"),
                                      "RA" = c("PCP4", "ADAMTS18"))))
   obj_int_avg1_filt = obj_int_avg1[genes,grep(reg, colnames(obj_int_avg1))]
   d = dist(obj_int_avg1_filt)
   hc = hclust(d, method="ward.D")
   plot(hc)
   hc = reorder(hc, dist = d)
   genes = hc$labels[hc$order]
   
   cells = Cells(obj_int_filt)[grepl(reg, obj_int_filt$cluster_int_sub2)]
   obj_int_filt_cur = subset(obj_int_filt, cells=cells)
   Idents(obj_int_filt_cur) = factor(Idents(obj_int_filt_cur), 
                                     levels= c("HVC_Glut-1", "HVC_Glut-4", "HVC_Glut-2", "HVC_Glut-5", "HVC_Glut-3",
                                               "RA_Glut-1", "RA_Glut-2", "RA_Glut-3"))
   gg = StackedVlnPlot(obj_int_filt_cur, genes) +
     theme(axis.text.x = element_text(size=6, angle=90, hjust=1)) &
     theme(axis.text.y= element_text(size=5),
           axis.title.y = element_text(angle=0, hjust = 0, size = 5, vjust=0.5),
           axis.line = element_line(size=.25),
           axis.ticks = element_line(size=.25)) 
   
   height = length(genes) / 10 + 1
   print(height)
 #  width = length(unique(Idents(obj_int_filt_cur))) 
   save_plot(file.path(figures_dir, sprintf("violin_%s.pdf", reg)), gg, base_height=height, base_width=3)
   save_plot(file.path(figures_dir, sprintf("violin_%s_slim.pdf", reg)), gg, base_height=height, base_width=2)
 })
 