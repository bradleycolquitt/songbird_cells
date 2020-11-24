

source("~/data2/rstudio/birds/utils/go.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
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
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggrepel)
theme_set(theme_cowplot())

ag_genes = get_axon_guidance_genes()
tx = import("~/data2/assembly/lonStrDom1/ncbi/lonStrDom1.ncbi.apollo6/merge_inter_apollo_processed_mt_gffread.gtf")
tx_mt = tx[seqnames(tx)=="MT"]


# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_X.qs")

## Input dir

dev_in_dir = file.path(dir_root, "devin_combined", "songbird_cells")
tree_dir = file.path(dev_in_dir, "trees")
script_name = "celltypes_hclust_hvc_ra_x_gaba_int_sub2_regress_zf_ortho"
tree_dir = file.path(tree_dir, script_name)
data_in_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
out_dir = file.path(dev_out_dir, "dataset_comparison")
script_name = "hvc_ra_x_gabaergic_ag_genes"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir)

## Markers
assay_to_use = "RNA"
markers_fname = file.path(tree_dir, sprintf("hvc_ra_x_gaba_rna_markers.rds"))


# Load data ---------------------------------------------------------------

obj_int_filt = qread(data_in_obj_fname)
res_to_use = "cluster_orig1"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)

# LGE only ----------------------------------------------------------------

markers_lge_fname = file.path(out_dir, sprintf("hvc_ra_x_gaba_lge_sct_markers.rds"))
cells = Cells(obj_int_filt)[obj_int_filt$cluster_orig1 %in% c("GABA-1", "0", "1", "2", "3", "4", "9")]
obj_int_filt_lge = subset(obj_int_filt, cells=cells)

obj_int_filt_lge$cluster_lge = case_when(obj_int_filt_lge$cluster_orig1 %in% c("0", "1", "2", "3", "4") ~ "MSN",
                                         obj_int_filt_lge$cluster_orig1 %in% c("9") ~ "PN",
                                         obj_int_filt_lge$cluster_orig1 %in% "GABA-1" ~ "GABA-1")
Idents(obj_int_filt_lge) = obj_int_filt_lge$cluster_lge
redo_markers = T
if (redo_markers) {
  DefaultAssay(obj_int_filt_lge) = "SCT"
  
  markers_lge = FindAllMarkers(obj_int_filt_lge, 
                               assay="SCT", 
                               test.use = "wilcox",
                               min.pct = .2, only.pos=T)
  saveRDS(markers_lge, markers_lge_fname)
} else {
  markers_lge = readRDS(markers_lge_fname)
}

obj_int_filt_lge_avg = AverageExpression(obj_int_filt_lge, assays="SCT")
# Heatmap, ag -------------------------------------------------------------


n_top = 15
to_plot1 = markers_lge %>% 
  filter(p_val_adj < .05) %>%
  filter(gene %in% ag_genes$external_gene_name)

to_plot = to_plot1 %>% 
  group_by(cluster) %>%
  top_n(n_top, avg_logFC)
#to_plot = c("EPHA4", "EFNA5", "EFNA3", "ROBO1", "ROBO2", "CXCR4", "ACKR3", "CXCL12", "NRG1", "NRG3", "ERBB4")
to_plot = unique(c(to_plot$gene, c("ERBB4", "NRG1")))




mat_a = log1p(obj_int_filt_lge_avg[[1]])
mat_a = mat_a[to_plot,]
mat_a_scale = t(apply(mat_a, 1, function(x) (x - min(x)) / (max(x) - min(x))))
#mat_a_scale = t(apply(mat_a, 1, function(x) scale(x, scale = F)) )
colnames(mat_a_scale) = c("MSN", "PN", "GABA-1")
mat_a_scale = mat_a_scale[,c("GABA-1", "PN", "MSN")]
blues = brewer_pal(palette="Blues")(9)
colors = colorRamp2(breaks=c(0,1), colors = blues[c(1,9)])


target_width = .25
target_height = target_width * nrow(mat_a) / ncol(mat_a) 


hm = Heatmap(mat_a_scale, 
             #col=colors,
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             clustering_method_rows = "ward.D",
             cluster_columns = F,
             show_row_dend = F,
             show_column_dend = F,
             #cell_fun = star_fun,
             row_names_gp = gpar(fontsize=5),
             column_names_gp = gpar(fontsize=5)
)
hm


target_height = .25
target_width = target_height * nrow(mat_a) / ncol(mat_a) 

hm = Heatmap(t(mat_a_scale), 
             col=colors,
             width = unit(target_width, "in"),
             height = unit(target_height, "in"),
             clustering_method_columns = "ward.D",
             cluster_rows = F,
             show_row_dend = F,
             show_column_dend = F,
             #cell_fun = star_fun,
             row_names_gp = gpar(fontsize=5),
             column_names_gp = gpar(fontsize=5)
)
hm

pdf(file.path(out_dir, "ag_heatmap_horiz.pdf"), width=10, height=4)
print(hm)
dev.off()
