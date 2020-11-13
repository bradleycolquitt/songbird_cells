library(Seurat)
library(tidyverse)
library(qs)

library(cowplot)
library(fgsea)
theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")
#source("~/data2/rstudio/birds/utils/common_aesthetics.R")


# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_int_sub2_glut_hvc_regress_umap_sub-IEG_whitney_zf_ortho"
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

## IEGs
gray_gmt_fname = "~/data/gene_lists/gray_2018_arg.xls.gmt"
gray_gmt = gmtPathways(gray_gmt_fname)
names(gray_gmt) = make.names(names(gray_gmt))
gray_gmt$ARGs = c(gray_gmt$ARGs, "FOSL2", "HOMER1")


## Singing
sing_fname = "~/data/gene_lists/whitney_science_2014/Supplementary Table S8 - Singing Regulated Transcripts v53.xlsx"
sing = readxl::read_xlsx(sing_fname, skip=2)
colnames(sing) = make.names(colnames(sing))
sing_hvc = sing %>% filter(HVC.p<.1) %>% 
  filter(!is.na(Symbol)) 

all_iegs = unique(c(gray_gmt$ARGs, sing_hvc$Symbol))

# Load data ---------------------------------------------------------------

redo = T
if (redo) {
  

  
  obj_int_filt = qread(data_fname)
  Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)

  cells = colnames(obj_int_filt)[grepl("Glut", Idents(obj_int_filt))]
  obj_int_filt = subset(obj_int_filt, cells = cells)
  obj_int_filt = subset(obj_int_filt, subset=position2=="hvc")

  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
  obj_int_filt = obj_int_filt %>% 
    FindVariableFeatures() 
  
  var_feats = VariableFeatures(obj_int_filt)
  var_feats_filt = var_feats[!(var_feats %in% all_iegs)]
  VariableFeatures(obj_int_filt) = var_feats_filt
  
  obj_int_filt = obj_int_filt %>%
    ScaleData()
    
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


# Average expression, ARG plot ------------------------------------------------------

mat_avg = AverageExpression(obj_int_filt, assay="SCT", slot="counts")
mat_avg = log1p(mat_avg[["SCT"]])
mat_avg_l = mat_avg %>% rownames_to_column(var = "gene") %>% pivot_longer(cols=contains("HVC"), names_to="cluster", values_to="value") %>%
  mutate(is_ieg = gene %in% gray_gmt$rapid.PRGs)

mat_avg_l_cur = mat_avg_l %>% filter(is_ieg)
gg = ggplot(mat_avg_l %>%filter(is_ieg), aes(cluster, value, color=cluster)) + 
  geom_point(position=position_jitter(width=.4) ) 
gg


# Dotplot, IEG ------------------------------------------------------------

iegs_cur = intersect(rownames(obj_int_filt), c(gray_gmt$rapid.PRGs, "BDNF", "FOSL2", "HOMER1"))
gg = DotPlot(obj_int_filt, features = iegs_cur, assay = "SCT") + coord_flip() + 
  theme(axis.text.x = element_text(angle=90, hjust=1))
gg

# Average expression, IEG plot ------------------------------------------------------

mat_avg = AverageExpression(obj_int_filt, assay="SCT")
mat_avg = log1p(mat_avg[["SCT"]])

ieg_manual = c("NR4A3", "EGR1", "EGR2", "")
mat_avg_l = mat_avg %>% rownames_to_column(var = "gene") %>% pivot_longer(cols=contains("HVC"), names_to="cluster", values_to="value") %>%
  mutate(is_ieg = gene %in% gray_gmt$ARGs)

gg = ggplot(mat_avg_l %>%filter(is_ieg), aes(cluster, value, color=cluster)) + 
  geom_point(position=position_jitter(width=.4) ) 
gg



# Plot gene expression ----------------------------------------------------


genes_to_plot = c("FOSL2", "BDNF", "HOMER1", "UCHL1", 
                  "GRM8", "FOXP1", "UTS2B", "NTS",
                  "GRIA4", "SLIT3", "MYO5B", "CHRM2", "CACNA1G",
                  "TESC", "CDH9", "ADARB2", "GPC5", "DPP10", "GFRA1", "FOXP2", "SCUBE1", "GFRA1",
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
