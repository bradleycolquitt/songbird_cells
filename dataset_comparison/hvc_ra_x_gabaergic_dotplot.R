
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
script_name = "hvc_ra_x_gabaergic_dotplot"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir)

## Markers
assay_to_use = "RNA"
markers_fname = file.path(tree_dir, sprintf("hvc_ra_x_gaba_rna_markers.rds"))

# Load finch data ---------------------------------------------------------

x_key = data.frame(cluster_orig = seq(0,22),
                   label = c("MSN (1)", "MSN (2)", "MSN (3)", "MSN (4)", "MSN (5)", "Astro (6)", "Astro (7)", "GLUT (8)", "Oligo (9)",
                             "PN (10)", "GLUT (11)", "Unk PPP1R1B (12)", "GABA SST (13)", "GABA PVALB (14)", "GABA PVALB (15)", "Micro (16)",
                             "GABA (17)", "Endo (18)", "GABA (19)", "OPC (20)", "GABA (21)", "GABA (22)", "GABA CHAT (23)"))
x_key_to_use = x_key %>% filter(grepl("MSN|GABA|PN", label))

redo = F
if (redo) {
  load(data_fname)
  obj_int = HVCX_seurat_REGION_res1.1
  
  cells = Cells(obj_int)[obj_int$cluster_orig %in% x_key_to_use$cluster_orig]
  cells2 = Cells(obj_int)[grepl("GABA-[0-9]", obj_int$cluster_orig)]
  
  cells = union(cells, cells2)
  
  obj_int_filt = subset(obj_int, cells=cells)
  table1 = table(obj_int_filt$cluster_orig)
  
  plan(multiprocess(workers=10))
  obj_int_filt = SCTransform(obj_int_filt, vars.to.regress = c("percent.mito", "dataset"))
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
  
}

res_to_use = "cluster_orig1"

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
#obj_int_filt_avg = AverageExpression(obj_int_filt, assays = "SCT", slot="data")
#obj_int_filt_avg = obj_int_filt_avg[[1]]


# Marker ID --------------------------------------------------------------

#markers_fname = file.path(data_dir, sprintf("%s_gaba_sct_markers.rds", data_fname_prefix))
redo_markers = F
if (redo_markers) {
  DefaultAssay(obj_int_filt) = "SCT"
  Idents(obj_int_filt) = obj_int_filt$cluster_orig
  markers = FindAllMarkers(obj_int_filt, 
                           assay="SCT", 
                           test.use = "wilcox",
                           min.pct = .2, only.pos=T)
  saveRDS(markers, markers_fname)
} else {
  markers = readRDS(markers_fname)
}

# LGE only ----------------------------------------------------------------

markers_lge_fname = file.path(out_dir, sprintf("hvc_ra_x_gaba_lge_sct_markers.rds"))
cells = Cells(obj_int_filt)[obj_int_filt$cluster_orig1 %in% c("GABA-1", "0", "1", "2", "3", "4", "9")]
obj_int_filt_lge = subset(obj_int_filt, cells=cells)

obj_int_filt_lge$cluster_lge = case_when(obj_int_filt_lge$cluster_orig1 %in% c("0", "1", "2", "3", "4") ~ "MSN",
                                         obj_int_filt_lge$cluster_orig1 %in% c("9") ~ "PN",
                                         obj_int_filt_lge$cluster_orig1 %in% "GABA-1" ~ "GABA-1")
Idents(obj_int_filt_lge) = obj_int_filt_lge$cluster_lge
redo_markers = F
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

write.table(markers_lge, file.path(out_dir, sprintf("hvc_ra_x_gaba_lge_sct_markers.txt")), quote=F, col.names=T, sep="\t", row.names=F)


# Markers,paired---------------------------------------------
markers_lge_fname = file.path(out_dir, sprintf("hvc_ra_x_gaba_lge_paired_sct_markers.rds", data_fname_prefix))

pairs = combn(c("GABA-1", "PN", "MSN"), 2, simplify=F)
names(pairs) = map_chr(pairs, ~paste(.x, collapse="_"))
redo_markers = T
if (redo_markers) {
  DefaultAssay(obj_int_filt_lge) = "SCT"
  
  markers_lge1 = map(pairs, function(p) {
    tmp = FindMarkers(obj_int_filt_lge, 
                      ident.1 = p[1],
                      ident.2 = p[2],
                      assay="SCT", 
                      test.use = "wilcox", 
                      logfc.threshold = .25, 
                      min.pct = .2, only.pos=F)
    tmp  %>% rownames_to_column(var="gene")
  }) %>% set_names(names(pairs)) %>% bind_rows(.id="pair")
  saveRDS(markers_lge1, markers_lge_fname)
} else {
  markers_lge1 = readRDS(markers_lge_fname)
}


# Dotplot -----------------------------------------------------------------
n_plot = 10
to_plot = markers_lge %>%
  group_by(cluster) %>%
  top_n(n_plot, avg_logFC) %>%
  ungroup() %>% 
  distinct(gene)

to_plot = to_plot$gene
to_plot = c("FOXP1", "FOXP2", "PPP1R1B", "DRD1", "DRD2", to_plot)
th = theme(axis.title = element_text(size=6),
           axis.text.x = element_text(size=5, angle=90, face="italic", hjust=1, vjust=.5 ),
           axis.text.y = element_text(size=5, hjust=1),
           axis.line = element_line(size=.25),
           axis.ticks = element_line(size=.25),
           legend.text = element_text(size=5),
           legend.title = element_text(size=5),
           legend.key.size = unit(.1, "in"),
           legend.direction = "horizontal")

Idents(obj_int_filt_lge) = obj_int_filt_lge$cluster_lge
Idents(obj_int_filt_lge) = factor(Idents(obj_int_filt_lge), levels=c( "MSN", "PN",  "GABA-1"))
cols = brewer_pal(palette="Reds")(9)
gg = DotPlot(obj_int_filt_lge, cols=c("white", cols[9]), features=to_plot,  assay = "SCT", dot.min = 0, dot.scale = 1) + 
  th + 
  labs(x="", y="")

gg

gg_nl = gg +
  theme(legend.position="none")

save_plot(file.path(out_dir, "dotplot.pdf"), gg_nl, base_height=1, base_width=2.5)

leg = get_legend(gg)
save_plot(file.path(out_dir, "dotplot_legend.pdf"), leg, base_height=3, base_width=3)


# Dotplot, manual ---------------------------------------------------------

n_plot = 10
to_plot = markers_lge %>%
  group_by(cluster) %>%
  top_n(n_plot, avg_logFC) %>%
  ungroup() %>% 
  distinct(gene, cluster) %>%
  mutate(cluster = factor(cluster, levels=c("MSN", "PN", "GABA-1"))) %>% 
  arrange(cluster)

to_plot = to_plot$gene
to_plot = unique(c("FOXP1", "FOXP2", "MEIS2", to_plot))

res_to_use = "cluster_lge"
DefaultAssay(obj_int_filt_lge) = "SCT"
obj_int_g = FetchData(obj_int_filt_lge, c(to_plot, res_to_use))
obj_int_ngg = obj_int_g %>% 
  gather(key = "gene", value = "value", -one_of(res_to_use)) %>%
  dplyr::rename(celltype = cluster_lge) 

thresh = 0
rgb2hex  <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
fill_func = function(x) rgb2hex(colorRamp(brewer_pal(palette="Reds")(9))(x))


obj_int_ngg_stat = obj_int_ngg %>% group_by(celltype, gene) %>%
  summarize(
    perc_expr = 100* sum(value>thresh) / n(),
    value_mean = mean(value)
  ) %>%
  ungroup() 

obj_int_ngg_stat = obj_int_ngg_stat %>%
  mutate(gene = factor(gene, levels=to_plot),
         celltype = factor(celltype, levels=c("MSN", "PN", "GABA-1")))
value_max = max(obj_int_ngg_stat$value_mean)
#value_min = .1
value_min = 0
expr_min = 0
expr_max = max(obj_int_ngg_stat$perc_expr)
#expr_thresh = .05
expr_thresh = 0
expr_scale_factor = .1

obj_int_ngg_stat = obj_int_ngg_stat %>% 
  group_by(gene) %>% 
  mutate(
    value_mean_scale = (value_mean - value_min) / (value_max - value_min),
    value_mean_scale = ifelse(value_mean_scale<0, 0, value_mean_scale),
    value_mean_fill =  fill_func(value_mean_scale)
  ) %>%
  mutate(
    perc_expr_scale = (perc_expr - expr_min) / (expr_max - expr_min),
    perc_expr_scale = ifelse(perc_expr_scale<expr_thresh, NA, perc_expr_scale),
    perc_expr_scale = perc_expr_scale * expr_scale_factor
  ) 


th = theme(axis.title = element_text(size=6),
           axis.text.x = element_text(size=5, angle=90, face="italic", hjust=1, vjust=.5 ),
           axis.text.y = element_text(size=5, hjust=1),
           axis.line = element_line(size=.25),
           axis.ticks = element_line(size=.25),
           legend.text = element_text(size=5),
           legend.title = element_text(size=5),
           legend.key.size = unit(.1, "in"),
           legend.direction = "horizontal")

gg = ggplot(obj_int_ngg_stat, aes(gene, celltype, color=value_mean, size=perc_expr)) + 
  geom_point(shape=20) + 
  scale_size_continuous(range=c(.002, 2)) + 
  scale_color_gradient(low="white", high=brewer_pal(palette="Reds")(9)[9]) + 
  th + 
  labs(x="", y="")
gg

gg_nl = gg +
  theme(legend.position="none")

save_plot(file.path(out_dir, "dotplot_manual.pdf"), gg_nl, base_height=1, base_width=2.75)

leg = get_legend(gg)
save_plot(file.path(out_dir, "dotplot_manual_legend.pdf"), leg, base_height=3, base_width=3)



# Species split -----------------------------------------------------------

to_plot= c("FOXP1", "FOXP2", "MEIS2", "PPP1R1B", "DRD1", "DRD2")
Idents(obj_int_filt_lge) = paste(obj_int_filt_lge$species, obj_int_filt_lge$cluster_lge, sep="_")
th = theme(axis.text.x = element_text(size=5, angle=90, face="italic", hjust=1, vjust=.5 ),
           axis.text.y = element_text(size=5, hjust=1),
           axis.line = element_line(size=.25),
           axis.ticks = element_line(size=.25),
           legend.text = element_text(size=5),
           legend.title = element_text(size=5),
           legend.key.size = unit(.1, "in"),
           legend.direction = "horizontal")

Idents(obj_int_filt_lge) = factor(Idents(obj_int_filt_lge), levels=c( "ZF_MSN", "ZF_PN", "ZF_GABA-1", "BF_GABA-1"))
cols = brewer_pal(palette="Reds")(9)
gg = DotPlot(obj_int_filt_lge, cols=c("white", cols[9]), features=to_plot,  assay = "SCT", dot.min = 0) + 
  th + 
  labs(x="", y="")
gg

# Dotplot, manual, ag genes ---------------------------------------------------------
ag_genes = get_axon_guidance_genes()

to_plot = markers_lge %>% 
  filter(p_val_adj < 1E-2) %>%
  filter(gene %in% ag_genes$external_gene_name)
#to_plot = c("EPHA4", "EFNA5", "EFNA3", "ROBO1", "ROBO2", "CXCR4", "ACKR3", "CXCL12", "NRG1", "NRG3", "ERBB4")
to_plot = unique(c(to_plot$gene, c("ERBB4", "NRG1", "NRG3")))
res_to_use = "cluster_lge"
DefaultAssay(obj_int_filt_lge) = "SCT"
obj_int_g = FetchData(obj_int_filt_lge, c(to_plot, res_to_use))
obj_int_ngg = obj_int_g %>% 
  gather(key = "gene", value = "value", -one_of(res_to_use)) %>%
  dplyr::rename(celltype = cluster_lge) 

thresh = 0
rgb2hex  <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
fill_func = function(x) rgb2hex(colorRamp(brewer_pal(palette="Reds")(9))(x))


obj_int_ngg_stat = obj_int_ngg %>% group_by(celltype, gene) %>%
  summarize(
    perc_expr = 100* sum(value>thresh) / n(),
    value_mean = mean(value)
  ) %>%
  ungroup() 

obj_int_ngg_stat = obj_int_ngg_stat %>%
  mutate(gene = factor(gene, levels=to_plot),
         celltype = factor(celltype, levels=c("MSN", "PN", "GABA-1")))
value_max = max(obj_int_ngg_stat$value_mean)
value_min = .1

expr_min = 10
expr_max = max(obj_int_ngg_stat$perc_expr)
expr_thresh = .05
expr_scale_factor = .1

obj_int_ngg_stat = obj_int_ngg_stat %>% 
  group_by(gene) %>% 
  mutate(
    value_mean_scale = (value_mean - value_min) / (value_max - value_min),
    value_mean_scale = ifelse(value_mean_scale<0, 0, value_mean_scale),
    value_mean_fill =  fill_func(value_mean_scale)
  ) %>%
  mutate(
    perc_expr_scale = (perc_expr - expr_min) / (expr_max - expr_min),
    perc_expr_scale = ifelse(perc_expr_scale<expr_thresh, NA, perc_expr_scale),
    perc_expr_scale = perc_expr_scale * expr_scale_factor
  ) 


th = theme(axis.title = element_text(size=6),
           axis.text.x = element_text(size=5, angle=90, face="italic", hjust=1, vjust=.5 ),
           axis.text.y = element_text(size=5, hjust=1),
           axis.line = element_line(size=.25),
           axis.ticks = element_line(size=.25),
           legend.text = element_text(size=5),
           legend.title = element_text(size=5),
           legend.key.size = unit(.1, "in"),
           legend.direction = "horizontal")

gg = ggplot(obj_int_ngg_stat, aes(gene, celltype, color=value_mean, size=perc_expr)) + 
  geom_point(shape=20) + 
  scale_size_continuous(range=c(.002, 2)) + 
  scale_color_gradient(low="white", high=brewer_pal(palette="Reds")(9)[9]) + 
  th + 
  labs(x="", y="") + 
  coord_flip()
gg

gg_nl = gg +
  theme(legend.position="none")

save_plot(file.path(out_dir, "dotplot_ag_manual.pdf"), gg_nl, base_height=1, base_width=2.5)

leg = get_legend(gg)
save_plot(file.path(out_dir, "dotplot_ag_manual_legend.pdf"), leg, base_height=3, base_width=3)

