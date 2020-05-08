library(Seurat)
library(tidyverse)
library(reshape2)
library(qs)
library(future)

library(rtracklayer)
library(fgsea)
library(msigdbr)

library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(scales)
library(colorspace)

theme_set(theme_cowplot())
source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")
source("~/data2/rstudio/birds/utils/stats.R")
source("~/data2/rstudio/birds/utils/go.R")

tx = import("~/data2/assembly/lonStrDom1/ncbi/lonStrDom1.ncbi.apollo6/merge_inter_apollo_processed_mt_gffread.gtf")
tx_mt = tx[seqnames(tx)=="MT"]

# Parameters --------------------------------------------------------------

res_to_use = "position2_cluster_int_sub2"
hclust_method = "ward.D"
sub_dir = sprintf("%s_%s", res_to_use, hclust_method)

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "marker_gene")
script_name = "combined_int_sub2_gaba_hvc_va_ra_markers"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive = T)

tree_dir = file.path(dev_dir, "trees")
script_data_dir = file.path(tree_dir, "celltypes_hclust_gaba_int_sub2")

data_out_obj_fname = file.path(script_data_dir, "obj_integrated_subclustered_gaba.qs")
data_out_avg_fname = file.path(script_data_dir, "average_expr.qs")

# Load data ---------------------------------------------------------------

res_to_use = "cluster_int_sub2"
obj_int_filt = qread(data_out_obj_fname)
obj_int_avg1 = qread(data_out_avg_fname)

Idents(obj_int_filt) = FetchData(obj_int_filt, res_to_use)
obj_int_filt = obj_int_filt %>% SCTransform(vars.to.regress="percent.mito")

# Calculate DEG -----------------------------------------------------------------

assay_to_use = "SCT"
DefaultAssay(obj_int_filt) = assay_to_use


markers_fname = file.path(out_dir, sprintf("marker_genes_gaba_%s_%s_hvc_vs_ra_seurat.rds", res_to_use, assay_to_use))

cts = unique(FetchData(obj_int_filt, res_to_use)[,1])
cts = cts[! (cts %in% c("GABA-6", "GABA-7", "GABA-8", "GABA-Pre"))]

ns = table(FetchData(obj_int_filt, c("position2", "cluster_int_sub2")))[,cts]
max_cells = min(c(ns))

markers = map(cts, function(ct) {
  print(ct)
  cells = Cells(obj_int_filt)[obj_int_filt$cluster_int_sub2==ct]
  tmp = subset(obj_int_filt, cells=cells)
  print(tmp)
  Idents(tmp) = tmp$position2

  res = FindMarkers(tmp, 
                    test.use = "wilcox",
                    ident.1 = "hvc",
                    ident.2 = "ra", 
                    max.cells.per.ident = max_cells, min.pct = 0
  ) 
  
  
  res %>% rownames_to_column(var="gene")
}) %>% set_names(cts) %>% bind_rows(.id="cluster_int_sub2")

markers = markers %>% 
  mutate(p_val_adj_signed = -1 * log10(p_val_adj))

genes_to_remove = c("LOC115495643", "LOC115492526")
markers_filt = markers %>%
  filter(!(gene %in% genes_to_remove )) %>% 
  filter(!(gene %in% tx_mt$gene_name)) %>%
  filter(!grepl("AZ255", gene)) %>%
  mutate(gene = case_when(gene=="LOC100223925" ~ "CDH6",
                          TRUE~gene))
saveRDS(markers_filt, markers_fname)


# Write out flat file -----------------------------------------------------

write.table(markers_filt, file.path(out_dir, "hvc_vs_ra_deg.txt"), sep="\t", row.names=F, quote=F)

# Plot Volcano ------------------------------------------------------------

markers_top = markers_filt %>%
  mutate(sign = sign(avg_logFC)) %>%
  group_by(cluster_int_sub2, sign) %>% 
  top_n(5, abs(avg_logFC)) %>%
  mutate(to_label = gene) %>%
  select(to_label, gene)

markers1 = markers_filt %>% left_join(markers_top)

gg = ggplot(markers1, aes(avg_logFC, p_val_adj_signed)) + 
  geom_point() + 
  geom_label_repel(aes(label=to_label)) + 
  facet_wrap(~cluster_int_sub2)
gg


# Plot heatmap ------------------------------------------------------------


nsig = 5
markers_int_top = markers %>% 
 # filter(avg_logFC>0) %>% 
  mutate(sign = avg_logFC>0) %>%
  group_by(cluster_int_sub2, sign) %>%
  top_n(-1 * nsig, p_val_adj) %>%
  top_n(nsig, abs(avg_logFC)) %>%
  filter(!(gene %in% tx_mt$gene_name)) %>%
  filter(!grepl("AZ255", gene)) %>% 
  ungroup() %>% 
  distinct(gene, .keep_all=T)

genes = markers_int_top$gene

Idents(obj_int_filt) = "position2_cluster_int_sub2"
obj_int_avg = AverageExpression(obj_int_filt, assays=c("SCT", "integrated", "RNA"), slot="data")
obj_int_avg1 = log1p(obj_int_avg[["SCT"]])
obj_int_avg_filt = obj_int_avg1[genes,]


# Heatmap -----------------------------------------------------------------
ag_genes = c("ERBB4", "CXCL12", "CXCR4")
obj_int_avg_filt = obj_int_avg1[ag_genes,]
blues = RColorBrewer::brewer.pal(9, "Blues")
obj_int_avg_filt_mat = as.matrix(obj_int_avg_filt)
obj_int_avg_filt_mat[is.nan(obj_int_avg_filt_mat)] = min(c(obj_int_avg_filt_mat), na.rm=T)

obj_int_avg_filt_mat_scale = t(apply(obj_int_avg_filt_mat, 1, scale))
colnames(obj_int_avg_filt_mat_scale) = colnames(obj_int_avg_filt_mat)

cnames = colnames(obj_int_avg_filt_mat_scale)

cols = colorRamp2(c(min(c(obj_int_avg_filt_mat_scale)), 0, max(c(obj_int_avg_filt_mat_scale))),
                  c("blue", "white", "red")
                  )

obj_int_avg_filt_mat_scale = obj_int_avg_filt_mat_scale[,c(paste(c("hvc", "ra"), rep(paste("GABA", 1:8, sep="-"), each=2),  sep="_"))]
hm = Heatmap(obj_int_avg_filt_mat_scale,
             cluster_columns = F,
             clustering_method_rows = "ward.D",
             show_row_dend = F,
             show_column_names = T, 
             column_names_side = "top",
             col = cols
)
hm
pdf(file.path(out_dir, "heatmap.pdf"), height=12, width=7)
hm
dev.off()

# Heatmap, fold change -----------------------------------------------------------------

res_mat = acast(markers, gene~cluster_int_sub2, value.var="avg_logFC", fill=0)
res_mat_filt = res_mat[rownames(res_mat) %in% genes,]

cols = colorRamp2(c(min(c(res_mat_filt)), 0, max(c(res_mat_Filt))),
                  c("blue", "white", "red")
)

hm = Heatmap(res_mat_filt,
             cluster_columns = T,
             clustering_method_columns = "ward.D",
             clustering_method_rows = "ward.D",
             show_row_dend = F,
             show_column_names = T, 
             column_names_side = "top",
             col = cols
)
hm
pdf(file.path(out_dir, "heatmap_fc.pdf"), height=12, width=7)
hm
dev.off()

# Scatter -----------------------------------------------------------------

walk(cts, function(ct) {
markers1_cur = markers1 %>% filter(cluster_int_sub2==ct) %>%
  mutate(to_color = case_when(sign>0 & !is.na(to_label) ~ "HVC",
                           sign<0 & !is.na(to_label) ~ "RA",
                           is.na(to_label) ~ "NONE")) %>%
  select(gene, to_label, to_color)
obj_int_avg1_df = obj_int_avg1 %>% rownames_to_column(var = "gene") %>% melt() 
colnames(obj_int_avg1_df) = c("gene", "cluster", "value")

tmp = obj_int_avg1_df %>% filter(grepl(ct, cluster)) %>%
  filter(!(gene %in% genes_to_remove)) %>% 
  spread(key = cluster, value = value)
colnames(tmp) = c("gene", "hvc", "ra")
tmp = tmp %>% left_join(markers1_cur) %>%
  mutate(to_color = if_else(is.na(to_color), "NONE", to_color)) %>%
  mutate(to_label = if_else(is.na(to_label), "", to_label))
gg = ggplot(tmp, aes(ra, hvc)) + 
  geom_abline(slope=1, intercept=0, linetype=2) +
  geom_point(aes(color=to_color)) +
  geom_text_repel(aes(label=to_label)) + 
  labs(x="RA - log2 average expression",
       y="HVC - log2 average expression",
       title=ct) + 
  scale_color_manual(values=c(NONE="black", position_pretty_colors2[c("HVC", "RA")])) + 
  theme(legend.position="none") + 
  
  coord_equal()
gg
save_plot(file.path(out_dir, sprintf("scatter_%s.pdf", ct)), gg, base_width=4, base_asp = 1.3)
})


# FGSEA -------------------------------------------------------------------

markers_full_fname = file.path(out_dir, "markers_full.rds")
redo = F
if (redo) {
  markers_full = map(cts, function(ct) {
    print(ct)
    cells = Cells(obj_int_filt)[obj_int_filt$cluster_int_sub2==ct]
    tmp = subset(obj_int_filt, cells=cells)
    print(tmp)
    Idents(tmp) = tmp$position2
    res = FindMarkers(tmp, 
                      test.use = "wilcox",
                      ident.1 = "hvc",
                      ident.2 = "ra", 
                      max.cells.per.ident = max_cells, min.pct = 0, 
                      logfc.threshold = 0, 
                      min.diff.pct = 0
                      
    ) 
    
    
    res %>% rownames_to_column(var="gene")
  }) %>% set_names(cts) %>% bind_rows(.id="cluster_int_sub2")
  saveRDS(markers_full, markers_full_fname)
} else {
  markers_full = readRDS(markers_full_fname)
}
table <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")#GO, BP
len <- table(table$gs_name)
len_filt = len[len>20 & len<200]
table1 = table %>% filter(gs_name %in% names(len_filt))

pathways = split(table1$human_gene_symbol, table1$gs_name)

gsea_res =  markers_full %>% group_by(cluster_int_sub2) %>% do({
  rnks = .$p_val_adj
  names(rnks) = .$gene
  rnks = sort(rnks)
  fgsea_res = fgsea(pathways = pathways,
                    stats = rnks,
                    minSize=15,
                    maxSize=500,
                    nperm=10000)
  fgsea_res
})

cts = unique(markers_full$cluster_int_sub2)
gsea_res_ind = map(cts, function(ct) {
  markers_cur = markers_full %>% filter(cluster_int_sub2==ct)
  rnks = markers_cur$p_val_adj
  names(rnks) = markers_cur$gene
  rnks = sort(rnks)
  gsea_res_cur = gsea_res %>% filter(cluster_int_sub2==ct)
  fgsea_filt = gsea_res_cur %>% filter(padj<.3) %>% arrange(pval)
  collapsedPathways = collapsePathways(fgsea_filt, pathways = pathways, stats = rnks)
  gsea_res_cur %>% filter(pathway %in% collapsedPathways$mainPathways) %>%
    arrange(-NES)
}) %>% bind_rows()
gsea_res_ind_filt = gsea_res_ind %>% filter(padj<.1)
gsea_res_filt = gsea_res %>% filter(padj<.1)


# FGSEA, plot -------------------------------------------------------------

gsea_res_cur = gsea_res %>% filter(pathway %in% gsea_res_filt$pathway)
gsea_mat = acast(gsea_res_cur, pathway~cluster_int_sub2, value.var="NES", fill=0)

fname = file.path(out_dir, "fgsea_heat.pdf")

scale_factor = .1
hm_height = nrow(gsea_mat) * scale_factor
hm_width = ncol(gsea_mat) * scale_factor
pdf(fname, height=4, width=10)
hm = Heatmap(gsea_mat, name="NES",
             width=unit(hm_width, "in"),
             height=unit(hm_height, "in"),
             row_names_gp = gpar(fontsize=7),
             column_names_gp = gpar(fontsize=7))
draw(hm, heatmap_legend_side = "left")
dev.off()
