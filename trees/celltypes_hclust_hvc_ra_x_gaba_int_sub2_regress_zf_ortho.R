library(Seurat)
library(tidyverse)
library(qs)
library(pvclust)
library(future)
library(dendextend)
library(ggdendro)
library(ggraph)
library(circlize)
library(cowplot)
plan(sequential)

source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")

# Parameters ---------------------------------------------------------------


hclust_method = "average"
nsig = 50
sub_dir = sprintf("%s_%s", hclust_method, nsig)

# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_X.qs")

## Output dir
dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
tree_dir = file.path(dev_out_dir, "trees")
script_name = "celltypes_hclust_hvc_ra_x_gaba_int_sub2_regress_zf_ortho"
tree_dir = file.path(tree_dir, script_name)
tree_sub_dir = file.path(tree_dir, sub_dir)
dir.create(tree_sub_dir, recursive = T)

data_out_obj_fname = file.path(tree_dir, "obj_integrated_subclustered_gaba.qs")
data_out_avg_fname = file.path(tree_dir, "average_expr.qs")
data_out_pv_fname = file.path(tree_sub_dir, "pvclust.qs")



# Load finch data ---------------------------------------------------------

x_key = data.frame(cluster_orig = seq(0,22),
                   label = c("MSN (1)", "MSN (2)", "MSN (3)", "MSN (4)", "MSN (5)",
                             "Astro (6)", "Astro (7)",
                             "GLUT (8)", 
                             "Oligo (9)",
                             "PN (10)", "GLUT (11)", "Unk PPP1R1B (12)", "GABA SST (13)", "GABA PVALB (14)", "GABA PVALB (15)", "Micro (16)",
                             "GABA (17)", "Endo (18)", "GABA (19)", "OPC (20)", "GABA (21)", "GABA (22)", "GABA CHAT (23)"))
x_key_to_use = x_key %>% filter(grepl("MSN|GABA|PN", label))

redo = F
if (redo) {
  #load(data_fname)
  #obj_int = HVCX_seurat_REGION_res1.1
  obj_int = qread(data_fname)
  cells = Cells(obj_int)[obj_int$cluster_orig %in% x_key_to_use$cluster_orig]
  cells2 = Cells(obj_int)[grepl("GABA-[0-9]", obj_int$cluster_orig)]
  
  cells = union(cells, cells2)
  
  obj_int_filt = subset(obj_int, cells=cells)
  table(obj_int_filt$cluster_orig)
  
  plan(multiprocess(workers=10))
  obj_int_filt = SCTransform(obj_int_filt, vars.to.regress = c("percent.mito", "dataset"))
  
  qsave(obj_int_filt, data_out_obj_fname)
} else {
  obj_int_filt = qread(data_out_obj_fname)
  
}


# Average data ------------------------------------------------------------


#obj_int_filt$region_cluster_orig = paste(obj_int_filt$region, obj_int_filt$cluster_orig, sep="_")

md = FetchData(obj_int_filt, c("cluster_orig")) %>%
  distinct(cluster_orig)
md = md %>% mutate(s = map(cluster_orig, function(x) unlist(str_split(x, "-")))) %>%
  mutate(s1 = map_chr(s, function(x) {
    case_when(length(x) == 1 ~ paste(x, collapse = "-"),
              length(x) == 2 ~ paste(x, collapse = "-"),
              length(x) == 3 ~ paste(x[1:2], collapse='-'))
    })) %>%
  mutate(cluster_orig1 = s1)

md_full = FetchData(obj_int_filt, c("cluster_orig")) %>%
  rownames_to_column()
md_full = md_full %>% left_join(md) %>%
  select(rowname, cluster_orig, cluster_orig1) %>%
  column_to_rownames()
obj_int_filt = AddMetaData(obj_int_filt, md_full)
obj_int_filt$region_cluster_orig = paste(obj_int_filt$region, obj_int_filt$cluster_orig1, sep="_")
qsave(obj_int_filt, data_out_obj_fname)

Idents(obj_int_filt) = FetchData(obj_int_filt, "region_cluster_orig")
obj_int_filt_avg = AverageExpression(obj_int_filt, assays = c("SCT", "RNA"), slot="counts")
obj_int_filt_avg1 = log1p(obj_int_filt_avg[["SCT"]])


# Marker ID --------------------------------------------------------------

markers_fname = file.path(tree_dir, sprintf("hvc_ra_x_gaba_rna_markers.rds"))
redo_markers = T
if (redo_markers) {
  Idents(obj_int_filt) = obj_int_filt$cluster_orig1
  markers = FindAllMarkers(obj_int_filt, 
                           assay="RNA", 
                           test.use = "wilcox",
                           min.pct = .2, only.pos=T)
  saveRDS(markers, markers_fname)
} else {
  markers = readRDS(markers_fname)
}

write.table(markers, file.path(tree_dir, "hvc_ra_x_gaba_rna_markers.txt"), quote=F, col.names=T, sep="\t")
# Cluster -----------------------------------------------------------------


nsig = 50
markers_int_top = markers %>% 
  mutate(sign = avg_logFC>0) %>%
  #group_by(cluster, sign) %>%
  group_by(cluster) %>% 
  #top_n(-1 * nsig, p_val_adj) %>%
  #top_n(nsig, abs(avg_logFC)) %>%
  top_n(nsig, avg_logFC) %>%
  distinct(gene, .keep_all=T)

obj_int_avg_filt = obj_int_filt_avg1[markers_int_top$gene,]


pv = pvclust(obj_int_avg_filt, method.hclust=hclust_method, nboot = 100, parallel=4L)
plot(pv)
qsave(pv, data_out_pv_fname)


# dendextend ------------------------------------------------------------------

dend = pv$hclust %>% as.dendrogram %>%
  set("branches_k_color", k=5) #%>% 
  #set("branches_lwd", c(1.5,1,1.5)) %>%
  #set("branches_lty", c(1,1,3,1,1,2)) %>%
  #set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
  #set("nodes_pch", 19) %>% set("nodes_col", c("orange", "black", "plum", NA))
# plot the dend in usual "base" plotting engine:
plot(dend)


# circlize ----------------------------------------------------------------


#plot(pv)
dend = as.dendrogram(pv)
dend = dend %>%
  dendextend::set("branches_k_color", value=c( "#2c9f58", "black"), k = 2) 

dend_gg = as.ggdend(dend)
dend_seg = dend_gg$segments %>%
  mutate(col = if_else(is.na(col), "black", col))

cols = unique(dend_seg$col)

pop_factor = 2
dend_seg = dend_seg %>%
  mutate(is_root = y==max(y)) %>% 
  mutate(y = if_else(col==cols[1] & !is_root,  y - pop_factor, y),
         yend = if_else(col==cols[1],  yend - pop_factor, yend))

dend_seg = dend_seg[-1,]

md = FetchData(obj_int_filt, c("region", "cluster_orig1", "region_cluster_orig"))  %>%
  distinct(region, cluster_orig1, .keep_all=T) %>% 
  mutate(region = sub("X", "Area X", region)) %>% 
  mutate(region_color = position_pretty_colors2[region])


dend_leaves <- dend_gg$labels %>%
  rename(region_cluster_orig = label) %>%
  left_join(md) %>%
  left_join(x_key %>% mutate(cluster_orig1 = as.character(cluster_orig))) %>%
  mutate(label = ifelse(is.na(label), cluster_orig1, as.character(label))) %>%
  # filter(region!="RA") %>%
  mutate(label1 = if_else(grepl("MSN|PN|GABA-1", label), label, ""))
gg = ggplot(data = dend_seg, aes(x=x, y=y, xend=xend, yend=yend, color=col)) + 
  geom_segment(size=.25, linejoin="mitre") + 
  geom_text(data=dend_leaves, aes(x=x, y=min(dend_seg$y)-1, label=label1, color=region_color), inherit.aes = F, size=6/2.8) + 
  scale_y_reverse() + 
  coord_polar(theta = "x") + 
  scale_color_identity() + 
  theme_void()
gg
save_plot(file.path(tree_sub_dir, "hclust_circ.pdf"), gg, base_height=2, base_width=2)

