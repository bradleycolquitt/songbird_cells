library(Seurat)
library(tidyverse)
library(qs)

library(future)
library(lisi)
library(qualpalr)

source("~/data2/rstudio/birds/utils/scRNA.R")

# Directories -------------------------------------------------------------


dir_root = "~/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT_regress", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

compare_dir = file.path(dev_dir, "dataset_comparison")
script_name = "hvc_ra_lisi"
compare_dir = file.path(compare_dir, script_name)
dir.create(compare_dir, recursive = T)

data_fname= file.path(dev_out_sub_dir, "obj_integrated_subclustered.qs")


# Load data ---------------------------------------------------------------

obj_int_filt = qread(data_fname)


# LISI --------------------------------------------------------------------

text_categories = c("Glutamatergic", "GABAergic", "Neurogenic", "Astrocyte", "Oligodendrocyte", 
                    "OPC", "Microglia", "Endothelial", "Blood")

text_palette = qualpal(n = length(text_categories), 
                       colorspace = list(h = c(0, 360), s = c(.6,.8), l = c(0.2, 0.6))
)
plot(text_palette)
text_colors = text_palette$hex
names(text_colors) = text_categories

ct_df = obj_int_filt@meta.data %>% distinct(cluster_int_sub2) %>%
  mutate(celltype = as.character(cluster_int_sub2)) %>%
  
  mutate(celltype_class = map(celltype, ~unlist(strsplit(.x, "-"))[[1]])) %>%
  mutate(celltype_class1 = case_when(grepl("Glut", celltype_class) & !grepl("Pre", celltype)  ~ "Glutamatergic",
                                     celltype_class == "GABA" & !grepl("Pre", celltype)  ~ "GABAergic",
                                     grepl("Pre", celltype) ~ "Neurogenic",
                                     celltype_class %in% c("Astro", "Epen")~ "Astrocyte",
                                     celltype_class %in% c("Oligo") ~ "Oligodendrocyte",
                                     celltype_class %in% c("OPC") ~ "OPC",
                                     celltype_class %in% c("Endo", "Mural", "VLMC") ~ "Endothelial",
                                     celltype_class %in% c("Micro") ~ "Microglia",
                                     celltype_class %in% c("RBC") ~ "Blood")) %>%
  mutate(text_colors = text_colors[celltype_class1])

md = FetchData(obj_int_filt, c("cluster_int_sub2", "position2")) %>% 
  rownames_to_column() %>%
  left_join(ct_df)%>%
  column_to_rownames()

obj_int_filt = AddMetaData(obj_int_filt, md)

# Lisi --------------------------------------------------------------------

groups = na.omit(unique(md$celltype_class1))

res = map(groups, function(gr) {
  print(gr)
  
  cells = Cells(obj_int_filt)[obj_int_filt$celltype_class1==gr]
  tmp = subset(obj_int_filt, cells=cells)
  md_tmp = FetchData(tmp, c("position2"))
  dat = Embeddings(tmp, reduction="umap")
  compute_lisi(dat, md_tmp, c("position2"))
})

names(res) = groups
res1 = bind_rows(res, .id="celltype_class1") %>%
  left_join(ct_df) %>%
  mutate(celltype_class1 = factor(celltype_class1, levels=text_categories))


# Plot --------------------------------------------------------------------

gg = ggplot(res1, aes(position2, fill=celltype_class1, color=celltype_class1)) + 
  geom_density() + 
  facet_grid(celltype_class1~., scales="free") + 
  scale_color_manual(name = "", values = text_colors) +
  scale_fill_manual(name = "", values = text_colors) +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0, hjust=0)) + 
# xlim(.8,2.2) + 
  scale_x_continuous(breaks=c(1,2), limits=c(.8,2.2)) + 
  scale_y_continuous(breaks=pretty_breaks(n=2)) + 
  labs(x="Region LISI")
gg

save_plot(file.path(compare_dir, "region_lisi.pdf"), gg, base_height=.5, base_width=4, nrow=length(groups))

