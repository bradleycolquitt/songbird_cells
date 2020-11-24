library(Seurat)
library(tidyverse)
library(qs)
library(qualpalr)
library(colorspace)
library(cowplot)
library(scales)
theme_set(theme_cowplot())

source("~/data2/rstudio/birds/utils/scRNA.R")
source("~/data2/rstudio/birds/utils/common_aesthetics.R")


# Directories -------------------------------------------------------------

dir_root = "~/sdd/data2/rstudio/birds/scRNA"
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dev_dir, "export")
data_fname= file.path(dev_data_dir, "HVC_RA_X.qs")

dir_out_root = "~/data2/rstudio/birds/scRNA"
dev_out_dir = file.path(dir_out_root, "devin_combined", "songbird_cells")
red_dir = file.path(dev_out_dir, "reduction_viz")
script_name = "combined_hvc_ra_x_umap_regress_zf_ortho"
red_dir = file.path(red_dir, script_name)
dir.create(red_dir, recursive = T)
figures_dir = file.path(red_dir)
dir.create(figures_dir)

data_out_obj_fname = file.path(red_dir, "obj_integrated_subclustered.qs")

res_to_use = "cluster_orig"

redo = T
dims = 30
n.neighbors = 50
min.dist = .5
reduction.name = sprintf("dims%snn%smindist%s", dims, n.neighbors, min.dist)

# Load data ---------------------------------------------------------------

if (redo) {
  

  
  obj_int_filt = qread(data_fname)
  #cells = colnames(obj_int_filt)[!is.na(obj_int_filt$position2)]
  #obj_int_filt = subset(obj_int_filt, cells = cells)
  
  x_key = data.frame(cluster_orig = as.character(seq(0,22)),
                     label = c("MSN (1)", "MSN (2)", "MSN (3)", "MSN (4)", "MSN (5)", "Astro (6)", "Astro (7)", "GLUT (8)", "Oligo (9)",
                               "PN (10)", "GLUT (11)", "Unk PPP1R1B (12)", "GABA SST (13)", "GABA PVALB (14)", "GABA PVALB (15)", "Micro (16)",
                               "GABA (17)", "Endo (18)", "GABA (19)", "OPC (20)", "GABA (21)", "GABA (22)", "GABA CHAT (23)"))
 # x_key_to_use = x_key %>% filter(grepl("MSN|GABA|PN", label))
  
  # UMAP --------------------------------------------------------------------
  
  DefaultAssay(obj_int_filt) = "integrated"
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


res_to_use = "label"
md = FetchData(obj_int_filt, "cluster_orig") %>%
  rownames_to_column() %>%
  left_join(x_key) %>%
  mutate(label = if_else(is.na(label), cluster_orig, as.character(label))) %>%
  column_to_rownames()

obj_int_filt = AddMetaData(obj_int_filt, md)
# Plot --------------------------------------------------------------------



cats = c("species", "region", res_to_use)
for ( ca in cats ) {
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=T, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none", panel.border = element_rect(fill=NA, color=1)) +
  labs(x="UMAP1", y="UMAP2") #+ 

gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", ca, reduction.name)), gg, base_height=7, base_asp =1 )
}

ca = "region"
pos_colors = position_pretty_colors2[c("HVC", "RA", "Area X")]
names(pos_colors) = c("HVC", "RA", "X")
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position="none", 
        panel.border = element_rect(fill=NA, color=1)) +
  scale_color_manual(values=pos_colors) + 
labs(x="", y="")
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=7, base_asp =1.1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=7, base_asp =1.1 )




# Colored by category -----------------------------------------------------
text_categories = c("Glutamatergic", "MSN", "GABAergic", "Neurogenic", "Astrocyte", "Oligodendrocyte", 
                    "OPC", "Microglia", "Endothelial", "Blood")
text_colors = qualitative_hcl(length(text_categories))
text_colors = brewer_pal(type="qual", palette=6)(length(text_categories))
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
                                     celltype_class == "MSN" ~ "MSN",
                                     grepl("Pre", celltype) ~ "Neurogenic",
                                     celltype_class %in% c("Astro", "Epen")~ "Astrocyte",
                                     celltype_class %in% c("Oligo") ~ "Oligodendrocyte",
                                     celltype_class %in% c("OPC") ~ "OPC",
                                     celltype_class %in% c("Endo", "Mural") ~ "Endothelial",
                                     celltype_class %in% c("Micro") ~ "Microglia",
                                     celltype_class %in% c("RBC") ~ "Blood")) %>%
  mutate(text_colors = text_colors[celltype_class1])

md = FetchData(obj_int_filt, "cluster_int_sub2") %>% 
  rownames_to_column() %>%
  left_join(ct_df)  %>%
  column_to_rownames()
obj_int_filt = AddMetaData(obj_int_filt, md)
em = Embeddings(obj_int_filt, reduction.name)

ct_df1 = ct_df %>% distinct(celltype_class1, .keep_all=T) %>%
  mutate(celltype_class1 = factor(celltype_class1, levels=text_categories)) %>%
  arrange(celltype_class1) %>% 
  mutate(
    xpos = min(em[,1]) - 0 * (max(em[,1]) - min(em[,1])),
    ypos = seq(max(em[,2]),by = -1 * .04 * (max(em[,2]) - min(em[,2])), length.out = n())
    )

ncells = length(Cells(obj_int_filt))
ca = "text_colors"
arrow_span = .1 * (max(em[,1]) - min(em[,1]))
arrow_obj = arrow(length = unit(0.03, "npc"), type="closed")
x_span = max(em[,1]) - min(em[,1])
y_span = max(em[,2]) - min(em[,2])
gg = DimPlot(obj_int_filt, reduction=reduction.name, group.by=ca, label=F, repel = T ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position="none"
        ) +
  scale_color_identity() +
  geom_text(
    data=ct_df1,
    aes(
      xpos, 
      ypos, 
      label=celltype_class1, 
      color=text_colors
    ),
    hjust=0,
    size=5
   ) + 
  annotate(x=max(em[,1]),
           y=min(em[,2]),
           label=sprintf("n=%s", format(ncells, big.mark = ",")),
           geom="text",
           hjust=1,
           size=5)+
coord_equal() 
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=6, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=6, base_asp =1 )

# Colored by species -----------------------------------------------------


md = FetchData(obj_int_filt, "species") %>% 
  rownames_to_column() %>%
  mutate(species_nice = case_when(species=="bf" ~ "Bengalese finch - cells",
                                  species=="zf" ~ "Zebra finch - nuclei")) %>% 
  mutate(species_nice = factor(species_nice, levels=c(  "Bengalese finch - cells", 
                                                        "Zebra finch - nuclei"
                                                    ))) %>% 
  column_to_rownames()
obj_int_filt = AddMetaData(obj_int_filt, md)
em = Embeddings(obj_int_filt, reduction.name)

ct_df1 = ct_df %>% distinct(celltype_class1, .keep_all=T) %>%
  mutate(celltype_class1 = factor(celltype_class1, levels=text_categories)) %>%
  arrange(celltype_class1) %>% 
  mutate(
    xpos = min(em[,1]) - .1 * (max(em[,1]) - min(em[,1])),
    ypos = seq(max(em[,2]),by = -1 * .04 * (max(em[,2]) - min(em[,2])), length.out = n())
  )

ncells = table(obj_int_filt$species)
names(ncells) = c("Bengalese", "Zebra")
ca = "species_nice"
arrow_span = .1 * (max(em[,1]) - min(em[,1]))
arrow_obj = arrow(length = unit(0.03, "npc"), type="closed")
x_span = max(em[,1]) - min(em[,1])
y_span = max(em[,2]) - min(em[,2])


gg = DimPlot(obj_int_filt,
             reduction=reduction.name, 
             group.by=ca,
             label=F
             ) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position=c(.6,.95)) + 
  annotate(x=max(em[,1]),
           y=min(em[,2]),
           label=sprintf("%s n=%s", names(ncells[1]), format(ncells[1], big.mark = ",")),
           geom="text",
           hjust=1,
           size=5)+
  annotate(x=max(em[,1]),
           y=min(em[,2] + y_span *.1),
           label=sprintf("%s n=%s", names(ncells[2]), format(ncells[2], big.mark = ",")),
           geom="text",
           hjust=1,
           size=5)
gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=6, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=6, base_asp =1 )


# Colored by species, split -----------------------------------------------------


md = FetchData(obj_int_filt, "species") %>% 
  rownames_to_column() %>%
  mutate(species_nice = case_when(species=="bf" ~ "Bengalese finch - cells",
                                  species=="zf" ~ "Zebra finch - nuclei")) %>% 
  mutate(species_nice = factor(species_nice, levels=c(  "Bengalese finch - cells", 
                                                        "Zebra finch - nuclei"
  ))) %>% 
  column_to_rownames()
obj_int_filt = AddMetaData(obj_int_filt, md)

ct_df1 = ct_df %>% distinct(celltype_class1, .keep_all=T) %>%
  mutate(celltype_class1 = factor(celltype_class1, levels=text_categories)) %>%
  arrange(celltype_class1) %>% 
  mutate(
    xpos = min(em[,1]) - .1 * (max(em[,1]) - min(em[,1])),
    ypos = seq(max(em[,2]),by = -1 * .04 * (max(em[,2]) - min(em[,2])), length.out = n())
  )

ncells = table(obj_int_filt$species)
names(ncells) = c("Bengalese", "Zebra")
ca = "species_nice"
arrow_span = .1 * (max(em[,1]) - min(em[,1]))
arrow_obj = arrow(length = unit(0.03, "npc"), type="closed")
x_span = max(em[,1]) - min(em[,1])
y_span = max(em[,2]) - min(em[,2])

sns = unique(obj_int_filt$species_nice)
species_colors = brewer_pal(type="qual")(2)
names(species_colors) = sns

em = Embeddings(obj_int_filt, reduction = reduction.name)
md = obj_int_filt@meta.data
em = cbind(em, md[match(rownames(em), rownames(md)),])

ggs = map(sns, function(sn) {
  
  colors_cur = species_colors
  colors_cur[names(colors_cur)!=sn] = "grey50"
  
  alphas = c(1, .1)
  names(alphas) = c(sn, setdiff(names(colors_cur), sn))
  
  gg = ggplot(em, aes_string("UMAP_1",
                             "UMAP_2")) + 
    geom_point(aes_string(color=ca, alpha=ca), size=.05) + 
    scale_color_manual(name="", values=colors_cur, guide=F) + 
    scale_alpha_manual(values=alphas, guide=F) + 
    labs(title=sn, subtitle = sprintf("n = %s", format(ncells[sn], big.mark = ","))) + 
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()
          ) 

  gg
})
pg = plot_grid(plotlist = ggs, align = "v", axis="lr", ncol=2)
pg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_split.pdf", reduction.name, ca)), pg, base_height=6, base_asp =1, ncol=2 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_split.png", reduction.name, ca)), pg, base_height=6, base_asp =1, ncol=2 )

pg_vert = plot_grid(plotlist = ggs, align = "v", axis="lr", nrow=2)
pg_vert
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_split_vert.pdf", reduction.name, ca)), pg_vert, base_height=6, base_asp =1, nrow=2 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s_split_vert.png", reduction.name, ca)), pg_vert, base_height=6, base_asp =1, nrow=2 )

# Colored by region -----------------------------------------------------


md = FetchData(obj_int_filt, "position2") %>% 
  rownames_to_column() %>%
  mutate(position_nice = toupper(position2))  %>% 
  column_to_rownames()
obj_int_filt = AddMetaData(obj_int_filt, md)

ncells = table(obj_int_filt$position2)
names(ncells) = toupper(names(ncells))
ca = "position_nice"
arrow_span = .1 * (max(em[,1]) - min(em[,1]))
arrow_obj = arrow(length = unit(0.03, "npc"), type="closed")
x_span = max(em[,1]) - min(em[,1])
y_span = max(em[,2]) - min(em[,2])

em = Embeddings(obj_int_filt, reduction = reduction.name)
md = obj_int_filt@meta.data
em = cbind(em, md[match(rownames(em), rownames(md)),])

pt.size = min(1583 / nrow(em), 1)
gg = ggplot(em, aes_string("UMAP_1",
                           "UMAP_2"))+
  geom_point(aes_string(color=ca), alpha=.2, size=pt.size) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position=c(.2,.95)
       ) +
  scale_color_manual(values=position_pretty_colors2[c("HVC", "RA")]) + 
annotate(x=max(em[,1]),
         y=min(em[,2]),
         label=sprintf("%s n=%s", names(ncells[1]), format(ncells[1], big.mark = ",")),
         geom="text",
         hjust=1,
         size=5)+
annotate(x=max(em[,1]),
         y=min(em[,2] + y_span *.1),
         label=sprintf("%s n=%s", names(ncells[2]), format(ncells[2], big.mark = ",")),
         geom="text",
         hjust=1,
         size=5)

gg
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.pdf", reduction.name, ca)), gg, base_height=6, base_asp =1 )
save_plot(file.path(figures_dir, sprintf("umap_%s_%s.png", reduction.name, ca)), gg, base_height=6, base_asp =1 )





