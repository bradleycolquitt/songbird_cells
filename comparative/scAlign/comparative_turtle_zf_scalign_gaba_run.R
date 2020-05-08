
library(Seurat)
library(tidyverse)
library(cowplot)
library(qs)
library(scAlign)

# Parameters -------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
perplexity = as.numeric(args[1])
batch.size = as.numeric(args[2])
architecture = args[3]
seed = as.numeric(args[4])

print(perplexity)
print(batch.size)

# Env variables -----------------------------------------------------------
brad_dir = Sys.getenv("BRAD_HOME")
sge_gpu = Sys.getenv("SGE_GPU")
print(sge_gpu)
if (sge_gpu == "") {
   sge_gpu = "0"
   device = "CPU"
} else {
  device = "GPU"
}

device = "CPU"
sge_gpu = ""
dir_root = file.path(brad_dir, "data2/rstudio/birds/scRNA")
dev_dir = file.path(dir_root, "devin_combined", "finch_cells")
dev_data_dir = file.path(dir_root, "devin_combined", "data")
dev_out = file.path(dev_dir, "preprocessing", "integrate", "zf_bf", "joint2", "SCT", "song")
dev_out_sub_dir = file.path(dev_out, sprintf("anchor%s_filter%s_score%s_maxfeatures%s_dims%s", 5, 200, 30, 200, 30))

out_dir = file.path(dev_dir, "comparative", "integrated", "gaba")
script_name = "comparative_turtle_zf_scalign_gaba_prep"
out_dir = file.path(out_dir, script_name)
dir.create(out_dir, recursive=T)

out_sub_dir = file.path(out_dir, sprintf("perplexity%s_batchsize%s_arch%s_%s", perplexity, batch.size, architecture, seed))
dir.create(out_sub_dir)

# Load data ---------------------------------------------------------------
cat("Loading data...\n")
fname_int = file.path(out_dir, "obj_scalign.qs")
fname_int_out = file.path(out_sub_dir, "obj_scalign.qs")
sca = qread(fname_int)

redo = T
if (redo) {

cat("Running scAlign...\n")
sca = scAlign(sca,
                   options=scAlignOptions(steps=15000,
                                          log.every=5000,
                                          norm=TRUE,
					  perplexity=perplexity,
                                          batch.size=batch.size,
                                          dropout.layer=T,
                                          batch.norm.layer=TRUE,
                                          early.stop=TRUE,
                                          architecture = "small",
					  gpu.device=sge_gpu,
					  seed=seed),
                   encoder.data="CCA",
                   decoder.data = "logcounts",
                   supervised='none',
                   run.encoder=TRUE,
                   run.decoder=FALSE,
                   log.dir=out_dir,
                   device=device)
  qsave(sca, fname_int_out)
} else {
  sca = qread(fname_int_out)
}

cn = colnames(sca)
cn[duplicated(cn)] = paste(cn[duplicated(cn)], 1:sum(duplicated(cn)), sep="_")
colnames(sca) = cn
sca_seu = as.Seurat(sca, counts = "counts", data = "logcounts")

reductions = c("ALIGNED.CCA", "CCA")
for (red in reductions) {

  sca_seu = RunUMAP(sca_seu,
                    reduction=red, 
                    dims=1:30, 
                    reduction.name=sprintf("umap%s", tolower(sub("\\.", "", red))), 
                    reduction.key = sprintf("umap%s_", tolower(sub("\\.", "", red))))
}

qsave(sca_seu, file.path(out_sub_dir, "obj_seurat.qs"))

cats = c("species", "dataset", "celltype_comb")
walk(cats, function(ca) {
  gg = DimPlot(sca_seu, reduction="umapalignedcca", group.by=ca, label=T) + 
theme(legend.position="none")
  save_plot(file.path(out_sub_dir, sprintf("aligned_%s.png", ca)), gg)
  #print(gg)
})

walk(cats, function(ca) {
  gg = DimPlot(sca_seu, reduction="umapcca", group.by=ca)  + 
theme(legend.position="none")
  save_plot(file.path(out_sub_dir, sprintf("cca_%s.png", ca)), gg)
  #print(gg)
})


