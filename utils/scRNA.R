#library(furrr)
library(future)
#library(geiger)
library(qs)

#exc_levels = c("Exc_NECTIN3", "Exc_BDNF", "Exc_C1QTNF4", "Exc_GRIA4", "Exc_IEG", "Exc_TAC1", "Exc_CDH9", "Exc_PLCB1", "Exc_CALM2", "Exc_NR4A2", "Exc_NPTXR", "Exc_NDNF")
#exc_levels = c("Exc_CDH9", "Exc_PLCB1", "Exc_MAT2B", "Exc_NECTIN3", "Exc_BDNF", "Exc_KCNK1", "Exc_IEG", "Exc_TAC1", "Exc_NPTN", "Exc_NEUROD6", "Exc_CALM2",  "Exc_PVALB", "Exc_VAMP2",  "Exc_NR4A2", "Exc_NPTXR", "Exc_ITM2B", "Exc_REST")

RunUMAP_naive = function(object, reduction.use="pca", dims.use=NULL, reduction.key="UMAP", reduction.name="umap",
                         n_neighbors = 30, n_components = 2, min_dist = 0.3, metric = "pearson") {
  require(umap)
  cells.use = colnames(object@data)
  
  dim.code = GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes = paste0(dim.code, dims.use)
  data.use = GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "cell.embeddings")
  data.use = data.use[cells.use, dim.codes, drop = FALSE]
  
  umap_config = umap.defaults
  umap_config$n_neighbors = n_neighbors
  umap_config$n_components = n_components
  umap_config$metric = metric
  umap_config$min_dist = min_dist
  umap_config$random_state = 42
  umap_config$transform_state = 42
  umap_out = umap(data.use, config=umap_config, method="naive")
  umap_output = umap_out$layout
  
  colnames(umap_output) = paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(umap_output) = cells.use
  object = SetDimReduction(object = object, reduction.type = reduction.name, 
                            slot = "cell.embeddings", new.data = as.matrix(x = umap_output))
  object = SetDimReduction(object = object, reduction.type = reduction.name, 
                            slot = "key", new.data = reduction.key)
  return(object)
}

RunUMAP_naive_v3 = function(object, reduction.use="pca", dims.use=NULL, reduction.key="UMAP", reduction.name="umap",
                         n_neighbors = 30, n_components = 2, min_dist = 0.3, metric = "pearson") {
  require(umap)
  
  cells.use = colnames(object)
  
  #dim.code = Get
  #dim.code = GetDimReduction(object = object, reduction.type = reduction.use, 
                             # slot = "key")
  dim.codes = paste0(reduction.key, dims.use)
  data.use = Embeddings(object, reduction = reduction.use)
  #data.use = GetDimReduction(object = object, reduction.type = reduction.use, 
   #                           slot = "cell.embeddings")
  data.use = data.use[cells.use, dims.use, drop = FALSE]
  
  umap_config = umap.defaults
  umap_config$n_neighbors = n_neighbors
  umap_config$n_components = n_components
  umap_config$metric = metric
  umap_config$min_dist = min_dist
  umap_config$random_state = 42
  umap_config$transform_state = 42
  umap_out = umap(data.use, config=umap_config, method="naive")
  umap_output = umap_out$layout
  
  colnames(umap_output) = paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(umap_output) = cells.use
  object@reductions[["umap"]] = CreateDimReducObject(embeddings = as.matrix(x=umap_output), assay = "RNA", key="UMAP_")
  #object = SetDimReduction(object = object, reduction.type = reduction.name, 
  #                          slot = "cell.embeddings", new.data = as.matrix(x = umap_output))
  #object = SetDimReduction(object = object, reduction.type = reduction.name, 
  #                          slot = "key", new.data = reduction.key)
  return(object)
}

createCleanedSeurat_v3 = function(scl, scale.factor = 10000, ...) {
  if (!is(scl, "SoupChannelList")) 
    stop("scl must be a SoupChannelList object")
  srat = CreateSeuratObject(scl$toc)
  cleaned = log(1 + scl$strainedExp * scale.factor)
  rownames(cleaned) = gsub("_", "-", rownames(cleaned))
  cleaned = as.matrix(cleaned)
  params = list(object = srat, assay.type = "RNA", normalization.method = "LogNormalize", 
                scale.factor = scale.factor, display.progress = TRUE)
  #srat@calc.params[["NormalizeData"]] = params
  #srat@calc.params[["NormalizeData"]]$object = NULL
  #srat@calc.params[["NormalizeData"]]$object2 = NULL
  #srat@calc.params[["NormalizeData"]]$time = Sys.time()
  srat = SetAssayData(object = srat, assay = "RNA", slot = "data", 
                      new.data = cleaned)
  srat@meta.data$nGene = colSums(srat[[DefaultAssay(srat)]]@data > 0)
  return(srat)
}

plot_umap_genes = function(obj, genes, assay="alra") {
  DefaultAssay(obj) = assay 
  md = obj@meta.data
  #md = md %>% rownames_to_column()
  obj_umap = Embeddings(obj, reduction = "umap")
  obj_umap = obj_umap %>% as.data.frame() %>% rownames_to_column(var="rowname") %>% 
    left_join(md) 
  expr_tmp = FetchData(obj, vars=genes)
  expr_tmp = as.data.frame(expr_tmp) %>% rownames_to_column(var="rowname") %>% left_join(obj_umap)

  colnames(expr_tmp) = make.names(colnames(expr_tmp))
  colnames(expr_tmp) = sub("RNA_", "", colnames(expr_tmp))
    
  ggs = map(make.names(genes), function(gene1) {
    gg = ggplot(expr_tmp, aes_string("UMAP_1", "UMAP_2", color=gene1))
    gg = gg + geom_point()
    gg = gg + scale_color_viridis_c(option="magma")
    gg = gg + theme(axis.text=element_blank(),
                    axis.line=element_blank(),
                    axis.ticks = element_blank(),
                    axis.title = element_blank())
    
    return(gg)
  })
  ggs

}

## Clustering 

de_pairs = function(obj, de_score_thresh=150, n_neighbors=2, pca_dist=NULL, to_return = "ids", mt_genes=NULL) {
  
  pairs = NULL
  if (is.null(pca_dist)) {
    ## calculate DE between tip pairs
    
    tree = Tool(obj, slot = "BuildClusterTree")
    edges = tree$edge
    tip_edges = edges[edges[,2] %in% as.numeric(tree$tip.label),]
    pairs = split(tip_edges[,2], tip_edges[,1])
    pairs = pairs[map_lgl(pairs, ~length(.x)==2)]
  } else {
    if (ncol(pca_dist) >= (1+n_neighbors)) {
      pairs_compact = lapply(rownames(pca_dist), function(i) {
        x_sorted = sort(pca_dist[i,])
        x_names = names(x_sorted)[2:(1+n_neighbors)]
        map(x_names, ~ c(i,.x))
      })
    pairs = unlist(pairs_compact, recursive = F)
    names(pairs) = 1:length(pairs)
    } else{
      pairs = list(`1`=c("1", "2"))
    }
  }
  
  marker_pairs = map(pairs, function(x)  {
    print(x)
    res = FindMarkers(obj, ident.1=x[1], ident.2=x[2], verbose=F)
    res %>% rownames_to_column(var="gene")
  }) %>% set_names(names(pairs)) %>% bind_rows(.id="pair")
  
  tips1 = map_chr(pairs, 1)
  tips2 = map_chr(pairs, 2)
  pairs_df = data.frame(pair = names(pairs), tip1 = tips1, tip2 = tips2)
  
  marker_pairs = marker_pairs %>% left_join(pairs_df)
  if (to_return == "markers") {
    return(marker_pairs)
  } else if (to_return == "ids") {
    
    de_score= function(x) {
      x_trans =  -1 * log10(x)
      x_trans =  if_else(x_trans>=20, 20, x_trans)
      sum(x_trans)
    }
    
    if (!is.null(mt_genes)) {
      marker_pairs = marker_pairs %>% filter(!(gene %in% mt_genes))

    }
    markers_pairs_stats = marker_pairs %>% group_by(pair) %>% summarize(n_sig = sum(p_val_adj<.05),
                                                                        de_score = de_score(p_val_adj))
    
    pairs_to_merge = markers_pairs_stats %>% filter(de_score < de_score_thresh)
    pairs_to_merge_list = pairs[pairs_to_merge$pair]
    ids = as.character(Idents(obj))
    
    if (length(pairs_to_merge_list) > 0 ) {
      pairs_to_merge_list = map(pairs_to_merge_list, ~sort(as.numeric(.x)))
      pairs_to_merge_list = rev(pairs_to_merge_list)
      for(pair in pairs_to_merge_list) {
        ids[ids==pair[2]] = pair[1]
      }
    }
    return(ids)
  }
}

relabel_clusters = function(obj) {
  
  cts_tb = sort(table(Idents(obj)), decreasing=T)
  new_ids_order = 1:length(cts_tb)
  Idents(obj) = plyr::mapvalues(x = Idents(obj), from = names(cts_tb), to = new_ids_order)
  Idents(obj) = factor(Idents(obj), levels=sort(as.numeric(unique(Idents(obj)))))
  return(obj)
}

calc_pca_dist = function(obj, dims) {
  my.lapply = lapply
  embeddings = Embeddings(object = obj, reduction = "pca")[, 
                                                            dims]
  data.dims = my.lapply(X = levels(x = obj), FUN = function(x) {
    cells = WhichCells(object = obj, idents = x)
    if (length(x = cells) == 1) {
      cells = c(cells, cells)
    }
    temp = colMeans(x = embeddings[cells, ])
  })
  data.dims = do.call(what = "cbind", args = data.dims)
  colnames(x = data.dims) = levels(x = obj)
  data.dist = dist(x = t(x = data.dims))
  return(as.matrix(data.dist))
}



combine_gene_fragments = function(obj_orig, features_to_transfer=NULL ) {
  require(Matrix.utils)
  dat = GetAssayData(obj_orig, assay = "RNA", slot="counts")
  rnames = rownames(dat)
  groups = sub("\\.[0-9]+$", "", rnames)
  groups = sub("\\-[0-9]+$", "", groups)
  dat1 = aggregate.Matrix(dat, groups, fun="sum" )
  rnames = rownames(dat1)
  rnames_to_include = which(rnames!="NA")
  dat1 = dat1[rnames_to_include,]
  
  obj = CreateSeuratObject(dat1)
  
  
  if (!is.null(features_to_transfer)) {
    features_to_transfer = c(features_to_transfer, "percent.mito")
  } else {
    features_to_transfer = colnames(obj_orig@meta.data)

  }
  obj = AddMetaData(obj, FetchData(obj_orig, features_to_transfer))
  return(obj)
}


merge_idents = function(obj, cluster_column, tree_dims=1:20, de_thresh, cluster_size_min = 30, comparison="neighbors", n_neighbors=2, niter=10, mt_genes=NULL) {
  DefaultAssay(obj) = "SCT"
  Idents(obj) = FetchData(obj, cluster_column)
  cts_tb = table(Idents(obj))
  
  #obj@reductions$umap@assay.used
  obj = BuildClusterTree(obj, dims=tree_dims)
  tree = Tool(obj, slot = "BuildClusterTree")
  tree$node.label = paste("node", 1:tree$Nnode, sep="_")
  plot(tree, show.node.label=T)
  
  ## Merge low cell clusters
  ct_to_merge = names(cts_tb)[cts_tb<=cluster_size_min]
  
  if (length(ct_to_merge) > 0) {
  edges = tree$edge
  node_to_merge = edges[edges[,2]==ct_to_merge,1]
  tips_to_merge = tips(tree, node_to_merge)
  tips_to_merge = tips_to_merge[tips_to_merge!=ct_to_merge]
  
  
  if (any(tips_to_merge > length(cts_tb))) {
    
    node_id = tips_to_merge[tips_to_merge > length(cts_tb)]
    sister_tips = edges[edges[,1]==node_id,2]
    ident_dists = as.matrix(calc_pca_dist(obj, dims.use))
    sample_probs = ident_dists[ct_to_merge,tips_to_merge]
    sample_probs = sample_probs/sum(sample_probs)
    dist_tip = sample(tips_to_merge, length(Idents(obj)[Idents(obj)==tips_to_merge]), replace=T, prob = sample_probs )
  } else {
    dest_tip = tips_to_merge[1]
  }
  Idents(obj)[Idents(obj)==ct_to_merge] = dest_tip
  
  obj = relabel_clusters(obj)
  }
  new_ids = Idents(obj)
  
  i = 1
  start = TRUE
  print("Starting ident merge...")
  while (start || (length(unique(new_ids)) < length(unique(Idents(obj))) && i < niter)) {
    start = FALSE
    print(i)
    Idents(obj) = new_ids
    obj = relabel_clusters(obj)
    obj = BuildClusterTree(obj, dims=tree_dims)
    tree = Tool(obj, slot = "BuildClusterTree")
    tree$node.label = paste("node", 1:tree$Nnode, sep="_")
    plot(tree, show.node.label=T)
    
    pca_dist = NULL
    if (comparison == "neighbors") {
      pca_dist = calc_pca_dist(obj, dims=tree_dims)
    }
    print("... Calculating DE ...")
    new_ids = de_pairs(obj, de_thresh,  n_neighbors=n_neighbors, pca_dist=pca_dist, mt_genes = mt_genes)
    i = i + 1
    print(i)
  }
  
  Idents(obj) = new_ids
  obj = relabel_clusters(obj)
  
  return(obj)
}

extract_dims_use = function(obj) {
  pvalues = obj@reductions$pca@jackstraw$overall.p.values
  dims.use = pvalues[pvalues[,2]<.01,1]
  dims.use
}

sctransform_umap = function(obj, vars_to_regress = "percent.mito", assay="RNA") {
  DefaultAssay(obj) = assay
  VariableFeatures(obj) = NULL
  obj = SCTransform(obj,
                           do.correct.umi = T,
                           return.only.var.genes = F,
                           vars.to.regress = vars_to_regress) %>%
    RunPCA()
  obj = JackStraw(object = obj, num.replicate = 100, reduction = "pca", dims = 20 )
  obj = ScoreJackStraw(obj, dims = 1:20)
  JackStrawPlot(object = obj, dims = 1:20)
  pvalues = obj@reductions$pca@jackstraw$overall.p.values
  dims.use = pvalues[pvalues[,2]<.01,1]
  obj = RunUMAP(obj,  reduction = "pca", dims = dims.use)
  # DefaultAssay(obj) = "SCT"
  # slot(obj, 'tools')[["RunALRA"]] = NULL
  # obj = RunALRA(obj)
  obj
}

seurat_jackstraw = function(obj) {
  obj = JackStraw(object = obj, num.replicate = 100, reduction = "pca", dims = 30 )
  obj = ScoreJackStraw(obj, dims = 1:30)
  obj
}


seurat_find_clusters = function(obj, ress, cluster_prefix, dims.use=NULL, algorithm=4, reduction ="pca") {

  if (is.null(dims.use)) {
    dims.use = extract_dims_use(obj)
  }
  

  obj = FindNeighbors(obj, dims = dims.use, reduction=reduction )

  for (res in ress) {
    print(res)
    obj = FindClusters(object = obj, 
                        algorithm = algorithm,
                        resolution = res, 
                        verbose=F)
    obj[[paste(cluster_prefix, res, sep = ".")]] = Idents(object = obj)
  }
  obj
}

seurat_find_cluster_batch = function(obj, dims.use_list, ress=seq(.1, .6, .1), cluster_prefix, reduction="pca") {
  for (i in seq_along(dims.use_list)) {
    du = dims.use_list[[i]]
    na = names(dims.use_list)[i]
    obj = seurat_find_clusters(obj, 
                               dims.use = du,
                               ress=ress,
                               cluster_prefix = sprintf("%s_dim%s", cluster_prefix, na), 
                               algorithm = 2,
                               reduction = reduction)

  }
  obj
}

seurat_run_umap_batch = function(obj, dims.use_list) {
  for (i in seq_along(dims.use_list)) {
    du = dims.use_list[[i]]
    na = names(dims.use_list)[i]
    obj = RunUMAP(obj, 
                  dims=du, 
                  reduction.name = sprintf("umap%s", na), 
                  reduction.key=sprintf("umap%s_", na))
  }
  obj
}

seurat_run_alra = function(obj, assay=NULL) {
  if (is.null(assay)) {
    DefaultAssay(obj) = "SCT"
  } else {
    DefaultAssay(obj) = assay
  }

  slot(obj, 'tools')[["RunALRA"]] = NULL
  obj = RunALRA(obj)
  obj
}

seurat_integrate = function(obj_list, params, plot_dir = NULL) {
  options(future.globals.maxSize= 7000 * 1024^2)
  obj_list = map(obj_list, function(x) {
    VariableFeatures(x) = NULL
    x
  })
  obj_features = SelectIntegrationFeatures(object.list = obj_list, assay=c("SCT", "SCT"), nfeatures = 3000, fvf.nfeatures=3000)
  
  obj_list = PrepSCTIntegration(object.list = obj_list, 
                                 anchor.features = obj_features, 
                                 verbose = T)
  
  obj_anchors = FindIntegrationAnchors(object.list = obj_list, 
                                        assay = c("SCT", "SCT"),
                                        reference = 1,
                                        normalization.method = "SCT",
                                        anchor.features = obj_features, 
                                        reduction = "cca",
                                        k.anchor = params[["k.anchor"]],
                                        k.filter = params[["k.filter"]],
                                        k.score = params[["k.score"]],
                                        verbose = T)
  obj_int = IntegrateData(anchorset = obj_anchors,
                                     normalization.method = "SCT", 
                           verbose = TRUE)
  obj_int = RunPCA(obj_int, verbose = FALSE)
  obj_int = RunUMAP(obj_int, dims = 1:30, verbose = FALSE) 
  
  qsave(obj_int, file.path(plot_dir, "obj_int.qs"))
  return(obj_int)
}

seurat_integrate_plots = function(obj_list, obj_int, plot_dir=NULL, res_to_use = "cluster_int_sub2", groups_to_plot = c("celltype_comb", "celltype2", "dataset")) {
    walk(groups_to_plot, function(gp) {
      plots = DimPlot(obj_int, group.by = c(gp), label=T)
      plots = plots + theme(legend.position = "none") + guides(color = guide_legend(nrow = 3, 
                                                                                    byrow = TRUE, override.aes = list(size = 3)))
      save_plot(file.path(plot_dir, sprintf("umap_dimplot_%s.pdf", gp)),plots, base_height=6, base_aspect_ratio = 1.1)
    })
    #cts = na.omit(unique(obj_integrated_ex$celltype_n_ex))
    #colors = c("grey", qualitative_hcl(length(cts)))
    #set.seed(2)
    #colors = c("grey", randomColor(length(cts)))
    #names(colors) = c("tasic_2018", cts)
    
    to_plot = FetchData(obj_list[[2]], res_to_use)
    to_plot = unique(to_plot[,1])
    print(to_plot)
    
    md = FetchData(obj_int, c("celltype_comb", "dataset"))
    md = md %>% mutate(celltype_comb_1 = if_else(dataset=="tasic_2018", "tasic_2018", celltype_comb))
    obj_int$celltype_comb_1 = md$celltype_comb_1
    
    walk(to_plot, function(ct_to_plot) {
      print(ct_to_plot)
      colors_cur = rep("grey", times=length(to_plot)+1)# randomColor(length(cts)))
      names(colors_cur) = c("tasic_2018", to_plot)
      colors_cur[ct_to_plot] = "red"
      
      
      md = obj_int@meta.data
      md = md %>%
        select(-rowname) %>%
        rownames_to_column()
      obj_umap = Embeddings(obj_int, reduction = "umap")
      obj_umap = obj_umap %>% as.data.frame() %>% 
        rownames_to_column(var="rowname") %>%
        left_join(md)
      
      
      gg_tas = ggplot(obj_umap) + 
        geom_point(aes(x = UMAP_1, y = UMAP_2, alpha=dataset, color=celltype_comb_1), size=.5) + 
        theme(legend.position="top",
              axis.text = element_blank(),
              axis.ticks= element_blank(),
              axis.line = element_blank()) + 
        labs(title="tasic et al. 2018", x="UMAP1", y="UMAP2") +
     #   scale_color_manual(values=colors_cur) + 
        scale_alpha_manual(values=c("tasic_2018" = .05, "bf"=1))
      save_plot(file.path(plot_dir, sprintf("bf_umap_%s.pdf", ct_to_plot)), gg_tas, base_height=6, base_aspect_ratio = 1.1)
      
      ##gg_legend = get_legend(gg_tas)
      # 
      # save_plot(file.path(exc_sub_dir, sprintf("exc_tas_bf_umap_%s_just_legend.pdf", ct_to_plot)), gg_legend, base_width = 10, base_height=6)
      # 
      # gg_tas = gg_tas +
      #   theme(legend.position =  "none")
      # save_plot(file.path(exc_sub_dir, sprintf("exc_tas_bf_umap_%s_no_legend.pdf", ct_to_plot)), gg_tas, base_height=6, base_aspect_ratio = 1.1)
    })
}
seurat_find_transfer_anchors = function(obj_list, obj_int, params, plot_dir=NULL) {
  obj_trx_anchors = FindTransferAnchors(reference = obj_list[[1]], 
                                        query =obj_list[[2]],
                                        reduction = "pcaproject",
                                        normalization.method = "SCT", 
                                        dims = 1:30,
                                        reference.assay = "SCT",
                                        query.assay = "SCT",
                                        features=VariableFeatures(obj_int), 
                                        verbose = T,
                                        k.anchor = params[["k.anchor"]],
                                        k.filter = params[["k.filter"]],
                                        k.score = params[["k.score"]])
  qsave(obj_trx_anchors, file.path(plot_dir, "transfer_anchors.qs"))
  return(obj_trx_anchors)
}

seurat_transfer_labels = function(obj_list, obj_trx_anchors, params_trx,  reference_label="cluster", plot_dir=NULL) {
  refdata = FetchData(obj_list[[1]], reference_label)
  predictions = TransferData(anchorset = obj_trx_anchors, 
                             refdata = as.character(refdata[,1]), 
                             dims = 1:30,
                             weight.reduction = "pcaproject",
                             k.weight = params_trx[["k.weight"]],
                             sd.weight = params_trx[["sd.weight"]])
  saveRDS(predictions, file.path(plot_dir, sprintf("predictions_%s.rds",reference_label)))
  return(predictions)
}

seurat_transfer_labels_plots = function(obj_list, predictions, reference_label, plot_dir = NULL) {
  
  ## Prediction max density plot
  pdf(file.path(plot_dir, sprintf("prediction_%s_density.pdf", reference_label)))
  plot(density(predictions$prediction.score.max))
  dev.off()
  
  ## Entropy of classification
    
}

