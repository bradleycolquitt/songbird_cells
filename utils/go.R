library(biomaRt)
library(topGO)
library(GO.db)
library(org.Hs.eg.db)

select = dplyr::select
rename = dplyr::rename
filter = dplyr::filter

entrez_host = "jul2015.archive.ensembl.org"
#entrez_host = "www.ensembl.org"
bm = "ENSEMBL_MART_ENSEMBL"

get_ensembl = function(ds = "hsapiens_gene_ensembl") {
  useMart(biomart=bm, dataset=ds, host=entrez_host)
}

get_entrez_gene_id = function(ds = "hsapiens_gene_ensembl") {
  ensembl = useMart(biomart=bm, dataset=ds,  host=entrez_host)
  getBM(c('entrezgene','hgnc_symbol'), mart=ensembl)
}

#entrez_gene_id = get_entrez_gene_id()

entrez_gene_id = na.omit(get_entrez_gene_id())

#' Return associated GO terms using topGO package
#'
#' @param data_list, list of data.frames containing gene names and 
#' some measure of importance/signficance. Must also have 'gene_id' column 
#' @param value_name, name of associated value to run GO test
#' @param thresh, value to threshold value_name
#' @return
#' @export
#'
#' @examples
get_go_results = function(data_list, value_name, thresh=0.1) {
  stopifnot(value_name %in% colnames(data_list[[1]]))
  gos1 = lapply(1:length(data_list), function(i) {
    print(names(data_list)[i])
    tmp = as.data.frame(data_list[[i]])
    
    if (nrow(tmp)==0)
      return(NULL)
    genelist = tmp[, value_name]

    names(genelist) = entrez_gene_id[match(tmp$gene_id, entrez_gene_id[,2]), 1]
    print(head(genelist))
    selectTop = function(d, th=thresh) {
      return(d<th)
    }
    
    go_sets = c("BP", "MF", "CC")
    gos = lapply(go_sets, function(go) {
      go_bp = new("topGOdata",
                  description = "", ontology = go,
                  allGenes = genelist , geneSel = selectTop,
                  nodeSize = 10,
                  annot = annFUN.org, mapping="org.Hs.eg.db", ID="entrez")
      resultFisher = runTest(go_bp, algorithm = "classic", statistic = "fisher")
      resultKS = runTest(go_bp, algorithm = "classic", statistic = "ks")
      resultWeight01 = runTest(go_bp, algorithm = "weight01", statistic = "fisher")
      GenTable(go_bp, classicFisher = resultFisher,
               classicKS = resultKS, weight01Fisher =  resultWeight01,
               orderBy = "weight01Fisher", ranksOf = c("classicFisher", "classicKS"), topNodes = 10)
    })
    names(gos) = go_sets
    gos
  })
  names(gos1) = names(data_list)
  gos1 
}

#' Return associated GO terms using topGO package
#'
#' @param data, data.frame containing gene names and 
#' some measure of importance/signficance. Must also have 'gene_id' column 
#' @param value_name, name of associated value to run GO test
#' @param thresh, value to threshold value_name
#' @return
#' @export
#'
#' @examples
get_go_results_single = function(data, value_name, thresh=0.1, outdir=NULL) {
  stopifnot(value_name %in% colnames(data))
  tmp = as.data.frame(data)
    
  if (nrow(tmp)==0)
    return(NULL)
  genelist = tmp[, value_name]
  print(head(genelist))
  names(genelist) = entrez_gene_id[match(tmp$gene_id, entrez_gene_id[,2]), 1]
  
  selectTop = function(d, th=thresh) {
    return(d<th)
  }
  
  go_sets = c("BP", "MF", "CC")
  gos = lapply(go_sets, function(go) {
    go_bp = new("topGOdata",
                description = "", ontology = go,
                allGenes = genelist , geneSel = selectTop,
                nodeSize = 15,
                annot = annFUN.org, mapping="org.Hs.eg.db", ID="entrez")
    resultFisher = runTest(go_bp, algorithm = "classic", statistic = "fisher")
    resultKS = runTest(go_bp, algorithm = "classic", statistic = "ks")
    resultWeight01 = runTest(go_bp, algorithm = "weight01", statistic = "fisher")
    resultWeight01.ks = runTest(go_bp, algorithm = "weight01", statistic = "ks")
    resultWeight01.t = runTest(go_bp, algorithm = "weight01", statistic = "t")
    gt = GenTable(go_bp, classicFisher = resultFisher,
             classicKS = resultKS, 
             weight01Fisher =  resultWeight01,
             weight01FisherKS =  resultWeight01.ks,
             weight01FisherT =  resultWeight01.t,
             orderBy = "weight01Fisher", ranksOf = c("classicFisher", "classicKS"), topNodes = 10)
    
    sigGenes = sigGenes(go_bp)
    return(list(gt=gt, sigGenes=sigGenes))
  })
  
  names(gos) = go_sets
  
  gos_gts = map(gos, 1)
  gos_genes = entrez_to_gene_symbol(gos[[1]][[2]])
  
  if (!is.null(outdir)) {
    if(!dir.exists(outdir))
      dir.create(outdir, recursive=T)
    lapply(1:length(gos_gts), function(i) {
      fname = paste(outdir, paste(names(gos)[i], ".txt", sep=""), sep="/")
      write.table(gos_gts[[i]], fname, quote=F, sep="\t", row.names=F)
    })
  }
  
  gos_gts = bind_rows(gos_gts, .id="set")
  gos_gts = gos_gts %>% mutate(genes = map(GO.ID, ~get_genes_by_go(.x, ensembl )))
  gos_gts = gos_gts %>% mutate(genes_sig = map(genes, ~.x[.x %in% gos_genes]))
  gos_gts
}

#' Return associated GO terms using topGO package
#'
#' @param gene_subset, character vector of genes to test
#' @param gene_full, background or full set of genes
#' @return list of GO results
#' 
#' @export
#'
#' @examples
#
get_go_results2 = function(gene_subset, gene_full) {
    genelist = vector("numeric", length=length(gene_full))
    genelist[gene_full %in% gene_subset] = 1
    print(table(genelist))
    names(genelist) = entrez_gene_id[match(gene_full, entrez_gene_id[,2]), 1]
    
    go_sets = c("BP", "MF", "CC")
    gos = lapply(go_sets, function(go) {
      go_bp = new("topGOdata",
                  description = "", ontology = go,
                  allGenes = genelist , geneSel = function(x) {x==1},
                  nodeSize = 10,
                  annot = annFUN.org, mapping="org.Hs.eg.db", ID="entrez")
      resultFisher = runTest(go_bp, algorithm = "classic", statistic = "fisher")
      resultKS = runTest(go_bp, algorithm = "classic", statistic = "ks")
      resultWeight01 = runTest(go_bp, algorithm = "weight01", statistic = "fisher")
      GenTable(go_bp, classicFisher = resultFisher,
               classicKS = resultKS, weight01Fisher =  resultWeight01,
               orderBy = "weight01Fisher", ranksOf = c("classicFisher", "classicKS"), topNodes = 10)
    })
    names(gos) = go_sets
    bind_rows(gos, .id="set")
    
  #})
  #names(gos1) = names(data_list)
  #gos1 
}


gene_symbol_to_entrez = function(ids, na_rm=T) {
  if (na_rm) {
    as.character(na.omit(entrez_gene_id[match(ids, entrez_gene_id[,2]),1]))
  } else {
    as.character(entrez_gene_id[match(ids, entrez_gene_id[,2]),1])
  }
}

entrez_to_gene_symbol = function(ids, na_rm=T) {
  if (na_rm) {
    as.character(na.omit(entrez_gene_id[match(ids, entrez_gene_id[,1]),2]))
  } else {
    as.character(entrez_gene_id[match(ids, entrez_gene_id[,1]),2])
  }
}
gene_symbol_to_ensembl = function(ids, na_rm=T) {
  if (na_rm) {
    as.character(na.omit(entrez_gene_id[match(ids, entrez_gene_id[,2]),1]))
  } else {
    as.character(entrez_gene_id[match(ids, entrez_gene_id[,2]),1])
  }
}

get_genes_by_go = function(go_term, ensembl) {
  ag = getBM(c('entrezgene','hgnc_symbol'), filters='go_id', values=go_term, mart=ensembl)
  return(ag[,2])
}


process_gene_ids = function(ids) {
  g = ids
  g = str_replace_all(g, "\\*", "")
  g = str_replace(g, "[\\.|\\-][1-3]", "")
  g = unlist(str_split(g, ","))
  g = unique(g)
  return(g)
}

get_go_by_id = function(go_id) {
  ensembl = get_ensembl()
  entry = getBM(c('entrezgene','external_gene_name'), filters='go_id', values=go_id, mart=ensembl)
  entry = entry %>% mutate(external_gene_name = toupper(external_gene_name))
  entry
}

get_tf_genes = function() {
  ensembl = get_ensembl()
  entry = getBM(c('entrezgene','external_gene_name'), filters='go_id', values=c("GO:0003700", "GO:0003702", "GO:0003709", "GO:0016563", "GO:0016564"), mart=ensembl)
  entry = entry %>% mutate(external_gene_name = toupper(external_gene_name))
  entry
}

## List from https://www.sciencedirect.com/science/article/pii/S0092867418301065?via%3Dihub#app2, Lambert, Cell, 2018
get_tf_genes_human = function() {
  dat = read_csv("~/data/gene_lists/lambert_tfs/tfs.csv")
  colnames(dat) = c("ensembl_id", "external_gene_name", "DBD", "TF")
  dat = dat %>% filter(TF=="Yes")
  dat
}
get_axon_guidance_genes = function() {
  ensembl = get_ensembl()
  entry = getBM(c('entrezgene','external_gene_name'), filters='go_id', values='GO:0007411', mart=ensembl)
  entry = entry %>% mutate(external_gene_name = toupper(external_gene_name))
  entry
}

axon_guidance_regex = "^NETO|^NRP|^SEMA|^FEZ|^EPH|^EFN|^DLG|^CNTN|^NRXN|^PLXN|^SLIT|^ROBO|^NCAM|L1CAM|RELN|^UNC|CDH|^FLRT|^DPYS|^ITGA|^NTN|^CDC42|ABL1|^ABLIM|CDK5|^CFL|^CHP|DCC"

get_voltage_gated_ion_channels_genes = function() {
  ensembl = get_ensembl()
  entry = getBM(c('entrezgene','external_gene_name'), filters='go_id', values='GO:0005216', mart=ensembl)
  entry = entry %>% mutate(external_gene_name = toupper(external_gene_name))
  entry
}


get_axon_guidance_genes_kegg = function() {
  ensembl = get_ensembl()
  entry = getBM(c('entrezgene','external_gene_name'), filters='go_id', values='GO:0007411', mart=ensembl)
  entry = entry %>% mutate(external_gene_name = toupper(external_gene_name))
  entry
}

get_axon_guidance_genes_reactome = function() {
  
}


gene_lists_df = data.frame(short_label = c("ribosome", 
                                           "microglia",
                                           "astrocyte", 
                                           "Alzheimers",
                                           "neuron", 
                                           "mitochondria",
                                           "postsynaptic", 
                                           "neuron",
                                           "temporal lobe",
                                           "neuron",
                                           "mitochondria",
                                           "astrocyte",
                                           "neuron",
                                           "endosome",
                                           "ribosome",
                                           "ribosome",
                                           "aging",
                                           "mitochondria",
                                           "presynaptic",
                                           "Alzheimers",
                                           "synapse",
                                           "mitochondria",
                                           "microglia",
                                           "nuclear",
                                           "cortex",
                                           "neuron",
                                           "neuron",
                                           "schizophrenia",
                                           "excitatory neuron",
                                           "excitatory neuron",
                                           "autism",
                                           "ribosome",
                                           "synapse",
                                           "synapse",
                                           "cytoplasmic",
                                           "excitatory neuron", 
                                           "ER/Golgi",
                                           "synapse",
                                           "cortex",
                                           "cytoskeleton",
                                           "facial motor nucleus",
                                           "trochlear nucleus",
                                           "trochlear nucleus",
                                           "parietal lobe",
                                           "cingulate neurons",
                                           "interneurons",
                                           "mitochondria",
                                           "autism",
                                           "mitochondria",
                                           "proteosome"),
                           UserDefinedCategories = c("salmon_M12_Ribosome__HumanMeta",
                                                     "upWithABeta_MGactivation_GSE772__MicrogialMarkers",
                                                     "brown_M3_Astrocytes__HumanMeta",
                                                     "DownWithAlzheimers_Liang__ADvsCT_inCA1", "Neuron_probable__Cahoy",
                                                     "green_M5_Mitochondria__HumanMeta",
                                                     "PostSynapticDensity_proteins__Bayes",
                                                     "Neuron_definite__Cahoy",
                                                     "Temporal Lobe_localMarker(top200)_IN_Cerebral Cortex__HBA",
                                                     "blue_M16_Neuron__CTX",
                                                     "cyan_M4_Mitochondria__MouseMeta", 
                                                     "brown_M15_Astrocyte__CTX",
                                                     "red_M11_Neuron__HumanMeta", 
                                                     "Recycling_endosome_TGN_Foster__MO", 
                                                     "midnightblue_M2_Ribosome__CTX",
                                                     "salmon_M12_Ribosome__MouseMeta", 
                                                     "blue_downAging_mitochondria_synapse__Lu_Aging",
                                                     "brown_downAD_mitochondrion__Blalock_AD",
                                                     "PresynapticCompartmentProteins_Morciano__MO",
                                                     "DownWithAlzheimers_Blalock__ADvsCT_inCA1",
                                                     "Synaptic__MitochondrialType",
                                                     "Mitochondria_Foster__MO",
                                                     "pink_M10_Microglia(Type1)__HumanMeta",
                                                     "turquoise_M14_Nucleus__MouseMeta",
                                                     "blue_Cortex__HumanChimp",
                                                     "tan_M13_Neuron__MouseMeta",
                                                     "Neuron__ABA",
                                                     "Schizophrenia_possible__DiseaseGenes",
                                                     "brown_pyramidalNeurons_Layer5/basolateralAmygdala__Sugino/Winden",
                                                     "GlutamatergicNeuronsInMouseCortex_Sugino__MO",
                                                     "Autism_associated_module_M12__Voineagu",
                                                     "yellow_noChangeAD_antigenProcessing_ribosome__Blalock_AD",
                                                     "pink_downAD_synapticTransmission__Blalock_AD", 
                                                     "purple_downAD_synapticTransmission__Blalock_AD",
                                                     "Cytoplasm_Foster__MO",
                                                     "pink_M14_GlutamatergicSynapticFunction__CTX",
                                                     "ER_Golgi_vesicles_Foster__MO",
                                                     "magenta_downAD_synapticTransmission__Blalock_AD",
                                                     "brown_Cortex__HumanChimp",
                                                     "turquoise_downAD_intracellularTransport_cytoskeleton__Blalock_AD",
                                                     "facial motor nucleus_localMarker(top200)_IN_Pontine Tegmentum__HBA",
                                                     "trochlear nucleus_localMarker(FC>2)_IN_Mesencephalon__HBA",
                                                     "trochlear nucleus_localMarker(top200)_IN_Mesencephalon__HBA",
                                                     "Parietal Lobe_localMarker(top200)_IN_Cerebral Cortex__HBA",
                                                     "turquoise_CingulateNeurons(allTypes)/Layer5__Sugino/Winden",
                                                     "blue_Interneurons_Hippocampus(Sst+)/Cingulate(Pvalb+)__Sugino/Winden",
                                                     "Somatic__MitochondrialType",
                                                     "Autism_differential_expression_across_at_least_one_comparison__Voineagu",
                                                     "cyan_M4_Mitochondria__HumanMeta",
                                                     "Proteasome_Foster__MO"))