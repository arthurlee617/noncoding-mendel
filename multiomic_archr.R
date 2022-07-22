library(ArchR)
setwd('~/multiome')
inputFiles <- c('~/multiome/multiome_R1_CN347SMN_e115/outs/atac_fragments.tsv.gz',
                '~/multiome/multiome_R2_CN347SMN_e115/outs/atac_fragments.tsv.gz')
names(inputFiles) <- c('multiome_R1',
                       'multiome_R2')

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 0, # keep ALL barcodes for now
  filterFrags = 0, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
addArchRGenome('mm10')
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchRProj_multiome",
  copyArrows = TRUE
)

saveArchRProject(proj)

#filter to valid barcodes only
md.1 <- read.table('~/multiome/multiome_R1_CN347SMN_e115/outs/per_barcode_metrics.csv', sep = ',', quote = '', header = T)
md.2 <- read.table('~/multiome/multiome_R2_CN347SMN_e115/outs/per_barcode_metrics.csv', sep = ',', quote = '', header = T)

md.1$archr_barcode <- paste('multiome_R1', md.1$barcode, sep = '#')
md.2$archr_barcode <- paste('multiome_R2', md.2$barcode, sep = '#')
md.1 <- md.1 %>% dplyr::filter(is_cell == 1)
md.2 <- md.2 %>% dplyr::filter(is_cell == 1)
proj <- proj[(proj$cellNames %in% md.1$archr_barcode) | (proj$cellNames %in% md.2$archr_barcode), ]
saveArchRProject(proj)
# Get gene integration for multiomic 
# compare p2g links for multiomic
# compare marker genes for multiomic vs archr project 

seRNA <- import10xFeatureMatrix( input = c('multiome_R1_CN347SMN_e115/outs/filtered_feature_bc_matrix.h5',
                                           'multiome_R2_CN347SMN_e115/outs/filtered_feature_bc_matrix.h5'),
                                 names = c('multiome_R1',
                                           'multiome_R2'))
 #https://github.com/GreenleafLab/ArchR/issues/507
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI"
)

peaks <- read.table('/ArchR/beds/shendure.bed', col.names =c('chr','start','end')) %>% makeGRangesFromDataFrame()
proj<- addPeakSet(proj, peaks)
proj<- addPeakMatrix(proj)

proj <- addPeak2GeneLinks(proj, reducedDims = 'IterativeLSI', useMatrix = 'GeneExpressionMatrix', maxDist = 500000)

saveArchRProject(proj)


multi <- readRDS('multiomics_cn347smn_e115_signac.rds') # load Signac object to get predicted cluster ID's 
pred.id <- multi$predicted.id
cells <- gsub('r2', 'multiome_R2', gsub('r1',"multiome_R1",gsub('_',"#", colnames(multi))))
names(pred.id) <- cells
pred.id <- pred.id[proj$cellNames]
proj$predictedCluster <- pred.id

cluster <- multi$seurat_clusters
cells <- gsub('r2', 'multiome_R2', gsub('r1',"multiome_R1",gsub('_',"#", colnames(multi))))
names(cluster) <- cells
cluster <- cluster[proj$cellNames]
proj$seuratCluster <- cluster

namedCluster <- multi$namedClusters
cells <- gsub('r2', 'multiome_R2', gsub('r1',"multiome_R1",gsub('_',"#", colnames(multi))))
names(namedCluster) <- cells
namedCluster <- namedCluster[proj$cellNames]
proj$namedCluster <- namedCluster

proj <- addUMAP(proj)

plotEmbedding(proj,name = 'namedCluster')


saveArchRProject(proj)




proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

proj <- addUMAP(proj, reducedDims = 'LSI_RNA', name='UMAP_RNA')

proj <- addCombinedDims(proj, reducedDims = c("IterativeLSI", "LSI_RNA"), name =  "LSI_Combined")
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined")

proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)
plotEmbedding(proj, name = 'Clusters', embedding = 'UMAP_Combined', discreteSet = 'ironMan' ) 

cM <- confusionMatrix(proj$Clusters, proj$predictedCluster)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
proj$namedCluster <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)
saveArchRProject(proj)
