library(cicero)
library(monocle3)
library(ArchR)
library(batchelor)

## Clusters 2 and 11 (CN7 based clusters)
proj <- loadArchRProject('/ArchR/subprojects/clusters_shendure_drop_SMN_CN12')
proj <- proj[proj$Sample %in% c('cluster2', 'cluster11')]
## Call new peaks on just this cluster; iterative sub-clsering - re-cluster (subclusters) + call peaks
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 5000,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "dissection")
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "dissection"
)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "dissection")


proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "subClusters",
  resolution = 0.8
)

plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "subClusters")
saveArchRProject(proj, outputDirectory = '/ArchR/subprojects/cluster2_11')


pathToMacs2 <- '/miniconda3/envs/ArchR/bin/macs2'
proj<- addGroupCoverages(ArchRProj = proj, groupBy = "subClusters")

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "subClusters", 
  pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)

## move over to cicero CDS

peaks <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
indata <- peaks@assays$data$PeakMatrix
cellinfo <- getCellColData(proj)
row.names(cellinfo) <- proj$cellNames
peakinfo <- data.frame(getPeakSet(proj))
peakinfo$site_name <- paste(peakinfo$seq, peakinfo$start, peakinfo$end, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
peakinfo <- data.frame(peakinfo)
cellinfo <- data.frame(cellinfo)

indata@x[indata@x > 0] <- 1
input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo,gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds@reducedDims$LSI <- proj@reducedDims$IterativeLSI$matSVD
input_cds@reducedDims$Aligned<- proj@reducedDims$Harmony$matDR
umap <- as.matrix(proj@embeddings$UMAP$df)
colnames(umap) <-NULL
input_cds@reducedDims$UMAP <- umap
plot_cells(input_cds, reduction_method = "UMAP", color_cells_by = "dissection", cell_size = 1)
input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds)


cds <- input_cds

ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

cell <- data.frame(getCellColData(proj))
cell$UMAPx <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_1`
cell$UMAPy <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_2`


ggplot(cell, aes(x = UMAPx, y =  UMAPy, color = Sample)) + geom_point(size = 0.5, stroke = .3, shape = 19) + theme_classic() + theme(legend.position="bottom")  +
geom_segment(aes_string(x="source_prin_graph_dim_1",
                        y="source_prin_graph_dim_2",
                        xend="target_prin_graph_dim_1",
                        yend="target_prin_graph_dim_2"), size=.75,color=I("grey28"), linetype="solid", na.rm=TRUE, data=edge_df)

saveRDS(edge_df, "cluster2_11_traj.rds")
saveRDS(cell, "cluster2_11_celldata.rds")
a <- getMatrixFromProject(proj, "GeneIntegrationMatrix")
saveRDS(a, "cluster2_11_gene.rds")



###
## CN3,4,7 based clusters
proj <- loadArchRProject('/ArchR/subprojects/clusters_shendure_drop_SMN_CN12')
proj <- proj[proj$Sample %in% c('cluster8', 'cluster19', 'cluster10', 'cluster2', 'cluster11')]
## Call new peaks on just this cell group; iterative sub-clustering - re-cluster (subclusters) + call peaks
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 5000,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "dissection")
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "dissection"
)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "Sample")


proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "subClusters",
  resolution = 0.8
)

plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "subClusters")
saveArchRProject(proj, outputDirectory = '/ArchR/subprojects/cluster8_19_10_2_11')
proj <- loadArchRProject("cluster8_19_10_2_11")

pathToMacs2 <- '/miniconda3/envs/ArchR/bin/macs2'
proj<- addGroupCoverages(ArchRProj = proj, groupBy = "subClusters")

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "subClusters", 
  pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)



peaks <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
indata <- peaks@assays$data$PeakMatrix
cellinfo <- getCellColData(proj)
row.names(cellinfo) <- proj$cellNames
peakinfo <- data.frame(getPeakSet(proj))
peakinfo$site_name <- paste(peakinfo$seq, peakinfo$start, peakinfo$end, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
peakinfo <- data.frame(peakinfo)
cellinfo <- data.frame(cellinfo)

indata@x[indata@x > 0] <- 1
input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo,gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds@reducedDims$LSI <- proj@reducedDims$IterativeLSI$matSVD
input_cds@reducedDims$Aligned<- proj@reducedDims$Harmony$matDR
umap <- as.matrix(proj@embeddings$UMAP$df)
colnames(umap) <-NULL
input_cds@reducedDims$UMAP <- umap
plot_cells(input_cds, reduction_method = "UMAP", color_cells_by = "dissection", cell_size = 1)
input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds)


cds <- input_cds

ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

cell <- data.frame(getCellColData(proj))
cell$UMAPx <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_1`
cell$UMAPy <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_2`


ggplot(cell, aes(x = UMAPx, y =  UMAPy, color = Sample)) + geom_point(size = 0.5, stroke = .3, shape = 19) + theme_classic() + theme(legend.position="bottom")  +
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"), size=.75,color=I("grey28"), linetype="solid", na.rm=TRUE, data=edge_df)

saveRDS(edge_df, "cluster8_19_10_2_11traj.rds")
saveRDS(cell, "cluster8_19_10_2_11celldata.rds")


#####


proj <- loadArchRProject('/ArchR/subprojects/clusters_shendure_drop_SMN_CN12')
proj <- proj[proj$Sample %in% c('cluster0', 'cluster7', 'cluster14', 'cluster15', )]
## Call new peaks on just this cluster; iterative sub-clsering - re-cluster (subclusters) + call peaks
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 5000,
  force = TRUE
)

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "dissection")
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "dissection"
)
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "Harmony", 
  name = "UMAP",
  force = TRUE
)
plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "Sample")


proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "subClusters",
  resolution = 0.8
)

plotEmbedding(proj, embedding = "UMAP", colorBy = "cellColData", name = "subClusters")
saveArchRProject(proj, outputDirectory = '/ArchR/subprojects/cluster8_19_10_2_11')


pathToMacs2 <- '/miniconda3/envs/ArchR/bin/macs2'
proj<- addGroupCoverages(ArchRProj = proj, groupBy = "subClusters")

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "subClusters", 
  pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)



peaks <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
indata <- peaks@assays$data$PeakMatrix
cellinfo <- getCellColData(proj)
row.names(cellinfo) <- proj$cellNames
peakinfo <- data.frame(getPeakSet(proj))
peakinfo$site_name <- paste(peakinfo$seq, peakinfo$start, peakinfo$end, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
peakinfo <- data.frame(peakinfo)
cellinfo <- data.frame(cellinfo)

indata@x[indata@x > 0] <- 1
input_cds <-  suppressWarnings(new_cell_data_set(indata,cell_metadata = cellinfo,gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

input_cds@reducedDims$LSI <- proj@reducedDims$IterativeLSI$matSVD[proj$cellNames,]
input_cds@reducedDims$Aligned<- proj@reducedDims$Harmony$matDR
umap <- as.matrix(proj@embeddings$UMAP$df)
colnames(umap) <-NULL
input_cds@reducedDims$UMAP <- umap
plot_cells(input_cds, reduction_method = "UMAP", color_cells_by = "dissection", cell_size = 1)
input_cds <- cluster_cells(input_cds)
input_cds <- learn_graph(input_cds)


cds <- input_cds

ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]

edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select_(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       source="sample_name",
                       source_prin_graph_dim_1="prin_graph_dim_1",
                       source_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select_(
                       target="sample_name",
                       target_prin_graph_dim_1="prin_graph_dim_1",
                       target_prin_graph_dim_2="prin_graph_dim_2"),
                   by = "target")

cell <- data.frame(getCellColData(proj))
cell$UMAPx <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_1`
cell$UMAPy <- proj@embeddings$UMAP$df$`Harmony#UMAP_Dimension_2`


ggplot(cell, aes(x = UMAPx, y =  UMAPy, color = Sample)) + geom_point(size = 0.5, stroke = .3, shape = 19) + theme_classic() + theme(legend.position="bottom")  +
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"), size=.75,color=I("grey28"), linetype="solid", na.rm=TRUE, data=edge_df)

saveRDS(edge_df, "cluster8_19_10_2_11traj.rds")
saveRDS(cell, "cluster8_19_10_2_11celldata.rds")


