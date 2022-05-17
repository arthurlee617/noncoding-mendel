library(ArchR)
library(Seurat)
library(dplyr)
setwd('~/lauren/ArchR/subprojects')



library(Seurat)
library(ArchR)
library(dplyr)
library(tidyr)
setwd('/home/adigioia/lauren/ArchR/subprojects')
library(stringr)
addArchRGenome("mm10")
## BE SURE TO SET FORCE = FALSE TO NOT OVERWRITE ARROW FILES
inputFiles <- c('/home/adigioia/scATAC/finalWorkflow/clusters/cluster0/cluster0.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster1/cluster1.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster2/cluster2.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster3/cluster3.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster4/cluster4.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster5/cluster5.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster6/cluster6.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster7/cluster7.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster8/cluster8.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster9/cluster9.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster10/cluster10.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster11/cluster11.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster12/cluster12.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster13/cluster13.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster14/cluster14.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster15/cluster15.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster16/cluster16.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster17/cluster17.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster18/cluster18.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster19/cluster19.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster20/cluster20.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster21/cluster21.bam',
                '/home/adigioia/scATAC/finalWorkflow/clusters/cluster22/cluster22.bam'
)
names(inputFiles) <- c('cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5', 'cluster6', 'cluster7', 'cluster8', 'cluster9', 'cluster10',
                       'cluster11', 'cluster12', 'cluster13', 'cluster14', 'cluster15', 'cluster16', 'cluster17', 'cluster18', 'cluster19', 'cluster20', 'cluster21', 'cluster22')

bamFlag = list(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,  hasUnmappedMate = NA,isMinusStrand = TRUE, isMateMinusStrand = NA, isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA, isNotPassingQualityControls = NA, isDuplicate = NA)
ArrowFiles <- createArrowFiles(inputFiles =inputFiles, sampleNames =names(inputFiles),filterTSS = 0, addTileMat = TRUE, 
                               addGeneScoreMat = TRUE, gsubExpression= ":.*",bamFlag = bamFlag, filterFrags = 0, force = FALSE)

proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "clusters_shendure", copyArrows = TRUE ) ### 85655 / 86089 cells

#conda activate umap
#python /home/adigioia/scATAC/finalWorkflow/scripts/runUmap.py /home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.SVDs.txt /home/adigioia/lauren/allSamples_filtered.umap2d.txt
#conda deactivate
### add UMAP coordinates generated from 
## do not filter doublets, etc; all QC / cell filtering has already happened at this point

umap <- read.table('/home/adigioia/lauren/allSamples_filtered.umap2d.txt')
colnames(umap) <- c('UMAPx', 'UMAPy')
cells <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.cells.txt')
umap$barcode <- cells$V1
celldata <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_finalFiltered.cellData.txt', header = TRUE)
celldata <- celldata %>% select(-c(UMAP1, UMAP2, UMAP3))

celldata <- merge(celldata, umap)
#celldata <- celldata %>% separate(barcode, c('barcode', NA), sep = '-', remove = TRUE)
celldata$cluster <- paste('cluster', celldata$cluster, sep = '')
celldata <- celldata %>% unite(barcode, c(cluster, barcode),  sep = '#')
celldata <- celldata %>% select(barcode, UMAPx, UMAPy)
rownames(celldata) <- celldata$barcode
celldata <- celldata %>% select(UMAPx, UMAPy)
colnames(celldata) <- c('shendure#UMAP_Dimension_1', 'shendure#UMAP_Dimension_2')
celldata <- celldata[proj$cellNames,]
proj@embeddings[['UMAP']] <- SimpleList(
  df = celldata, 
  params = NULL
)



### REDUCED DIMS FOR ARCHR PROJECT - get SVD from shendure workflow
# adigioia@rdt01343:~/lauren$ cp '/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.SVDs.txt' .
# adigioia@rdt01343:~/lauren$ cut -f2-50 allSamples_filtered.SVDs.txt > allSamples_filtered_trimmed.SVDs.txt
svds <- read.table('/home/adigioia/lauren/allSamples_filtered_trimmed.SVDs.txt')
svds <- as.matrix(svds)
cells <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.cells.txt')
colnames(cells) <- 'barcode'
celldata <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.cellData.txt', header = TRUE)
celldata <- merge(celldata, cells)
#celldata <- celldata %>% separate(barcode, c('barcode', NA), sep = '-', remove = TRUE)
celldata$cluster <- paste('cluster', celldata$cluster, sep = '')
celldata <- celldata %>% unite(barcode, c(cluster, barcode),  sep = '#')
celldata <- celldata %>% select(barcode)
rownames(svds) <- celldata$barcode
svds <- svds[celldata$barcode %in% proj$cellNames,]
svds <- svds[proj$cellNames,]
colnames(svds) <- gsub('V', 'LSI', colnames(svds))
proj@reducedDims[['IterativeLSI']] <- SimpleList(
  matSVD = svds,
  nDimensions = 49,
  outliers = NULL,
  scaleDims = NA,
  corToDepth = NA
)

cells <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.cells.txt')
colnames(cells) <- 'barcode'
celldata <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_filtered.cellData.txt', header = TRUE)
celldata <- merge(celldata, cells)
#celldata <- celldata %>% separate(barcode, c('barcode', NA), sep = '-', remove = TRUE)
celldata$cluster <- paste('cluster', celldata$cluster, sep = '')
celldata <- celldata %>% unite(barcode, c(cluster, barcode),  sep = '#')
rownames(celldata) <- celldata$barcode
celldata <- celldata[proj$cellNames,]
proj$Clusters <- proj$Sample
proj$dissection <- celldata$prefix
proj$celltype <- celldata$tissue
proj$bioSample <- paste(celldata$tissue, celldata$time, sep = '_')
proj$GFP <- grepl('neg', celldata$tissue)
proj$development <- as.factor(celldata$time)

plotEmbedding(proj, colorBy = 'cellColData', name = 'Sample')
saveArchRProject(proj)

##Add peakset, add peakset annotations
library(regioneR)
SamplePeaks <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_finalFiltered.peaks.sorted.bed')
colnames(SamplePeaks) <- c('chr', 'start', 'end')
SamplePeaks <- toGRanges(SamplePeaks, genome="mm10")

proj <- addPeakSet(proj, peakSet = SamplePeaks)
proj <- addPeakMatrix(proj)

saveArchRProject(proj) ### peaks added

### add annotations

peaks <- getPeakSet(proj)
GenomeAnnotation <- getGenomeAnnotation(proj)
GeneAnnotation <- getGeneAnnotation(proj)

fastAnnoPeaks <- function(
  peaks = NULL, 
  BSgenome = NULL, 
  geneAnnotation = NULL, 
  promoterRegion = c(2000, 100),
  logFile = NULL
){
  
  #Validate
  peakSummits <- resize(peaks,1,"center")
  BSgenome <- validBSgenome(BSgenome)
  
  #First Lets Get Distance to Nearest Gene Start
  distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
  mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
  promoters <- extendGR(resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
  type <- rep("Distal", length(peaks))
  type[which(og & oe)] <- "Exonic"
  type[which(og & !oe)] <- "Intronic"
  type[which(op)] <- "Promoter"
  mcols(peaks)$peakType <- type
  
  #First Lets Get Distance to Nearest TSS's
  distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToTSS <- mcols(distTSS)$distance
  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
  }
  
  #Get NucleoTide Content
  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  peaks
  
}

peaks <- fastAnnoPeaks(peaks, BSgenome = "mm10", geneAnnotation = GeneAnnotation, 
                       promoterRegion = c(2000, 100), logFile = NULL)


proj@peakSet$GC <- peaks$GC
proj@peakSet$nearestTSS <- peaks$nearestTSS
proj@peakSet$peakType <- peaks$peakType
proj@peakSet$nearestGene <- peaks$nearestGene
proj@peakSet$distToGeneStart <- peaks$distToGeneStart
proj@peakSet$distToTSS <- peaks$distToTSS

saveArchRProject(proj)





proj <- loadArchRProject("clusters_shendure")

#Cluster3: CN12
#Cluster4: SMN10
#Cluster5:SMN11
#Cluster9: SMN11

proj <- proj[!(proj$Sample %in% c('cluster3', 'cluster4', 'cluster5', 'cluster9')),]
saveArchRProject(proj, outputDirectory = 'clusters_shendure_drop_SMN_CN12')
proj <- loadArchRProject("clusters_shendure_drop_SMN_CN12")

integrated<- readRDS('~/lauren/scRNA/scRNA_integrated_transgeneregression.R')
proj <- loadArchRProject('~/lauren/ArchR/subprojects/clusters_shendure_drop_SMN_CN12')
### USE THIS DIM REDUCTION FOR GENE INTEGRATION - manually added dim reduction doesn't seem to be working
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI_ArchR", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

proj<- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI_ArchR", 
  name = "UMAP_ArchR", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)



proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_ArchR",
  seRNA = integrated,
  addToArrow = TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  force = TRUE
)

proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI_ArchR"
)

saveArchRProject(proj)


cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.0,
  resolution = 1,
  returnLoops = FALSE
)

peaks <- data.frame(metadata(cA)[[1]])
peaks$queryHits <- rownames(peaks)
colnames(peaks)[c(1,2,3)] <- c('chr_a', 'start_a', 'end_a')
peaks2 <- data.frame(metadata(cA)[[1]])
peaks2$subjectHits <- rownames(peaks2)
colnames(peaks2)[c(1,2,3)] <- c('chr_b', 'start_b', 'end_b')
cA <- merge(cA, peaks, by = 'queryHits')
cA <- merge(cA, peaks2, by = 'subjectHits')
cA <- data.frame(cA)
cA <- cA %>% select(chr_a, start_a, end_a, chr_b, start_b, end_b, correlation)




plotEmbedding(proj, embedding = 'UMAP', colorBy = 'cellColData', name ='')


groups <- proj$predictedGroup_Un
names(groups) <- proj$cellNames

proj2 <- loadArchRProject('~/lauren/ArchR/subprojects/cluster2')
proj2$predictedGroup <- groups[proj2$cellNames]
plotEmbedding(proj2, name = 'predictedGroup', pal = paletteDiscrete(unique(integrated$seurat_clusters), set = "stallion"))


markerGenes <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markerGenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
for (i in 1:19){
  clust <- names(markerList@listData)[i]
  a <- data.frame(markerList@listData[i])
  write.table(a,paste('~/lauren/ArchR/subprojects/markergenes/Gene_Integration_markergenes_', clust, '.tsv', sep = ''), quote = F, row.names = F, col.names = TRUE, sep = '\t')
  }

### All gene integration analysis uses '~/lauren/ArchR/subprojects/clusters_shendure_drop_SMN_CN12'

proj <- addPeak2GeneLinks(proj, reducedDims = "IterativeLSI_ArchR")

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj$CNclusters <- if_else(proj$Sample %in% c('cluster2', 'cluster8', 'cluster11'), 'CN7', 
                           if_else(proj$Sample %in% c('cluster7', 'cluster10', 'cluster19'), 'CN34', 
                                   if_else(proj$Sample %in% c('cluster12', 'cluster6'), 'CN6', proj$Sample)))


markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
motifs <- gsub(' .*', '', heatmapEM@column_names_param$labels)
### plot footprints of all of these motifs
proj <- loadArchRProject('clusters_shendure_drop_SMN_CN12')
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "Sample")


seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[motifs], 
  groupBy = "Sample"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Enriched-Footprints",
  addDOC = FALSE,
  smoothWindow = 5
)




markersPeaks2 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "CNclusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs2 <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")


proj <- addBgdPeaks(proj, method = 'ArchR')

proj<- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
saveArchRProject(proj)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

motifs <- c( 'Onecut1', 'Onecut2', 'Tcf21', 'Hoxb1', 'Tead4', 'Myod1')
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Sample", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)


# Identify deviant TF motifs: averaged by clusters here

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Sample")

# Subset to get just deviation z scores; identify maximum delta in z score between all clusters
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

# Identify TF's whose motif accesibillity is correlated with their own gene activity (use gene integration matrix and ArchR reduced dims)
corGSM_MM <- correlateMatrices(
  ArchRProj = proj,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI_ArchR"
)

# Add max delta deviation to correlation data frame - annotate each motif w/ maximum delta observed between clusters
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# Identify positive TF regulators: TFs whose correlation between motif and gene expression is > 0.5 (p adj < .01, max inter-clluster difference in dev Z score in top quartile)


corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

## Dot plot

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  ) 

p


### Plot deviations



motifs <- sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- markerMotifs[!(grepl('Srebf', markerMotifs))]

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

p <- plotGroups(ArchRProj = proj, 
                groupBy = "Sample", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow =3, rel_widths = c(2, rep(1, length(p2) - 1))),p2))



### Plotting motif enrichment
proj <- addImputeWeights(proj, reducedDims = 'IterativeLSI_ArchR')
markerMotifs <- c("z:Ebf1_90",     "z:Elf1_285",    "z:Elf4_283",    "z:Erg_287",     "z:Ets1_284",    "z:Ets2_828",    "z:Etv2_270",    "z:Fli1_277",
                  "z:Gata2_383",   "z:Hoxa2_420",  
                  "z:Hoxa3_584" ,  "z:Lhx4_454",    "z:Meis1_425",   "z:Meis2_461",   "z:Meis3_525",    "z:Neurod4_787", "z:Nhlh1_86",    "z:Nhlh2_84",    "z:Onecut1_826",
                  "z:Onecut2_827", "z:Onecut3_245", "z:Pax5_694",    "z:Pax8_854",    "z:Pou2f2_616",    "z:Sox2_750",    "z:Sp1_142",    
                     "z:Tead3_867",   "z:Tead4_868" )
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = markerMotifs[1:10], 
     embedding = "UMAP",
    imputeWeights = NULL
  )
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = c('Ebf1', 'Elf1', 'Elf4', 'Erg', 'Ets1', 'Ets2', 'Etv2', 'Fli1', 'Gata2', 'Hoxa2'), 
  embedding = "UMAP",
  imputeWeights = NULL
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs[11:20]), 
  embedding = "UMAP",
  imputeWeights = NULL
)
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = c('Hoxa3', 'Lhx4', 'Meis1', 'Meis2', 'Meis3', 'Neurod4', 'Nhlh1', 'Nhlh2', 'Onecut1', 'Onecut2'), 
  embedding = "UMAP",
  imputeWeights = NULL
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs[21:28]), 
  embedding = "UMAP",
  imputeWeights = NULL
)
p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = c('Onecut3', 'Pax5', 'Pax8','Pou2f2','Sox2','Sp1','Tead3','Tead4'), 
  embedding = "UMAP",
  imputeWeights = NULL
)


p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))


##Motif footprinting

motifs <- c('Onecut3', 'Pax5', 'Pax8','Pou2f2','Sox2','Sp1','Tead3','Tead4','Hoxa3', 'Lhx4', 'Meis1', 'Meis2', 'Meis3', 'Neurod4', 'Nhlh1', 'Nhlh2', 'Onecut1', 'Onecut2',
            'Ebf1', 'Elf1', 'Elf4', 'Erg', 'Ets1', 'Ets2', 'Etv2', 'Fli1', 'Gata2', 'Hoxa2')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "PosTF-Footprints",
  addDOC = FALSE,
  smoothWindow = 5
)




ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

df <- data.frame(TF = rownames(CN6motifs), mlog10Padj = assay(CN6motifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp2 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

df <- data.frame(TF = rownames(CN34motifs), mlog10Padj = assay(CN34motifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp3 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf('CN_motifs.pdf')
ggUp + ggtitle('CN7')
ggUp2 + ggtitle('CN6')
ggUp3 + ggtitle('CN34')
dev.off()

splitondash <- function(x){
  x <- strsplit(x, "_")[[1]][1]
}
splitonspace <- function(x){
  x <- strsplit(x, " ")[[1]][1]
}
motifs <- heatmapEM2@column_names_param$labels
motifs <- lapply(motifs, splitonspace)
motifs <- unlist(motifs)
availGenes <- getFeatures(proj, useMatrix = "GeneIntegrationMatrix")
motifs <- motifs[!(grepl( paste0(notAvail, collapse = "|"), motifs))]
notAvail <- motifs[!(motifs %in% availGenes)]
motifs <- motifs[motifs %in% availGenes]
motifs <- c("Meis1", "Pbx1", "Hoxb1", "Lhx6")
motifs_matrix <- getFeatures(proj, select = paste(motifs, collapse = "|"), useMatrix = "MotifMatrix")
motifs_matrix <- grep("z:", motifs_matrix, value = TRUE)

p <- plotEmbedding(proj, colorBy = "MotifMatrix", name = motifs_matrix)

pdf(file = "clusterEnrichedMotifs_MotifMatrix.pdf")
for (i in seq(1, length(p), 9)){
  p1 <- p[i:(i+8)]
  p1c <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  
  do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
}
dev.off()



p <- plotEmbedding(proj, colorBy = "GeneIntegrationMatrix", name = motifs)

pdf(file = "clusterEnrichedMotifs_GeneExpression.pdf")
for (i in seq(1, length(p), 9)){
  p1 <- p[i:(i+8)]
  p1c <- lapply(p1, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  
  do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
}
dev.off()





motifs <- c('Nkx61')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

mat <- assays(enrichMotifs)[['mlog10Padj']]
mat <- mat[markerMotifs,]
ArchR:::.ArchRHeatmap(
mat =mat,
scale = FALSE,
limits = c(0, max(mat)),
color = paletteContinuous(set = "comet", n = 100),
clusterCols = FALSE,
clusterRows = FALSE,
labelRows = TRUE,
useRaster = TRUE,
fontSizeCols = 6,
borderColor = FALSE,
customColLabel = seq_len(ncol(mat)),
showRowDendrogram = FALSE,
draw = FALSE,
name = "Norm. Enrichment -log10(P-adj) [0-Max]"
) 

motifPositions <- getPositions(proj)
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype", minReplicates = 2, maxReplicates = 8, force = TRUE)
motifs <- c('Foxl1', 'Gata2', 'Gata3')
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))



markerMotifs
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "celltype"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Divide",
  plotName = "Pbx1",
  addDOC = FALSE,
  smoothWindow = 5
)

motifs <- c('Alx4', "Ets2", "Erg", "Irf3", "Mafa", "Maz", "Nkx", "Phox2a", "Phox2b","Rfx2", "Rfx4")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "celltype"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Divide",
  plotName = "elementE_select_motifs",
  addDOC = FALSE,
  smoothWindow = 5
)

proj$GFPclust <- if_else(proj$Sample %in% c('cluster0', 'cluster1', 'cluster14', 'cluster16', 'cluster18', 'cluster20', 'cluster21', 'cluster22'), 
                         'GFPneg', if_else(proj$Sample %in% c('cluster2', 'cluster6', 'cluster8', 'cluster10', 'cluster11', 'cluster12', 'cluster15', 'cluster19'),
                                           'GFPpos', proj$Sample), proj$Sample)


markersPeaks2 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "GFPclust",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs2 <- peakAnnoEnrichment(
  seMarker = markersPeaks2,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 20, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")


### SUPERCLUSTER MOTIF ENRICHMENT
proj$superclusters <- if_else(proj$Sample %in% c('cluster19', 'cluster10', 'cluster7'), 'A',
                              if_else(proj$Sample %in% c('cluster6', 'cluster12'), 'B',
                                      if_else(proj$Sample %in% c('cluster11', 'cluster8', 'cluster2'), 'C',
                                              if_else(proj$Sample %in% c('cluster3'), 'D',
                                                      if_else(proj$Sample %in% c('cluster17', 'cluster16', 'cluster15', 'cluster5'), 'E',
                                                              if_else(proj$Sample %in% c('cluster13'), 'F',
                                                                      if_else(proj$Sample %in% c('cluster22'), 'WW',
                                                                              if_else(proj$Sample %in% c('cluster21', 'cluster9'), 'XX',
                                                                                      if_else(proj$Sample %in% c('cluster14'), 'YY',
                                                                                              if_else(proj$Sample %in% c('cluster20', 'cluster18'), 'ZZ', 'NONE'))))))))))

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "superclusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs2 <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 20, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")



motifs <- gsub(' .*', '', heatmapEM@column_names_param$labels)
### plot footprints of all of these motifs
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "superclusters", force = TRUE)


seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[motifs], 
  groupBy = "superclusters"
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Superclusters-Enriched-Footprints",
  addDOC = FALSE,
  smoothWindow = 5
)


seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[motifs], 
  groupBy = "superclusters",
  useGroups = c('A', 'B', 'C', 'E','F','WW','YY','ZZ')
)


plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Superclusters-Enriched-Footprints-removeXX",
  addDOC = FALSE,
  smoothWindow = 5
)



### PLOTS
##'~/lauren/scRNA/scRNA_integrated_transgeneregression.R'
DimPlot(integrated, cols = paletteDiscrete(unique(integrated$seurat_clusters), set = "stallion"), label = TRUE) # RNA_integrated_clusters.pdf
DimPlot(integrated, group.by = 'CN') #RNA_integrated_CN
FeaturePlot(integrated, features = 'Isl1', min.cutoff = 'q1')

 #https://github.com/satijalab/seurat/issues/2087



early_2 <- subset(integrated, subset= seurat_clusters %in% c("0", "7", "3", "2", "24","22","1", "14", "28", "20", "25") )
early_2$major_cluser <- early_2$seurat_clusters
early_2 <- FindVariableFeatures(early_2, assay = 'RNA', nfeatures = 3000)
early_2 <- ScaleData(early_2, verbose = FALSE, features = early_2[['RNA']]@var.features)
early_2 <- RunPCA(early_2, npcs = 50, verbose = FALSE)
early_2 <- RunUMAP(early_2, reduction = "pca", dims = 1:50)
early_2 <- FindNeighbors(early_2, reduction = "pca", dims = 1:50)
early_2 <- FindClusters(early_2, resolution = 1.8)

DimPlot(early_2, group.by = 'major_cluser',cols = paletteDiscrete(unique(integrated$seurat_clusters))) + DimPlot(early_2,
                                                                                                                 cols = paletteDiscrete(unique(early_2$seurat_clusters), set = 'bear')) + DimPlot(early_2, group.by = 'CN')
all.markers <- FindAllMarkers(early_2)

plotEmbedding(proj, colorBy ='cellColData', name = 'predictedGroup_Un', pal =  paletteDiscrete(unique(integrated$seurat_clusters)))

heatmap(prop.table(table( proj$Clusters,proj$predictedGroup_Un)))


### Get marker genes for all subclusters
proj <- loadArchRProject('clusters_shendure_drop_SMN_CN12')
celldata <- read.table('/home/adigioia/scATAC/finalWorkflow/finalSamples/allSamples_finalFiltered.cellData.txt', header = TRUE)
celldata$cluster <- paste('cluster', celldata$cluster, sep = '')
celldata <- celldata %>% unite(barcode, c(cluster, barcode),  sep = '#')
rownames(celldata) <- celldata$barcode
celldata <- celldata[proj$cellNames,]
proj$subcluster <- celldata$subclusterNum

setwd('~/lauren/ArchR/subprojects/clusters_projects')
for (cluster in unique(proj$Sample)){
  print(cluster)
  subproj <- proj[proj$Sample == cluster,]

  markersGI <- getMarkerFeatures(
    ArchRProj = subproj, 
    useMatrix = "GeneIntegrationMatrix", 
    groupBy = "subcluster",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersGI, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
  
  for(subclust in unique(subproj$subcluster)){
    print(paste('On ', cluster, subclust, sep = ' '))
    a <- data.frame(markerList[subclust][[1]])
    filename <- paste(cluster,  subclust, 'markergenes.tsv', sep = "_")
    a <- a %>% select(name, Log2FC)
    write.table(a, filename, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  
}


pal <- c(ArchRPalettes$ironMan, ArchRPalettes$grove)
pal <- pal[c(1:23)]
names(pal) <- sort(unique(proj$Clusters))
 
pdf('umap_bysample_final.pdf')
plotEmbedding(proj, name = "dissection", labelMeans = F, discreteSet = "circus") + ggtitle("Collored by Dissection")
plotEmbedding(proj, name = "bioSample", labelMeans = F, discreteSet = "bear") + ggtitle("Colored By Celltype / Time")
plotEmbedding(proj, name = "Clusters", pal= pal, labelMeans = F) + ggtitle("Colored By Cluster")
plotEmbedding(proj, name = "celltype", discreteSet = "grove", labelMeans = F) + ggtitle("Colored By Celltype")
dev.off()

proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)


#### MOTIF ENRICHMENT


### TF_enrichment_bycluster.pdf

proj$celltype <- as.character(proj$celltype)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
pdf('test.pdf')
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()


## GFP positive celltypes
CN7 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CN7",
  bgdGroups = c('CN12', "CN6","CN3_CN4","SMN")
)

CN6 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CN6",
  bgdGroups = c('CN12', "CN7","CN3_CN4","SMN")
)
CN34 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CN3_CN4",
  bgdGroups = c('CN12', "CN7","CN6","SMN")
)

CN7motifs <- peakAnnoEnrichment(
  seMarker = CN7,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
CN6motifs <- peakAnnoEnrichment(
  seMarker = CN6,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
CN34motifs <- peakAnnoEnrichment(
  seMarker = CN34,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(CN7motifs), mlog10Padj = assay(CN7motifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))


proj$occular <- if_else(proj$celltype %in% c("CN3_CN4", "CN6"), "occular", proj$celltype)
## Occular vs lower MN
Occular <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "occular",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = 'occular',
  bgdGroups = c('CN12', "SMN", "CN7")
)

occularMotifs<- peakAnnoEnrichment(
  seMarker = Occular,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(occularMotifs), mlog10Padj = assay(occularMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))



#Midbrain vs hindbrain
proj$CN34 <- if_else(proj$celltype %in% c("CN3_CN4", "CN3_CN4-neg"), "CN34", proj$celltype)
grouped <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "CN34",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = 'CN34',
  bgdGroups = c('CN12',"CN7", "CN6", "CN12-neg", "CN7-neg", "CN6-neg")
)

groupedMotifs<- peakAnnoEnrichment(
  seMarker = grouped,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)


df <- data.frame(TF = rownames(groupedMotifs), mlog10Padj = assay(groupedMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp2 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))



proj$celltype <- as.character(proj$celltype)
grouped_neg <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = 'CN3_CN4-neg',
  bgdGroups = c("CN12-neg","CN7-neg", "CN6-neg")
)

groupednegMotifs<- peakAnnoEnrichment(
  seMarker = grouped_neg,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(groupednegMotifs), mlog10Padj = assay(groupednegMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp3 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))



#branchial
CN7 <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CN7",
  bgdGroups = c('CN12', "CN6","CN3_CN4")
)
branchial <- peakAnnoEnrichment(
  seMarker = CN7,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)




df <- data.frame(TF = rownames(branchial), mlog10Padj = assay(branchial)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp4 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))



proj$celltype <- as.character(proj$celltype)
grouped_pos <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = 'CN3_CN4',
  bgdGroups = c("CN12","CN7", "CN6")
)

groupedposMotifs<- peakAnnoEnrichment(
  seMarker = grouped_pos,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(groupedposMotifs), mlog10Padj = assay(groupedposMotifs)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp5 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))





pdf('Motif_enrichment_pairwise_comparisons.pdf')
ggUp + ggtitle("Occular vs lower MNs: \nCN3,4,6 GFP+ vs CN7,12,SMN GFP+")
ggUp2 + ggtitle("Midbrain vs Hindbrain: \nCN3,4 Pos/Neg vs CN6,7,12 Pos/Neg")
ggUp5 + ggtitle("Midbrain vs Hindbrain: \nCN3,4 Pos vs CN6,7,12 Pos")
ggUp3 + ggtitle("Midbrain vs Hindbrain: \nCN3,4 Neg vs CN6,7,12 Neg")
ggUp4 + ggtitle("Somatic motor vs Branchial motor: \nCN3,4,6,12 GFP+ vs CN7 GFP+")
dev.off()





markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

tokeep <- 'Zfp148|Sp2|Klf7|Tcf21|Myod1|Myog|Ferd3l|Ascl1|Atoh8|Myf5|Tcfap4|Msc|Nhlh2|Zfp263|Tead4|Foxl1|Foxd1|Elf1|Sfpi1|Ehf|Spic|Fli1|Hoxb1|Ebf1|Pax9|Zeb1|Dbx2|Ahctf1|Onecut3|Bbx|Hoxb6|Gata|Smarcc1|Fos|Bach2|Jund|Alx3|Hoxc5|Lmx1b|Hlx|Prop1|Sebox|Nkx62|Vax1|Lhx3|Meis2|Pknox1|Gsc2|Pitx2|Obox1|Otx1|Sox2|Pou5f1|Sox9|Irf1|Sry'

allmot <- rownames(enrichMotifs)[grep(tokeep, rownames(enrichMotifs))]
cleanmot <- gsub('_.*','',allmot) #Zfp263 appears 4 times; Foxl1 appears 28! times; all others are OK. Foxl1-355 appears in original heatmap; Zfp263-876 appears in original heatmap; use these
tokeep <- 'Zfp148|Sp2|Klf7|Tcf21|Myod1|Myog|Ferd3l|Ascl1|Atoh8|Myf5|Tcfap4|Msc|Nhlh2|Zfp263_876|Tead4|Foxl1_355|Foxd1|Elf1|Sfpi1|Ehf|Spic|Fli1|Hoxb1|Ebf1|Pax9|Zeb1|Dbx2|Ahctf1|Onecut3|Bbx|Hoxb6|Gata|Smarcc1|Fos|Bach2|Jund|Alx3|Hoxc5|Lmx1b|Hlx|Prop1|Sebox|Nkx62|Vax1|Lhx3|Meis2|Pknox1|Gsc2|Pitx2|Obox1|Otx1|Sox2|Pou5f1|Sox9|Irf1|Sry'
tokeep <- rownames(enrichMotifs)[grep(tokeep, rownames(enrichMotifs))]
tokeep <- tokeep[-c(11,13,14)]

enrichMotifs <- enrichMotifs[, c('cluster7', 'cluster10', 'cluster19', 'cluster6', 'cluster12', 'cluster2', 'cluster8', 'cluster11', 'cluster3', 'cluster4', 'cluster15', 'cluster16', 'cluster17', 'cluster13', 'cluster1', 'cluster0', 'cluster22', 'cluster9', 'cluster21', 'cluster14', 'cluster18', 'cluster20')]
heatmapEM2 <- plotEnrichHeatmapCustom(enrichMotifs, tokeep = tokeep, transpose = TRUE, clusterCols = F)
pdf('test.pdf', width = 10)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

heatmapEM2 <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE, returnMatrix = T)
pdf('test-2.pdf', width = 15, par(mar=c(4,4,0.5,0.5)))
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()





### DORCS


proj <- loadArchRProject('clusters_shendure_drop_SMN_CN12')

proj <- addPeak2GeneLinks(proj, reducedDims = 'IterativeLSI_ArchR', maxDist = 50000)
p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF <- data.frame(p2geneDF) %>% dplyr::filter(FDR < .0001) %>% dplyr::filter(Correlation > .1)
genes <- getGeneAnnotation(proj)$genes
TSS <- getGeneAnnotation(proj)$TSS


peakset <- data.frame(getPeakSet(proj))
pdf('dist_to_tss_hist.pdf')
hist(peakset$distToTSS, breaks = 100)
dev.off()
quantile(peakset$distToTSS, seq(0,1,.1))
peakset <- peakset %>% tidyr::unite(peakName, seqnames, start, end, sep = "_")
peakset <- peakset %>% dplyr::select(peakName, peakType, width, nearestGene, distToGeneStart, distToTSS)
p2geneDF <- p2geneDF %>% dplyr::select(peakName, geneName, Correlation, FDR)
p2geneDF <- merge(p2geneDF, peakset, all.x = T, by = 'peakName')
p2geneDF_trim <- unique(setDT(p2geneDF)[order(FDR)], by = "peakName") # for any peaks that link to multiple genes, only select peak with smallest FDR
p2geneDF_grouped <- p2geneDF_trim %>% group_by(geneName) %>% summarise(n_peaks = n_distinct(peakName))
p2geneDF_grouped <- p2geneDF_grouped[order(p2geneDF_grouped$n_peaks),]
p2geneDF_grouped$x <- seq(1,11943)
p2geneDF_grouped <- data.frame(p2geneDF_grouped)
p2geneDF_grouped$label <- if_else(p2geneDF_grouped$n_peaks >10, T, F)
pdf('test.pdf')
ggplot(p2geneDF_grouped, aes(x=x, y = n_peaks)) + geom_text_repel(data = p2geneDF_grouped[p2geneDF_grouped$label == T,], aes(label = geneName)) + labs(title = "P2G Associations: \nIdentifying DORCs") + geom_point() + theme_minimal() + geom_hline(yintercept = 10.5, linetype = 'dashed') + xlab('') + ylab('Number of Peak-Gene Links')
dev.off()


# Cell read normalization: peak x cell matrix 
pm <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
mat <- pm@assays@data$PeakMatrix
mat_norm <- mat %*% Matrix::Diagonal(x = 1 / sqrt(Matrix::colSums(mat^2))) # normalize by cell total count
rm(mat)
row_names <- data.frame(pm@rowRanges)
row_names <- row_names %>% unite(peakName, seqnames, start, end)
peakNames <- row_names$peakName
rownames(mat_norm) <- peakNames


# recalculate 500kb

#proj <- addPeak2GeneLinks(proj, reducedDims = 'IterativeLSI_ArchR', maxDist = 500000)
p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF <- data.frame(p2geneDF) %>% dplyr::filter(FDR < .0001) %>% dplyr::filter(Correlation > .1)

#identify DORCS
DORC <- p2geneDF_grouped[p2geneDF_grouped$label == T,]$geneName
p2geneDF <- p2geneDF %>% dplyr::filter(geneName %in% DORC)

#Make DORC mat for one DORC
DORC_gata <- p2geneDF %>% dplyr::filter(geneName == 'Gata2')
DORC_mat <- mat_norm[DORC_gata$peakName,]
DORC_mat <- Matrix::colSums(DORC_mat) #DORC_mat is matrix of cellwise DORC scores
proj$Gata2_DORC <-DORC_mat[proj$cellNames]

write.table(p2geneDF_grouped, 'DORC_genes.tsv', sep = '\t', quote = F, row.names = F, col.names = T)

## Relationships between peak properties and p2g corr 
# 
# peakset <- data.frame(getPeakSet(proj))
# peakset <- peakset %>% unite(peakName, seqnames, start, end)
# p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
# p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
# p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
# p2geneDF <- data.frame(p2geneDF) %>% dplyr::filter(FDR < .0001) %>% dplyr::filter(Correlation > .1)
# p2geneDF <- merge(p2geneDF, peakset, by = "peakName", all.x = T)
# p2geneDF_trim <- unique(setDT(p2geneDF)[order(FDR)], by = "peakName") # for any peaks that link to multiple genes, only select peak with smallest FDR
# 
# pdf('GC-vs-correlation.pdf')
# ggplot(p2geneDF_trim, aes(x = GC, y = Correlation)) + geom_point() + theme_minimal()
# dev.off()
# pdf('distToTSS-vs-correlation.pdf')
# ggplot(p2geneDF_trim, aes(x = distToTSS, y = Correlation)) + geom_point() + theme_minimal()
# dev.off()
# 
# corr <- cor.test(p2geneDF_trim$Correlation, p2geneDF_trim$distToTSS, method = 'spearman')
# 
# 
# pm <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
# mat <- pm@assays@data$PeakMatrix
# mat <- Matrix::rowSums(mat)
# row_names <- data.frame(pm@rowRanges)
# row_names <- row_names %>% unite(peakName, seqnames, start, end)
# peakNames <- row_names$peakName
# names(mat) <- peakNames
# mat <- mat[p2geneDF_trim$peakName]
# p2geneDF_trim$coverage <- mat
# pdf('coverage-vs-correlation.pdf')
# ggplot(p2geneDF_trim, aes(x = coverage, y = Correlation)) + geom_point() + theme_minimal()
# dev.off()
# 


g <- getMatrixFromProject(proj, useMatrix = "GeneIntegrationMatrix")
mat <- g@assays@data$GeneIntegrationMatrix
rownames(mat) <- g@elementMetadata$name
colnames(mat) <- g@colData@rownames
mat <- mat[DORC,]


cluster2 <- loadArchRProject('cluster2')
cluster2$DORC_gata <- DORC_mat[cluster2$cellNames]
pdf('cluster2_DORCgata.pdf')
plotEmbedding(cluster2, name = 'DORC_gata')
dev.off()
pdf('cluster2_expression_gata.pdf')
plotEmbedding(cluster2, name = 'Gata2', colorBy = 'geneintegrationmatrix')
dev.off()
cluster2 <- addImputeWeights(cluster2)
smoothed_DORC <- imputeMatrix(mat = t(as.matrix(DORC_mat[cluster2$cellNames])), imputeWeights = getImputeWeights(cluster2))
cluster2$smoothed_DORC <- smoothed_DORC
pdf('cluster2_DORCgata.pdf')
plotEmbedding(cluster2, name = 'smoothed_DORC')
dev.off()
pdf('cluster2_expression_gata.pdf')
plotEmbedding(cluster2, name = 'Gata2', colorBy = 'geneintegrationmatrix')
dev.off()




### feature plots for
proj <- loadArchRProject('clusters_shendure_drop_SMN_CN12')
proj<- addImputeWeights(proj, reducedDims = 'IterativeLSI_ArchR')
pdf('UMAP_geneint_lab_meeting_12_18_22.pdf')
plotEmbedding(proj, colorBy = 'geneintegrationmatrix',name="Prrxl1")
plotEmbedding(proj, colorBy = 'geneintegrationmatrix',name="Sox10")
plotEmbedding(proj, colorBy = 'geneintegrationmatrix',name="Mbp")
plotEmbedding(proj, colorBy = 'geneintegrationmatrix',name="Cnp")
dev.off()

proj <- addPeak2GeneLinks(proj, "IterativeLSI_ArchR", maxDist = 1000000)

p2g <- getPeak2GeneLinks(proj, corCutOff = .1)

ol<- findOverlaps(toGRanges('chr13:115498270-117109688'),p2g$Peak2GeneLinks)
p2g$Peak2GeneLinks<- p2g$Peak2GeneLinks[ol@to]
df <- data.frame(p2g$Peak2GeneLinks)


p<- plotBrowserTrack(proj, groupBy = 'Sample', geneSymbol = "Isl1", upstream = 800000, downstream = 800000, loops = p2g)

pdf('browser.pdf')
grid::grid.draw(p$Isl1)
dev.off()

p2g$Peak2GeneLinks <- p2g$Peak2GeneLinks[which((df$start == 116309688)|(df$end ==116309688))]
p<- plotBrowserTrack(proj, groupBy = 'Sample', geneSymbol = "Isl1", upstream = 800000, downstream = 800000, loops = p2g, plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
                     sizes = c(10, 1.5, 3, 4), useMatrix = "GeneIntegrationMatrix")
pdf('cluster_isl_broswerplot.pdf')
grid::grid.draw(p$Isl1)
dev.off()

proj$celltype <- as.character(proj$celltype)
p<- plotBrowserTrack(proj, groupBy = 'celltype', useGroups = c('CN3_CN4','CN6','CN7', 'CN3_CN4-neg', 'CN6-neg','CN7-neg'), geneSymbol = "Isl1", upstream = 800000, downstream = 800000, loops = p2g, plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
                     sizes = c(10, 1.5, 3, 4), useMatrix = "GeneIntegrationMatrix")
pdf('celltype_isl_browserplot.pdf')
grid::grid.draw(p$Isl1)
dev.off()
