setwd('~/multiome')
library(ArchR)
library(parallel)
addArchRGenome("mm10")

proj <- loadArchRProject('/ArchR/clusters_shendure')



saveArchRProject(proj, outputDirectory = 'clusters_shendure_copy')





#####
mat <- getMatrixFromProject(proj, useMatrix = 'PeakMatrix')
mat_matrix <- mat@assays@data$PeakMatrix
rownames(mat_matrix) <- paste(mat@rowRanges@seqnames, mat@rowRanges@ranges, sep = '-')
seu <- CreateChromatinAssay(counts= mat_matrix)
ccd <- data.frame(getCellColData(proj))
seuATAC <- CreateSeuratObject(seu, assay = "ATAC", meta.data=ccd)


seuATAC <- FindTopFeatures(seuATAC, min.cutoff = 10)
seuATAC <- RunTFIDF(seuATAC)
seuATAC <- RunSVD(seuATAC)
seuATAC <- RunUMAP(seuATAC, reduction.model = 'LSI')

seuATAC$dataset <- "ATAC"
saveRDS(seuATAC, 'atac-clusters-signac.rds')
