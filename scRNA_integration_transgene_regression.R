library(Seurat)
library(dplyr)
library(patchwork)
setwd('/scRNA')

#####
# Integrated analysis
#####
data.dirs = c('CN6_E10-5A/mm10_cur','CN6_E11-5A/mm10_cur', 
              'scRNA_CN7_e10-5/rep1_6-12-17/negOnlyCells', 
              'scRNA_CN7_e10-5/rep1_6-12-17/posOnlyCells', 
              'scRNA_CN7_e10-5/rep3_9-11-19', 
              'scRNA_CN7_e11-5/rep1_5-23-17', 
              'CN34_E10_5B/mm10_cur', 
              'CN34_E10-5A/mm10_cur', 
              'CN34_E10-5C_PosNeg/mm10_cur', 
              'CN34_E11-5A/mm10_cur', 
              'CN34_E11_5B/mm10_cur')

dir1 <- 'CN6_E10-5A/mm10_cur'
dir2 <- 'CN6_E11-5A/mm10_cur'
dir3 <- 'scRNA_CN7_e10-5/rep1_6-12-17/negOnlyCells'
dir4 <- 'scRNA_CN7_e10-5/rep1_6-12-17/posOnlyCells'
dir5 <- 'scRNA_CN7_e10-5/rep3_9-11-19'
dir6 <- 'scRNA_CN7_e11-5/rep1_5-23-17'
dir7 <- 'CN34_E10_5B/mm10_cur'
dir8 <- 'CN34_E10-5A/mm10_cur'
dir9 <- 'CN34_E10-5C_PosNeg/mm10_cur'
dir10 <- 'CN34_E11-5A/mm10_cur'
dir11 <- 'CN34_E11_5B/mm10_cur'

ob1 <- Read10X(data.dir = dir1)
ob2 <- Read10X(data.dir = dir2)
ob3 <- Read10X(data.dir = dir3)
ob4 <- Read10X(data.dir = dir4)
ob5 <- Read10X(data.dir = dir5)
ob6 <- Read10X(data.dir = dir6)
ob7 <- Read10X(data.dir = dir7)
ob8 <- Read10X(data.dir = dir8)
ob9 <- Read10X(data.dir = dir9)
ob10 <- Read10X(data.dir = dir10)
ob11 <- Read10X(data.dir = dir11)

seu1 <- CreateSeuratObject(counts = ob1, project = "CN6-10", min.cells = 3, min.features = 200)
seu2 <- CreateSeuratObject(counts = ob2, project = "CN6-11", min.cells = 3, min.features = 200)
seu3 <- CreateSeuratObject(counts = ob3, project = "CN7-10-neg", min.cells = 3, min.features = 200)
seu4 <- CreateSeuratObject(counts = ob4, project = "CN7-10-pos", min.cells = 3, min.features = 200)
seu5 <- CreateSeuratObject(counts = ob5, project = "CN7-10-r3", min.cells = 3, min.features = 200)
seu6 <- CreateSeuratObject(counts = ob6, project = "CN7-11", min.cells = 3, min.features = 200)
seu7 <- CreateSeuratObject(counts = ob7, project = "CN34-10-B", min.cells = 3, min.features = 200)
seu8 <- CreateSeuratObject(counts = ob8, project = "CN34-10-A", min.cells = 3, min.features = 200)
seu9 <- CreateSeuratObject(counts = ob9, project = "CN34-10-C", min.cells = 3, min.features = 200)
seu10 <- CreateSeuratObject(counts = ob10, project = "CN34-11-A", min.cells = 3, min.features = 200)
seu11 <- CreateSeuratObject(counts = ob11, project = "CN34-11-B", min.cells = 3, min.features = 200)

obj.list <- c(seu1, seu2, seu3, seu4, seu5, seu6, seu7, seu8, seu9, seu10, seu11)

obj.list <- lapply(X = obj.list, FUN = function(x) {
  ###REMOVE geGFP and tdTomato if they exist #https://github.com/satijalab/seurat/issues/2610
  counts <- GetAssayData(x, assay = "RNA")
  counts <- counts[-(which(rownames(counts) %in% c('g.eGFP', 'tdTomato'))),]
  x <- subset(x, features = rownames(counts))
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  x <- subset(x, subset = percent.mt < 5)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- ScaleData(x)
})




#feats <- c("nFeature_RNA", "nCount_RNA")
#VlnPlot(alldata, features = "nFeature_RNA")
#VlnPlot(alldata, features = "nCount_RNA")
obj.list[8][[1]] <- subset(obj.list[8][[1]], subset = nFeature_RNA > 4000) # Cellranger web summary: Low Fraction Reads in Cells	56.4%	Ideal > 70%. Application performance may be affected. Many of the reads were not assigned to cell-associated barcodes. This could be caused by high levels of ambient RNA or by a significant population of cells with a low RNA content, which the algorithm did not call as cells. The latter case can be addressed by inspecting the data to determine the appropriate cell count and using --force-cells.
# Cut off tail distribution of low nfeature_RNA cells; 661 samples -> 586 
obj.list[9][[1]]<- subset(obj.list[9][[1]], subset = nFeature_RNA > 1500) #No noted QC issues. cut off tail dist of low nfeature_RNA cells. 1884 -> 1871
obj.list[10][[1]] <- subset(obj.list[10][[1]], subset = nFeature_RNA > 3000) #801->754
obj.list[5][[1]] <- subset(obj.list[5][[1]], subset = nFeature_RNA > 3000) # 4511 -> 4219

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 10000)
  x <- ScaleData(x)
})

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 10000)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchorset = anchors)

## stop here 

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated , verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:50)
integrated$CN <- if_else(grepl('CN34', integrated$orig.ident), 'CN34', if_else(grepl('CN6', integrated$orig.ident), 'CN6', if_else(grepl('CN7', integrated$orig.ident), 'CN7', 'none')))
integrated <- RunTSNE(integrated, dims =1:50, perplexity = 50)

integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:50)
integrated <- FindClusters(integrated, resolution = 1.2)
saveRDS(integrated, '/scRNA/scRNA_integrated_transgeneregression.rds')
