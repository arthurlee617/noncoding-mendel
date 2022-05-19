library(Seurat)
setwd('~/lauren/ArchR/scRNAseq')

dir <- 'CN6_E11-5A/mm10_cur'
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN6-e115')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)
hist(as.numeric(seu['g.eGFP',]@assays$RNA@counts))
FeaturePlot(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 20)], features = 'g.eGFP')
seu <- subset(seu, cells = colnames(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 20)]))
cn6 <- AverageExpression(seu, assays = 'RNA')
a <- cn6$RNA
write.table(a, 'CN6_e115_expression.tsv', sep = '\t', quote = F, row.names = T, col.names = F)


dir <- 'CN34_E11_5B/mm10_cur'
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN34-e115')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)
hist(as.numeric(seu['g.eGFP',]@assays$RNA@counts))
FeaturePlot(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 1)], features = 'g.eGFP')
seu <- subset(seu, cells = colnames(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 1)]))
cn34 <- AverageExpression(seu, assays = 'RNA')
a <- cn34$RNA
write.table(a, 'CN34_e115_expression.tsv', sep = '\t', quote = F, row.names = T, col.names = F)


dir <- 'scRNA_CN7_e11-5/rep1_5-23-17'
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN7-e115')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)
hist(as.numeric(seu['g.eGFP',]@assays$RNA@counts))
FeaturePlot(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 5)], features = 'g.eGFP')
seu <- subset(seu, cells = colnames(seu[,which(as.numeric(seu['g.eGFP',]@assays$RNA@counts) > 5)]))
cn711 <- AverageExpression(seu, assays = 'RNA')
a <- cn711$RNA
write.table(a, 'CN7_e115_expression.tsv', sep = '\t', quote = F, row.names = T, col.names = F)







dir <- 'scRNA_CN7_e10-5/rep1_6-12-17/posOnlyCells'
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN7-e105')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)
cn710 <- AverageExpression(seu, assays = 'RNA')
a <- cn710$RNA
write.table(a, 'CN7_e105_expression.tsv', sep = '\t', quote = F, row.names = T, col.names = F)



dir <- 'scRNA_CN7_e10-5/rep2_8-15-19' #no GFP 
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN7-e115')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)

dir <- 'scRNA_CN7_e10-5/rep3_9-11-19' #no GFP 
ob <- Read10X(dir)
seu <- CreateSeuratObject(ob, project = 'CN7-e115')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:50)
