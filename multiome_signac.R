library(Signac)
library(Seurat)

setwd("~/multiome")
seuMulti <- readRDS('multiome-signac-noArchR.rds')
#### R2 atac

fragpath <- "~/nextseq/Multiome_CMN347SMN_e11-5_R2/outs/fragments.tsv.gz"
metadata <- read.csv(
  file = "~/nextseq/Multiome_CMN347SMN_e11-5_R2/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
) %>% dplyr::filter(is__cell_barcode == 1)
atac.frags <- CreateFragmentObject(path = fragpath, cells = rownames(metadata))

library(future)
plan('multisession', workers = 32)
counts <- FeatureMatrix(
  fragments = atac.frags,
  features = granges(seuMulti),
  cells = rownames(metadata)
)

atac.assay <- CreateChromatinAssay(
  counts = counts,
  min.features = 0,
  fragments = atac.frags
)

r2 <- CreateSeuratObject(counts = atac.assay, assay = "peaks")

# compute LSI
r2 <- FindTopFeatures(r2, min.cutoff = 10)
r2 <- RunTFIDF(r2)
r2 <- RunSVD(r2)
r2$dataset <- 'multiome_R2'
seuMulti$dataset <- 'multiome_R1'


## add predicted ID info for replicate 2
seuATAC <- RunUMAP(seuATAC, reduction = "lsi", dims = 2:30, return.model = TRUE)
r2 <- RunUMAP(r2, reduction = "lsi", dims = 2:30, return.model = TRUE)
seuATAC$seurat_cluster <-seuATAC$Sample
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = seuATAC,
  query = r2,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

# map query onto the reference dataset
r2 <- MapQuery(
  anchorset = transfer.anchors,
  reference = seuATAC,
  query = r2,
  refdata = seuATAC$seurat_cluster ,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)

DimPlot(r2, reduction = 'umap', group.by = 'predicted.id', label = T)

merged <- merge(seuMulti, r2)

merged <- FindTopFeatures(merged, min.cutoff = 10)
merged <- RunTFIDF(merged)
merged <- RunSVD(merged)
merged <- RunUMAP(merged, reduction = 'lsi', dims = 2:30)
p1 <- DimPlot(merged, group.by = "dataset")


DimPlot(merged, group.by = 'predicted.id', size = 1)