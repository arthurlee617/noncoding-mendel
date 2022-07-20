library(GenomicRanges)
library(tidyr)
library(dplyr)
library(regioneR)
setwd('~/hail/hail-files/SV')
inh <- read.table('annot.inhsubset.filter.bed', header=TRUE, comment.char = '')
colnames(inh)[1] <- 'chr'
inhGR <- toGRanges(inh)

args = commandArgs(trailingOnly=TRUE)
bed <- args[1]
CN <- args[2]
pheno <- args[3]
outfile <- paste('SV-results/', CN, '-SV-inher-exploded.tsv', sep = '')
cn <- read.table(bed)
colnames(cn) <- c('chr', 'start', 'end', 'cn')
cn <- toGRanges(cn)
findOverlaps(inhGR, cn)


overlaps <- findOverlapPairs(inhGR, cn)
SV <- overlaps@first
cns <- overlaps@second

SV <- data.frame(SV)
cns <- data.frame(cns)
SV$peak <- paste(cns$seqnames, ':', cns$start, '-', cns$end, sep='')
SV$cn <- cns$cn


SV <- SV %>% mutate(samples = strsplit(as.character(samples), ",")) %>% unnest(samples)
SV <- SV %>% separate(samples, c(NA, 'famID', NA), sep = 'S|-', remove=FALSE)

pheno <- read.table('../engle_ccdd_032619.ped', header=TRUE)
pheno <- data.frame(famID = pheno$FAMILY, pheno = pheno$PHENOTYPE)
pheno <- unique(pheno)

SV <- merge(SV, pheno, by = 'famID')

# if (phenotype == 'cfeom'){
#   phenotype <- 'CFEOM'} else if 
# (phenotype == 'cfp'){
#   phenotype <- 'CFP'} else if 
# (phenotype == 'ptosis'){
#   phenotype <- 'Ptosis'} else if 
# (phenotype == 'mgjw'){
#   phenotype <- 'PtosisMGJWS'} else if 
# (phenotype == 'drs'){
#   phenotype <- 'DRS/HGP'} else if 
# (phenotype == 'fnp'){
#   phenotype <- 'CFEOM/FNP'}else 
#     {phenotype <- 'Moebius'}

  
#SV <- SV %>% filter(pheno == phenotype)
SV <- SV %>% unite("SVloc", seqnames:start, remove = TRUE, sep = ":")
SV <- SV %>% unite("SVloc", SVloc:end, remove = TRUE, sep = "-")

write.table(SV, outfile, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
