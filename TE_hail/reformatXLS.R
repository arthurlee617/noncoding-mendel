# Take .tsv of TE calls and convert from tall to wide format for hail matrix table import

#de novos
setwd('/hail-TE')
TE <- read.csv('TE-denovos.csv', header=TRUE)
R <- colnames(TE)
TE$Chr <- str_remove(TE$Chr,'chr')
TE$locus <- paste(TE$Chr, TE$Breakpoint, sep=':')
library(stringr)
#TE$famID <- str_split(TE$Sample, '-', simplify = TRUE)
#TE$famID.2 <- NULL
#TE$famID.3 <- NULL
#TE$famID.4 <- NULL
TE$Breakpoint <- NULL
TE <- unique(TE)
library(tidyr)

TE_wide <- spread(TE, Sample, Features)
TE_wide$Features  <- NULL
TE_wide <- unique(TE_wide)
write.table(TE_wide, file='TE-denovos.tsv', quote=FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)

setwd('/mt-searches')
TE <- read.table('20201001_CCDD_TE_Raw_genome-wide_rare_ins_BZ.tsv', header=FALSE, sep = "\t")
setwd('/hail-TE')
colnames(TE) <- R
TE$Chr <- str_remove(TE$Chr,'chr')
TE$locus <- paste(TE$Chr, TE$Breakpoint, sep=':')
library(stringr)
#TE$famID <- str_split(TE$Sample, '-', simplify = TRUE)
#TE$famID.2 <- NULL
#TE$famID.3 <- NULL
#TE$famID.4 <- NULL
TE$Breakpoint <- NULL
TE <- unique(TE)
library(tidyr)

TE_wide <- spread(TE, Sample, Features)
TE_wide$Features  <- NULL
TE_wide <- unique(TE_wide)
write.table(TE_wide, file='TE-900Genomes.tsv', quote=FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
