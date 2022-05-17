library(dplyr)
library(tidyr)
setwd('~/lauren/hail/hail-files/SV')

args = commandArgs(trailingOnly=TRUE)
cn <- args[1]
infile <- paste('SV-results/', cn, '-SV-inher-aggregated', sep ='')
results <- read.table(paste(infile, '.tsv', sep =''), header=TRUE, sep = '\t')

bed <- results %>% dplyr::select(SVloc, name)
bed <- bed %>% separate(SVloc, c('chr', 'ranges'), sep = ':')
bed <- bed %>% separate(ranges, c('start', 'end'), sep = '-')
outfile <- paste(infile, '.bed', sep = '')
write.table(bed, outfile, sep = '\t', quote = FALSE, row.names=FALSE, col.names=FALSE)