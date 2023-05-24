## general goal: clusters of mainly GFP-positive cells are more pure than clusters of mainly GFP-negative cells

############# combine time points but keep GFP+/- separate
## this one probably works best (after 90 degree transpose)
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")

# get rid of all columns except the cluster and class
keeps <- c("cluster", "tissue")
b <- a[keeps]

b$tissue <- factor(b$tissue, c("CN3_CN4", "CN6", "CN7", "CN12", "SMN", "CN3_CN4-neg", "CN6-neg", "CN7-neg", "CN12-neg", "SMN-neg")) 

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$tissue, b$cluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("heatmap.by.tissue.gfp.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap(ccmatrixpercent, 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "cluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
dev.off()

############# separate time points and keep GFP+/- separate
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")

# get rid of all columns except the cluster and class
keeps <- c("cluster", "class")
b <- a[keeps]

b$class <- factor(b$class, c("CN3_CN4-10.5", "CN3_CN4-11.5", "CN6-10.5", "CN6-11.5", "CN7-10.5", "CN7-11.5", "CN12-10.5", "CN12-11.5", "SMN-10.5", "SMN-11.5", "CN3_CN4-neg-10.5", "CN3_CN4-neg-11.5", "CN6-neg-10.5", "CN6-neg-11.5", "CN7-neg-10.5", "CN7-neg-11.5", "CN12-neg-10.5", "CN12-neg-11.5", "SMN-neg-10.5", "SMN-neg-11.5")) 


# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$class, b$cluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("heatmap.by.tissue.time.gfp.ordered.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap((ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "cluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
dev.off()

############# separate all replicates
## cluster membership is not determined by replicate; consistent with not batch effects
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")

# get rid of all columns except the cluster and class
keeps <- c("cluster", "prefix")
b <- a[keeps]

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$prefix, b$cluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("heatmap.by.replicate.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap(ccmatrixpercent, 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "cluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
dev.off()


############# separate all replicates and plot against subcluster
## cluster membership is not determined by replicate; consistent with not batch effects
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
#setwd('/Users/cassia/Desktop')
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")
a$subcluster <- paste(a$cluster, a$subclusterNum, sep=".")

# get rid of all columns except the cluster and class
keeps <- c("subcluster", "prefix")
b <- a[keeps]

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$prefix, b$subcluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("subcluster.replicate.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap(t(ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        column_order = order(as.numeric(rownames(ccmatrixpercent))), 
        row_order = get_order(o))
dev.off()

Heatmap((ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))

############# tissue+time vs. subcluster (order tissues)
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
#setwd('/Users/cassia/Desktop')
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")
a$subcluster <- paste(a$cluster, a$subclusterNum, sep=".")

# get rid of all columns except the cluster and class
keeps <- c("subcluster", "class")
b <- a[keeps]

b$class <- factor(b$class, c("CN3_CN4-10.5", "CN3_CN4-11.5", "CN6_10.5", "CN6-11.5", "CN7-10.5", "CN7-11.5", "CN12-10.5", "CN12-11.5", "SMN-10.5", "SMN-11.5", "CN3_CN4-neg-10.5", "CN3_CN4-neg-11.5", "CN6-neg-10.5", "CN6-neg-11.5", "CN7-neg-10.5", "CN7-neg-11.5", "CN12-neg-10.5", "CN12-neg-11.5", "SMN-neg-10.5", "SMN-neg-11.5")) 

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$class, b$subcluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("subcluster.tissue.time.gfp.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap(ccmatrixpercent, 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
dev.off()

############# tissue+time vs. subcluster (order tissues), flip axes
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
#setwd('/Users/cassia/Desktop')
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")
a$subcluster <- paste(a$cluster, a$subclusterNum, sep=".")

# get rid of all columns except the cluster and class
keeps <- c("subcluster", "class")
b <- a[keeps]

b$class <- factor(b$class, c("CN3_CN4-10.5", "CN3_CN4-11.5", "CN6-10.5", "CN6-11.5", "CN7-10.5", "CN7-11.5", "CN12-10.5", "CN12-11.5", "SMN-10.5", "SMN-11.5", "CN3_CN4-neg-10.5", "CN3_CN4-neg-11.5", "CN6-neg-10.5", "CN6-neg-11.5", "CN7-neg-10.5", "CN7-neg-11.5", "CN12-neg-10.5", "CN12-neg-11.5", "SMN-neg-10.5", "SMN-neg-11.5")) 

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$class, b$subcluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

pdf("subcluster.tissue.time.gfp.pdf", useDingbats=F)
# generate heatmap
library(ComplexHeatmap)
Heatmap((ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
dev.off()


############# time vs. cluster
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
#setwd('/Users/cassia/Desktop')
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")
a$subcluster <- paste(a$cluster, a$subclusterNum, sep=".")

# get rid of all columns except the cluster and class
keeps <- c("cluster", "time")
b <- a[keeps]

b$class <- factor(b$class, c("CN3_CN4-10.5", "CN3_CN4-11.5", "CN6-10.5", "CN6-11.5", "CN7-10.5", "CN7-11.5", "CN12-10.5", "CN12-11.5", "SMN-10.5", "SMN-11.5", "CN3_CN4-neg-10.5", "CN3_CN4-neg-11.5", "CN6-neg-10.5", "CN6-neg-11.5", "CN7-neg-10.5", "CN7-neg-11.5", "CN12-neg-10.5", "CN12-neg-11.5", "SMN-neg-10.5", "SMN-neg-11.5")) 

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$time, b$cluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

# generate heatmap
library(ComplexHeatmap)
Heatmap((ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))

############# time vs. subcluster
#Generate a heatmap of cluster vs tissue+gfp+time
#all cells are represented

#install the complexheatmap package from github
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")

#read data table
library(dplyr)
a <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
a$class <- paste(a$tissue, a$time, sep="-")
a$subcluster <- paste(a$cluster, a$subclusterNum, sep=".")

# get rid of all columns except the cluster and class
keeps <- c("subcluster", "time")
b <- a[keeps]

b$class <- factor(b$class, c("CN3_CN4-10.5", "CN3_CN4-11.5", "CN6-10.5", "CN6-11.5", "CN7-10.5", "CN7-11.5", "CN12-10.5", "CN12-11.5", "SMN-10.5", "SMN-11.5", "CN3_CN4-neg-10.5", "CN3_CN4-neg-11.5", "CN6-neg-10.5", "CN6-neg-11.5", "CN7-neg-10.5", "CN7-neg-11.5", "CN12-neg-10.5", "CN12-neg-11.5", "SMN-neg-10.5", "SMN-neg-11.5")) 

# put data from b into a matrix of percentages 
# % of cluster that is each class
clusterclasstable <- table(b$time, b$subcluster) #table with number in each cluster/class
cctablepercent <- prop.table(clusterclasstable, 2) #table of percentages of columns
ccmatrixpercent <- unclass(cctablepercent) #turn table into matrix for heatmap functions

# choose color scheme
library(circlize)
col_whiblu = colorRamp2(c(0, 0.5, 1), c("white", "cornflowerblue", "blue"))

# order the clusters (columns)
library(seriation)
o = seriate(dist(t(ccmatrixpercent)), method = "GW")

# generate heatmap
library(ComplexHeatmap)
Heatmap((ccmatrixpercent), 
        name="all cells", 
        col = col_whiblu,
        border=TRUE, rect_gp = gpar(col = "azure4", lwd = 1),
        row_title = "tissue and time", column_title = "subcluster",
        row_order = order(as.numeric(rownames(ccmatrixpercent))), 
        column_order = get_order(o))
