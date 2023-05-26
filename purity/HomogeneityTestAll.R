#Tests of homogeneity and completeness
#in this version, "tissue" is dependent on gfp neg/pos
#e.g. "CN12-neg" is different than "CN12"(pos)

#Test 1: cluster, Tissue+time
#Test 2: subcluster, Tissue+time
#Test 3a: cluster, Tissue+time for only GFP pos
#Test 3b: cluster, Tissue+time for only GFP neg
#Test 4a: subcluster, Tissue+time for only GFP pos
#Test 4b: subcluster, Tissue+time for only GFP neg
#Test 5a: cluster, Tissue+time for only GFP pos at e10.5
#Test 5b: cluster, Tissue+time for only GFP neg at e10.5
#Test 6a: subcluster, Tissue+time for only GFP pos at e10.5
#Test 6b: subcluster, Tissue+time for only GFP neg at e10.5
#Test 7a: cluster, Tissue+time for only GFP pos at e11.5
#Test 7b: cluster, Tissue+time for only GFP neg at e11.5
#Test 8a: subcluster, Tissue+time for only GFP pos at e11.5
#Test 8b: subcluster, Tissue+time for only GFP neg at e11.5


## install.packages("sabre")
# full original data set
d <- read.table("allSamples_finalFiltered.cellData_112019.txt", header=T, sep="\t")
#create column for "class" which is tissue-time
d$class <- paste(d$tissue, d$time, sep="-")
# create subcluster column
d$subcluster <- paste(d$cluster, d$subclusterNum, sep=".")


#filter for GFP positive or negative (in tissue name)
gfp.neg <- d[grep("-neg$", d$tissue),]
gfp.pos <- d[-grep("-neg$", d$tissue),]

#filter for embryonic day 10.5
gfp.pos.e10 <- gfp.pos[grep("-10.5$", gfp.pos$class),]
gfp.neg.e10 <- gfp.neg[grep("-10.5$", gfp.neg$class),]

#filter for embryonic day 11.5
gfp.pos.e11 <- gfp.pos[grep("-11.5$", gfp.pos$class),]
gfp.neg.e11 <- gfp.neg[grep("-11.5$", gfp.neg$class),]

## store the class values into x
x1 <- (d$class)
x2 <- (d$class)
x3a <- (gfp.pos$class)
x3b <- (gfp.neg$class)
x4a <- (gfp.pos$class)
x4b <- (gfp.neg$class)
x5a <- (gfp.pos.e10$class)
x5b <- (gfp.neg.e10$class)
x6a <- (gfp.pos.e10$class)
x6b <- (gfp.neg.e10$class)
x7a <- (gfp.pos.e11$class)
x7b <- (gfp.neg.e11$class)
x8a <- (gfp.pos.e11$class)
x8b <- (gfp.neg.e11$class)

## store the cluster values into y
y1 <- (d$cluster)
y2 <- (d$subcluster)
y3a <- (gfp.pos$cluster)
y3b <- (gfp.neg$cluster)
y4a <- (gfp.pos$subcluster)
y4b <- (gfp.neg$subcluster)
y5a <- (gfp.pos.e10$cluster)
y5b <- (gfp.neg.e10$cluster)
y6a <- (gfp.pos.e10$subcluster)
y6b <- (gfp.neg.e10$subcluster)
y7a <- (gfp.pos.e11$cluster)
y7b <- (gfp.neg.e11$cluster)
y8a <- (gfp.pos.e11$subcluster)
y8b <- (gfp.neg.e11$subcluster)


library(sabre)
#calculate vmeasure, homogeneity, and completeness
v1 <- vmeasure(x1, y1, z = NULL, B = 1)
v2 <- vmeasure(x2, y2, z = NULL, B = 1)
v3a <- vmeasure(x3a, y3a, z = NULL, B = 1)
v3b <- vmeasure(x3b, y3b, z = NULL, B = 1)
v4a <- vmeasure(x4a, y4a, z = NULL, B = 1)
v4b <- vmeasure(x4b, y4b, z = NULL, B = 1)
v5a <- vmeasure(x5a, y5a, z = NULL, B = 1)
v5b <- vmeasure(x5b, y5b, z = NULL, B = 1)
v6a <- vmeasure(x6a, y6a, z = NULL, B = 1)
v6b <- vmeasure(x6b, y6b, z = NULL, B = 1)
v7a <- vmeasure(x7a, y7a, z = NULL, B = 1)
v7b <- vmeasure(x7b, y7b, z = NULL, B = 1)
v8a <- vmeasure(x8a, y8a, z = NULL, B = 1)
v8b <- vmeasure(x8b, y8b, z = NULL, B = 1)

print("test 1")
v1
print("test 2")
v2
print("test 3a")
v3a
print("test 3b")
v3b
print("test 4a")
v4a
print("test 4b")
v4b
print("test 5a")
v5a
print("test 5b")
v5b
print("test 6a")
v6a
print("test 6b")
v6b
print("test 7a")
v7a
print("test 7b")
v7b
print("test 8a")
v8a
print("test 8b")
v8b

