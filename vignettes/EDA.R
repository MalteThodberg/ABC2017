## ------------------------------------------------------------------------
library(ABC2017)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(edgeR)

data(zebrafish)

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(zebrafish$Expression > 1)
trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

## ------------------------------------------------------------------------
# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)

## ------------------------------------------------------------------------
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))

## ------------------------------------------------------------------------
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))

## ---- eval=FALSE---------------------------------------------------------
#  pheatmap::pheatmap(mat=EM, annotation_col=zebrafish$Design, scale="row")
#  # Running this takes waaaaay to long!

## ---- tidy=TRUE----------------------------------------------------------
pheatmap::pheatmap(mat=EM, kmeans_k=10, annotation_col=zebrafish$Design, scale="row")

## ---- tidy=TRUE----------------------------------------------------------
pheatmap::pheatmap(mat=EM, kmeans_k=100, annotation_col=zebrafish$Design, scale="row")

## ---- tidy=TRUE----------------------------------------------------------
pheatmap::pheatmap(mat=EM, kmeans_k=1000, annotation_col=zebrafish$Design, scale="row")

## ------------------------------------------------------------------------
# Transpose for samples-wise distances
distmat <- dist(t(EM)) #

## ---- tidy=TRUE----------------------------------------------------------
pheatmap::pheatmap(mat=as.matrix(distmat), color=brewer.pal(name="RdPu", n=9), clustering_distance_rows=distmat, clustering_distance_cols=distmat, annotation_col=zebrafish$Design, annotation_row=zebrafish$Design)

## ---- tidy=TRUE----------------------------------------------------------
class(UScitiesD)
UScitiesD

## ---- tidy=TRUE----------------------------------------------------------
# Call cmdscale
mds <- cmdscale(UScitiesD, k=2)
mds

## ---- tidy=TRUE----------------------------------------------------------
qplot(x=mds[,1], y=mds[,2], label=rownames(mds), geom="text")

## ---- tidy=TRUE----------------------------------------------------------
mds <- cmdscale(distmat, k=2)

## ---- tidy=TRUE----------------------------------------------------------
P <- data.frame(zebrafish$Design, mds)

## ---- tidy=TRUE----------------------------------------------------------
qplot(x=X1, y=X2, label=rownames(P), color=gallein, geom="text", data=P)

## ---- tidy=TRUE----------------------------------------------------------
plotMDS(EM, top=1000)

## ---- tidy=TRUE----------------------------------------------------------
plotMDS(EM, top=10)

## ---- tidy=TRUE----------------------------------------------------------
plotRLDF(EM, nprobes=100, labels.y=zebrafish$Design$gallein, trend=TRUE, robust=TRUE)

## ---- tidy=TRUE----------------------------------------------------------
# Perfrom PLSDA
plsda <- DiscriMiner::plsDA(variables=scale(t(EM)), group=zebrafish$Design$gallein, autosel=FALSE, comps=2)

# Lots of output, we are interested in the components
summary(plsda)

## ------------------------------------------------------------------------
qplot(data=as.data.frame(plsda$components), x=t1, y=t2, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))

## ------------------------------------------------------------------------
data(pasilla)

## ------------------------------------------------------------------------
data(tissues)

