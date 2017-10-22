---
title: "Exploratory Data Analysis (EDA)"
author: "Malte Thodberg"
date: "2017-10-20"
output:
  ioslides_presentation:
    smaller: true
    highlight: tango
    transition: faster
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{EDA}
  %\usepackage[UTF-8]{inputenc}
---

## Introduction

This presentation will provide a short introduction to different ways of visualizing a genomics dataset:

Dimensionality reductions (projections)

- PCA
- LDA and PLSDA
- MDS

Heatmap like:

- Correlograms
- Heatmaps with prior k-means clustering.

Clusterings:

- Agglomerative
- Divisive

## Setup
Packages:

```r
library(ABC2017)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(edgeR)
```

```
## Loading required package: limma
```

```r
data(zebrafish)
```

We start from the trimmed and normalized EM from before:

```r
# Trim
above_one <- rowSums(zebrafish$Expression > 1)
trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)
```

# Principal Components Anaysis (PCA)

## PCA

You already know PCA from BoHTA. PCA decomposes the (scaled) EM into principle components representing orthogonal rotations that maximizes total variance. 

This is a powerful way of visualizing data since it will highlight the major but mutually exclusive patterns, as well as quantifiying the contribution of each pattern to the total through explained variance of the components.

PCA in R is simple (Recap from BoHTA): 

```r
# Perfrom PCA
pca <- prcomp(x=t(EM), scale=TRUE, center=TRUE)

# Inspect components
summary(pca)
```

```
## Importance of components%s:
##                            PC1     PC2     PC3     PC4     PC5       PC6
## Standard deviation     82.4185 66.7112 57.2301 53.4794 48.8824 2.234e-13
## Proportion of Variance  0.3436  0.2251  0.1657  0.1447  0.1209 0.000e+00
## Cumulative Proportion   0.3436  0.5688  0.7344  0.8791  1.0000 1.000e+00
```

## PCA


```r
qplot(data=as.data.frame(pca$x), x=PC1, y=PC2, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

## PCA


```r
qplot(data=as.data.frame(pca$x), x=PC1, y=PC3, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

# Heatmaps

## Aim

Heatmaps are classic way of visulizaing datasets as a clustered matrix with color-reprensentation of magnitude.

In addition to visualization of the EM, it also features row and column wise hierachical clusterings, and can be further annotated with groupings along the margins.

Heatmaps are normally most useful for small dataset, otherwise rows, columns and trees tend to get smeared.

A way to get around this is to group genes prior to plotting, which allows for plotting of much bigger data but loosing resolution of individual features.

## Heatmap of large datasets

```r
pheatmap::pheatmap(mat=EM, annotation_col=zebrafish$Design, scale="row")
# Running this takes waaaaay to long!
```

## Heatmap of large datasets

```r
pheatmap::pheatmap(mat = EM, kmeans_k = 10, annotation_col = zebrafish$Design, 
    scale = "row")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

## Heatmap of large datasets

```r
pheatmap::pheatmap(mat = EM, kmeans_k = 100, annotation_col = zebrafish$Design, 
    scale = "row")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

## Heatmap of large datasets

```r
pheatmap::pheatmap(mat = EM, kmeans_k = 1000, annotation_col = zebrafish$Design, 
    scale = "row")
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

# Distance matrices

## Aim

In many cases we are not (yet) interested in individual genes or clusters, but only on how the samples are different.

A simple way a quantifiying this is to calculate the distance (i.e. euclidian distance) between samples, and plot these values in a heatmap.

Following is examples of doing this. Let's first define the distances:

```r
# Transpose for samples-wise distances
distmat <- dist(t(EM)) #
```

## Heatmap of distance matrix

```r
pheatmap::pheatmap(mat = as.matrix(distmat), color = brewer.pal(name = "RdPu", 
    n = 9), clustering_distance_rows = distmat, clustering_distance_cols = distmat, 
    annotation_col = zebrafish$Design, annotation_row = zebrafish$Design)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

# Multi-Dimensional Scaling (MDS)

## Aim

We can view the distance matrix as a new high-dimensional space, where measurements are now distances to other samples.

This means we can use dimensionality reduction to try an represent the distance matrix in fewer dimensions!

Multi-Dimensional Scaling finds a low-dimensional representation (usually 2D) of a high-dimensional distance matrix, preserving the distances between samples as best as possible.

Let's see how this works first by using the example data `UScitiesD`!

## MDS: Small example

The UScitiesD is `dist` object holding distances between 9 major US cities:


```r
class(UScitiesD)
```

```
## [1] "dist"
```

```r
UScitiesD
```

```
##               Atlanta Chicago Denver Houston LosAngeles Miami NewYork
## Chicago           587                                                
## Denver           1212     920                                        
## Houston           701     940    879                                 
## LosAngeles       1936    1745    831    1374                         
## Miami             604    1188   1726     968       2339              
## NewYork           748     713   1631    1420       2451  1092        
## SanFrancisco     2139    1858    949    1645        347  2594    2571
## Seattle          2182    1737   1021    1891        959  2734    2408
## Washington.DC     543     597   1494    1220       2300   923     205
##               SanFrancisco Seattle
## Chicago                           
## Denver                            
## Houston                           
## LosAngeles                        
## Miami                             
## NewYork                           
## SanFrancisco                      
## Seattle                678        
## Washington.DC         2442    2329
```

## MDS: Small example

We can use `cmdscale` function to reduce the dimensionality of the distance matrix:


```r
# Call cmdscale
mds <- cmdscale(UScitiesD, k = 2)
mds
```

```
##                     [,1]       [,2]
## Atlanta        -718.7594  142.99427
## Chicago        -382.0558 -340.83962
## Denver          481.6023  -25.28504
## Houston        -161.4663  572.76991
## LosAngeles     1203.7380  390.10029
## Miami         -1133.5271  581.90731
## NewYork       -1072.2357 -519.02423
## SanFrancisco   1420.6033  112.58920
## Seattle        1341.7225 -579.73928
## Washington.DC  -979.6220 -335.47281
```

We can then easily plot this 2D representation: How does it look?

## MDS: Small example


```r
qplot(x = mds[, 1], y = mds[, 2], label = rownames(mds), geom = "text")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

## MDS: Real example

Let's use the distance matrix from before and do MDS:


```r
mds <- cmdscale(distmat, k = 2)
```

Gather all the info we need for plotting:

```r
P <- data.frame(zebrafish$Design, mds)
```

How does the plot look compared to the PCA-plot?

## MDS: Real example


```r
qplot(x = X1, y = X2, label = rownames(P), color = gallein, geom = "text", data = P)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

## MDS: edgeR and limma

One of the advantages of MDS is that the distance matrix can be calculated in any way possible. While it might not always be possible to do PCA (i.e. missing values), it is usually possible to define some distance measure between samples.

MDS is the built-in visualization method in edgeR (and limma), where distances are calculated in a different way:

Instead of calculating distance based on all genes, it only used topN features: Either defined by overall differences or pairwise differences.

This goes back to the assumption that most genes are not changing, and therefore there's little to gain in including them in the analysis.

## MDS: edgeR and limma

`plotMDS` will do all the work for us (can also be called directly on a `DGEList`):


```r
plotMDS(EM, top = 1000)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)


## MDS: edgeR and limma

The same plot, but using far fewer features:


```r
plotMDS(EM, top = 10)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

# Supervised projections

## Aim

So far we have look at unsupervised projections: We have not taken anything known about the samples into account.

Another approach is supervised projections, where we try to create a low-dimensional representation that best captures differences between our known groups.

A popular methods for this is Linear Discriminant Analysis (LDA). Unfortunately, LDA runs into problems with genomics data due to high multi-colinearity of variables.

limma has a version of LDA for high-dimensional data that implements a few tricks to get around this issue. Again, this method uses only a subset of top genes.

It should be noted though, that supervised projections almost always produce nicely looking plots (due to the high number of input features). That means they should be interpretated with extra caution!

## LDA in limma

We will get back to what many of these settings mean!

```r
plotRLDF(EM, nprobes = 100, labels.y = zebrafish$Design$gallein, trend = TRUE, 
    robust = TRUE)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

## PLSDA

For those interested, another (more general, but much slower) alternative to LDA for genomics data is Partial Least Squares Discriminant Analysis (PLSDA):


```r
# Perfrom PLSDA
plsda <- DiscriMiner::plsDA(variables = scale(t(EM)), group = zebrafish$Design$gallein, 
    autosel = FALSE, comps = 2)

# Lots of output, we are interested in the components
summary(plsda)
```

```
##                Length Class  Mode   
## functions      39538  -none- numeric
## confusion          4  table  numeric
## scores            12  -none- numeric
## loadings       39536  -none- numeric
## y.loadings         4  -none- numeric
## classification     6  factor numeric
## error_rate         1  -none- numeric
## components        12  -none- numeric
## Q2                 6  -none- numeric
## R2                 8  -none- numeric
## VIP            59304  -none- numeric
## comp_vars      39536  -none- numeric
## comp_group         4  -none- numeric
## specs              6  -none- list
```

## PLSDA

```r
qplot(data=as.data.frame(plsda$components), x=t1, y=t2, geom=c("text"), 
			color=zebrafish$Design$gallein, label=rownames(zebrafish$Design))
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

## Exercise

The `zebrafish` dataset is very simple! It's few samples and only limited extra information.

Your exercise now is to perform EDA on the other datasets in `ABC2017`:

```r
data(pasilla)
```

If you are super quick, you can also have a look at a much more complex dataset:

```r
data(tissues)
```

You are free to do whatever you like - **but** you should be able to motivate the choice of algorithm/package/tool etc.

We will go through some example solutions at the end of today.

Happy exploring!


