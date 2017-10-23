---
title: "Normalization and RNA-composition"
author: "Malte Thodberg"
date: "2017-10-23"
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

This presentations presents some background on EM-normalization for library size and RNA-composition,
as wells as some examples on how this is applied in R using the package `edgeR`.

Density curves and log-log plots will be used to explore the effects of different normalization methods.

# RNA-composition and DE

## TPM normalization
Setup simple EM:

```r
sample1 = c(10, 20, 30, 10, 10, 10) # Library size of 100 counts
sample2 = 2 + sample1 * 2 # Double library size
sample3 = 1 + sample1 * 3 # Triple library size
EM = data.frame(sample1, sample2, sample3)

EM
```

```
##   sample1 sample2 sample3
## 1      10      22      31
## 2      20      42      61
## 3      30      62      91
## 4      10      22      31
## 5      10      22      31
## 6      10      22      31
```

Note the different library sizes:

```r
colSums(EM)
```

```
## sample1 sample2 sample3 
##      90     192     276
```

## TPM normalization
TPM scaling:

```r
scale(EM, center=FALSE, scale=colSums(EM)) # Lets forget the M-part for now...
```

```
##        sample1   sample2   sample3
## [1,] 0.1111111 0.1145833 0.1123188
## [2,] 0.2222222 0.2187500 0.2210145
## [3,] 0.3333333 0.3229167 0.3297101
## [4,] 0.1111111 0.1145833 0.1123188
## [5,] 0.1111111 0.1145833 0.1123188
## [6,] 0.1111111 0.1145833 0.1123188
## attr(,"scaled:scale")
## sample1 sample2 sample3 
##      90     192     276
```

Samples can now be compared directly for analysis!

## CPM normalization
Introduce DE for some TCs

```r
EM.DE = EM
EM.DE[4:6,2] = EM.DE[4:6,2] * 5
EM.DE[4:6,3] = EM.DE[4:6,3] * 4

EM.DE
```

```
##   sample1 sample2 sample3
## 1      10      22      31
## 2      20      42      61
## 3      30      62      91
## 4      10     110     124
## 5      10     110     124
## 6      10     110     124
```
The total RNA content of sample2+3 has increased!

## CPM normalization
TPM scaling

```r
scale(EM.DE, center=FALSE, scale=colSums(EM.DE))
```

```
##        sample1    sample2    sample3
## [1,] 0.1111111 0.04824561 0.05585586
## [2,] 0.2222222 0.09210526 0.10990991
## [3,] 0.3333333 0.13596491 0.16396396
## [4,] 0.1111111 0.24122807 0.22342342
## [5,] 0.1111111 0.24122807 0.22342342
## [6,] 0.1111111 0.24122807 0.22342342
## attr(,"scaled:scale")
## sample1 sample2 sample3 
##      90     456     555
```
Non-DE genes are now under-sampled!

## CPM normalization
This can affect downstream analysis i.e. distance matrix calculations.

```r
dist(t(scale(EM, center=FALSE, scale=colSums(EM))))
```

```
##             sample1     sample2
## sample2 0.012991866            
## sample3 0.004518910 0.008472956
```

```r
dist(t(scale(EM.DE, center=FALSE, scale=colSums(EM.DE))))
```

```
##            sample1    sample2
## sample2 0.33260796           
## sample3 0.28669731 0.04593348
```

# Advanced normalization

## Setup

Packages needed for the analysis:

```r
library(ABC2017)
library(edgeR)
library(ggplot2)
theme_set(theme_minimal()) # Make ggplots prettier
```

We will use the small `zebrafish` dataset:

```r
data(zebrafish)
```

The dataset is a list which contains:

- Expression: EM
- Annotation: Annotation
- Design: Information about the samples.

The same format is used for the remaining datasets in the `ABC2017` package

## Setup
EM:

```r
head(zebrafish$Expression)
```

```
##                    Ctl1 Ctl3 Ctl5 Trt9 Trt11 Trt13
## ENSDARG00000000001  304  129  339  102    16   617
## ENSDARG00000000002  605  637  406   82   230  1245
## ENSDARG00000000018  391  235  217  554   451   565
## ENSDARG00000000019 2979 4729 7002 7309  9395  3349
## ENSDARG00000000068   89  356   41  149    45    44
## ENSDARG00000000069  312  184  844  269   513   243
```

Annotation:

```r
head(zebrafish$Design)
```

```
##       gallein
## Ctl1  control
## Ctl3  control
## Ctl5  control
## Trt9  treated
## Trt11 treated
## Trt13 treated
```

## Plotting distributions

edgeR (via limma) provides the `plotDensities` function for exploring the effect of normalization 


```r
plotDensities(zebrafish$Expression, legend="topright")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

## Plotting distributions

That did not look to good! Since the data spans multiple orders of magnitude, we can try with	 a log-scale instead.

This however brings up the problem of 0 counts - for which `log` is not defined.

The get around this probelm a small pseudo-count can be added to all counts in the EM. This does not necesarily have to be an integer, and is usually chosen to be between 0.1 and 1.0.

## Plotting distributions


```r
# Pseoducount and log
plotDensities(log(zebrafish$Expression+1), legend="topright")
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

## Plotting distributions

Notice how the lower quartile is zero - this means that we have a large number of genes with very low counts.

Counts with 1-3 counts are not very interesting, since they are likely to be either noise or expressed at biologically irrelevant levels. It's customary to perform som ad-hoc trimming or filtering to remove these prior to analysis.

Here we only keep genes with at least 2 counts in at least 4 samples:


```r
# Trim
above_one <- rowSums(zebrafish$Expression > 1)

trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Pseoducount and log
log_trimmed_em <-  log(trimmed_em + 1)
```

## Plotting distributions


```r
plotDensities(log_trimmed_em, legend="topright")
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

## Plotting distributions

Now we have a clearer picture of the distribution of counts within each sample. The large difference in distributions shows the need for normalization, before the samples can be compared.

As with everything in R, we do not have to recode everything from scratch. The edgeR package has a function `cpm` which has implented a large number of normalization methods and log-transformation.

edgeR does this by implementing the use of normalization factors, which is use to rescale the actual library sizes to take into account differences in RNA-composition.

## Using edgeR to normalize

Using edgeR is simple, but first we must save the EM as a `DGEList`:

```r
# Create DGEList-object from the trimmed em
dge <- DGEList(trimmed_em)

# Use edgeR to calculate normalization factors
dge <- calcNormFactors(object=dge, method="TMM")

# calculate log cpm values
TMM_em <- cpm(x=dge, log=TRUE, prior.count=1.0)

head(TMM_em)
```

```
##                        Ctl1     Ctl3       Ctl5     Trt9      Trt11
## ENSDARG00000000001 2.863841 1.070888  2.9026167 1.876326 -0.4314598
## ENSDARG00000000002 3.854225 3.360984  3.1620111 1.563859  3.3678160
## ENSDARG00000000018 3.225830 1.928374  2.2617491 4.309565  4.3375880
## ENSDARG00000000019 6.152048 6.250098  7.2664067 8.029591  8.7165741
## ENSDARG00000000068 1.103626 2.524350 -0.1101648 2.419952  1.0285423
## ENSDARG00000000069 2.901188 1.578056  4.2156712 3.269221  4.5232027
##                       Trt13
## ENSDARG00000000001 3.987476
## ENSDARG00000000002 4.999126
## ENSDARG00000000018 3.860666
## ENSDARG00000000019 6.425998
## ENSDARG00000000068 0.207237
## ENSDARG00000000069 2.646674
```

## Using edgeR to normalize
The resulting plot shows a nicer alignment of the main peak:

```r
plotDensities(TMM_em)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

## log-log plots

Another way of visualizing normalization is via a log-log plot. This is simply a scatterplot with paired expression values for two samples.

Although it only allows for pairwise comparison, it is a nice way to see the effect of normalization and the variance of expression at different levels.

## log-log plots
First we consider the (trimmed) log(counts+1)

```r
qplot(data=log_trimmed_em, x=Ctl3, y=Trt13, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

## log-log plots
Compare this with edgeR's TMM normalization:

```r
qplot(data=as.data.frame(TMM_em), x=Ctl3, y=Trt13, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)

## Exercises

Team up and do the following exercises:

1. Make a plot of density curves showing the effect of all the different normalization methods implemented by `edgeR`
2. Make a plot that compares the estimated normalization factors for the `TMM` and `RLE` normalization methods.
3. Investigate the effect of different pseudocount values using log-log plots.

*HINT for 1:* Use `apply`-family and facets to compare multiple datasets. Remember to include `method="none"`.

*HINT for 2:* Read the `calcNormFactors` help file to see where the normalization factors are stored

## Question 1 code

```r
# Convert to a DGElist
dge <- DGEList(trimmed_em)

# Normalize using each of four methods
edgeR_methods <- c("none", "TMM", "RLE", "upperquartile")
dges <- lapply(edgeR_methods, calcNormFactors, object=dge)

# Calculate CPMs
norms <- lapply(dges, cpm, log=TRUE)
```

## Question 1 plot


```r
par(mfrow=c(2,2))
mapply(plotDensities, norms, edgeR_methods, MoreArgs=list(group=NULL, col=NULL, legend=FALSE))
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

```
##   [,1]         [,2]         [,3]         [,4]        
## X Numeric,3072 Numeric,3072 Numeric,3072 Numeric,3072
## Y Numeric,3072 Numeric,3072 Numeric,3072 Numeric,3072
```

## Question 2 code

```r
# Extract the normalization factors
norm_factors <- sapply(dges, function(x) x$samples$norm.factors)
colnames(norm_factors) <- edgeR_methods
```

## Question 2 plot

```r
plot(as.data.frame(norm_factors))
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

## Question 3 code

```r
# Create DGEList-object from the trimmed em
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")

# calculate log cpm values
TMM_v <- cpm(x=dge, log=TRUE, prior.count=0.1)
TMM_w <- cpm(x=dge, log=TRUE, prior.count=1.0)
TMM_x <- cpm(x=dge, log=TRUE, prior.count=5.0)
TMM_y <- cpm(x=dge, log=TRUE, prior.count=10.0)
TMM_z <- cpm(x=dge, log=TRUE, prior.count=20.0)
```

## Question 3 plot v

```r
qplot(data=as.data.frame(TMM_v), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)

## Question 3 plot w

```r
qplot(data=as.data.frame(TMM_w), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)

## Question 3 plot x

```r
qplot(data=as.data.frame(TMM_y), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-1.png)

## Question 3 plot y

```r
qplot(data=as.data.frame(TMM_x), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)

## Question 3 plot z

```r
qplot(data=as.data.frame(TMM_z), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28-1.png)
