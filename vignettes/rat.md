---
title: "Analysis of the rat dataset"
author: "Malte Thodberg"
date: "2017-10-23"
output:
  ioslides_presentation:
    smaller: true
    highlight: tango
    transition: faster
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{rat}
  %\usepackage[UTF-8]{inputenc}
---

## Libraries and data


```r
# Libraries
library(ABC2017)
library(ggplot2)
theme_set(theme_minimal())
library(edgeR)
```

```
## Loading required package: limma
```

```r
# data
data("rat")
```

## Design


```r
rat$Design
```

```
##                   sample.id num.tech.reps protocol         strain     Time
## SRX020102         SRX020102             1  control Sprague Dawley 2 months
## SRX020103         SRX020103             2  control Sprague Dawley 2 months
## SRX020104         SRX020104             1   L5 SNL Sprague Dawley 2 months
## SRX020105         SRX020105             2   L5 SNL Sprague Dawley 2 months
## SRX020091-3     SRX020091-3             1  control Sprague Dawley  2 weeks
## SRX020088-90   SRX020088-90             2  control Sprague Dawley  2 weeks
## SRX020094-7     SRX020094-7             1   L5 SNL Sprague Dawley  2 weeks
## SRX020098-101 SRX020098-101             2   L5 SNL Sprague Dawley  2 weeks
```

## Trim an normalize the data


```r
# Trim 
EM <- subset(rat$Expression, rowSums(rat$Expression >= 3) >= 2)

# Calculate normalization factors
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Normalized values
logEM <- cpm(dge, prior.count=1, log=TRUE)
```

## Check normalization


```r
plotDensities(logEM, group=rat$Design$Time)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

## MDS


```r
# Perform MDS, but only save output
mds <-plotMDS(dge, top=1000, gene.selection="pairwise", plot=FALSE)

# Save as a data.frame
P <- data.frame(rat$Design, mds$cmdscale.out)
```

## MDS


```r
qplot(data=P, x=X1, y=X2, color=protocol, shape=Time, label=rownames(P))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

## Model matrix


```r
# Relevel to make "untreated" intercept
rat$Design$Time <- relevel(factor(rat$Design$Time), "2 weeks")
mod <- model.matrix(~protocol*Time, data=rat$Design)
mod
```

```
##               (Intercept) protocolL5 SNL Time2 months
## SRX020102               1              0            1
## SRX020103               1              0            1
## SRX020104               1              1            1
## SRX020105               1              1            1
## SRX020091-3             1              0            0
## SRX020088-90            1              0            0
## SRX020094-7             1              1            0
## SRX020098-101           1              1            0
##               protocolL5 SNL:Time2 months
## SRX020102                               0
## SRX020103                               0
## SRX020104                               1
## SRX020105                               1
## SRX020091-3                             0
## SRX020088-90                            0
## SRX020094-7                             0
## SRX020098-101                           0
## attr(,"assign")
## [1] 0 1 2 3
## attr(,"contrasts")
## attr(,"contrasts")$protocol
## [1] "contr.treatment"
## 
## attr(,"contrasts")$Time
## [1] "contr.treatment"
```

## LDA with model matrix


```r
# Perform MDS, but only save output
ldf <-plotRLDF(logEM, design=mod, nprobes=1000, plot=FALSE, trend=TRUE, robust=TRUE)

# Save as a data.frame
P <- data.frame(rat$Design, ldf$training)
```

## LDA with model matrix


```r
qplot(data=P, x=X1, y=X2, color=protocol, shape=Time, label=rownames(P))
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

## Fitting models and testing


```r
# Calculate dispersion
disp <- estimateDisp(dge, design=mod, robust=TRUE)

# Fit models
fit <- glmFit(disp, design=mod)

# Perform the tests for the coefficients
time_effect <- glmLRT(fit, coef="Time2 months")
SNL_effect <- glmLRT(fit, coef="protocolL5 SNL")
interaction_effect <- glmLRT(fit, coef="protocolL5 SNL:Time2 months")
```

## BCV


```r
plotBCV(disp)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

## Number of hits


```r
nDE <- lapply(list(time=time_effect, SNL=SNL_effect, interaction=interaction_effect), decideTestsDGE)
nDE <- sapply(nDE, summary)
rownames(nDE) <- c("Down", "Stable", "Up")
nDE
```

```
##         time  SNL interaction
## Down     316 3517          15
## Stable 15401 8469       15868
## Up       182 3913          16
```

## MA-plot


```r
plotSmear(SNL_effect, de.tags=rownames(SNL_effect)[decideTestsDGE(SNL_effect) != 0])
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

## Top genes


```r
topTags(SNL_effect)
```

```
## Coefficient:  protocolL5 SNL 
##                        logFC   logCPM       LR        PValue           FDR
## ENSRNOG00000004805  4.119781 6.855273 874.6111 3.244438e-192 5.158332e-188
## ENSRNOG00000020136  5.600370 4.722513 663.7008 2.341540e-146 1.861408e-142
## ENSRNOG00000018808  3.983517 6.184542 599.6148 2.030450e-132 1.076071e-128
## ENSRNOG00000001476  4.676480 4.630042 583.2013 7.546380e-129 2.999497e-125
## ENSRNOG00000015156  3.487241 6.911991 484.7044 2.023453e-107 6.434175e-104
## ENSRNOG00000004874  2.609712 7.009175 476.2087 1.428032e-105 3.784047e-102
## ENSRNOG00000004731 -2.528212 6.453337 400.1812  5.029081e-89  1.142248e-85
## ENSRNOG00000013496  4.444739 3.812216 397.7378  1.711583e-88  3.401557e-85
## ENSRNOG00000014327  3.407460 4.828255 393.3446  1.547973e-87  2.734580e-84
## ENSRNOG00000021029  6.067354 2.998231 388.3760  1.868202e-86  2.970254e-83
```


