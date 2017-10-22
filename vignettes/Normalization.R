## ------------------------------------------------------------------------
sample1 = c(10, 20, 30, 10, 10, 10) # Library size of 100 counts
sample2 = 2 + sample1 * 2 # Double library size
sample3 = 1 + sample1 * 3 # Triple library size
EM = data.frame(sample1, sample2, sample3)

EM

## ------------------------------------------------------------------------
colSums(EM)

## ------------------------------------------------------------------------
scale(EM, center=FALSE, scale=colSums(EM)) # Lets forget the M-part for now...

## ------------------------------------------------------------------------
EM.DE = EM
EM.DE[4:6,2] = EM.DE[4:6,2] * 5
EM.DE[4:6,3] = EM.DE[4:6,3] * 4

EM.DE

## ------------------------------------------------------------------------
scale(EM.DE, center=FALSE, scale=colSums(EM.DE))

## ------------------------------------------------------------------------
dist(t(scale(EM, center=FALSE, scale=colSums(EM))))
dist(t(scale(EM.DE, center=FALSE, scale=colSums(EM.DE))))

## ------------------------------------------------------------------------
library(ABC2017)
library(edgeR)
library(ggplot2)
theme_set(theme_minimal()) # Make ggplots prettier

## ------------------------------------------------------------------------
data(zebrafish)

## ------------------------------------------------------------------------
head(zebrafish$Expression)

## ------------------------------------------------------------------------
head(zebrafish$Design)

## ------------------------------------------------------------------------
plotDensities(zebrafish$Expression, legend="topright")

## ------------------------------------------------------------------------
# Pseoducount and log
plotDensities(log(zebrafish$Expression+1), legend="topright")

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(zebrafish$Expression > 1)

trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Pseoducount and log
log_trimmed_em <-  log(trimmed_em + 1)

## ------------------------------------------------------------------------
plotDensities(log_trimmed_em, legend="topright")

## ------------------------------------------------------------------------
# Create DGEList-object from the trimmed em
dge <- DGEList(trimmed_em)

# Use edgeR to calculate normalization factors
dge <- calcNormFactors(object=dge, method="TMM")

# calculate log cpm values
TMM_em <- cpm(x=dge, log=TRUE, prior.count=1.0)

head(TMM_em)

## ------------------------------------------------------------------------
plotDensities(TMM_em)

## ------------------------------------------------------------------------
qplot(data=log_trimmed_em, x=Ctl3, y=Trt13, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_em), x=Ctl3, y=Trt13, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
# Convert to a DGElist
dge <- DGEList(trimmed_em)

# Normalize using each of four methods
edgeR_methods <- c("none", "TMM", "RLE", "upperquartile")
dges <- lapply(edgeR_methods, calcNormFactors, object=dge)

# Calculate CPMs
norms <- lapply(dges, cpm, log=TRUE)

## ------------------------------------------------------------------------
par(mfrow=c(2,2))
mapply(plotDensities, norms, edgeR_methods, MoreArgs=list(group=NULL, col=NULL, legend=FALSE))

## ------------------------------------------------------------------------
# Extract the normalization factors
norm_factors <- sapply(dges, function(x) x$samples$norm.factors)
colnames(norm_factors) <- edgeR_methods

## ------------------------------------------------------------------------
plot(as.data.frame(norm_factors))

## ------------------------------------------------------------------------
# Create DGEList-object from the trimmed em
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")

# calculate log cpm values
TMM_v <- cpm(x=dge, log=TRUE, prior.count=0.1)
TMM_w <- cpm(x=dge, log=TRUE, prior.count=1.0)
TMM_x <- cpm(x=dge, log=TRUE, prior.count=5.0)
TMM_y <- cpm(x=dge, log=TRUE, prior.count=10.0)
TMM_z <- cpm(x=dge, log=TRUE, prior.count=20.0)

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_v), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_w), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_y), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_x), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

## ------------------------------------------------------------------------
qplot(data=as.data.frame(TMM_z), x=Ctl3, y=Ctl5, alpha=I(0.1)) + geom_smooth(method="gam") + geom_abline(color="red")

