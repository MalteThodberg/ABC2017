## ------------------------------------------------------------------------
# Libraries
library(ABC2017)
library(ggplot2)
theme_set(theme_minimal())
library(edgeR)

# data
data("rat")

## ------------------------------------------------------------------------
rat$Design

## ------------------------------------------------------------------------
# Trim 
EM <- subset(rat$Expression, rowSums(rat$Expression >= 3) >= 2)

# Calculate normalization factors
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Normalized values
logEM <- cpm(dge, prior.count=1, log=TRUE)

## ------------------------------------------------------------------------
plotDensities(logEM, group=rat$Design$Time)

## ------------------------------------------------------------------------
# Perform MDS, but only save output
mds <-plotMDS(dge, top=1000, gene.selection="pairwise", plot=FALSE)

# Save as a data.frame
P <- data.frame(rat$Design, mds$cmdscale.out)

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, color=protocol, shape=Time, label=rownames(P))

## ------------------------------------------------------------------------
# Relevel to make "untreated" intercept
rat$Design$Time <- relevel(factor(rat$Design$Time), "2 weeks")
mod <- model.matrix(~protocol*Time, data=rat$Design)
mod

## ------------------------------------------------------------------------
# Perform MDS, but only save output
ldf <-plotRLDF(logEM, design=mod, nprobes=1000, plot=FALSE, trend=TRUE, robust=TRUE)

# Save as a data.frame
P <- data.frame(rat$Design, ldf$training)

## ------------------------------------------------------------------------
qplot(data=P, x=X1, y=X2, color=protocol, shape=Time, label=rownames(P))

## ------------------------------------------------------------------------
# Calculate dispersion
disp <- estimateDisp(dge, design=mod, robust=TRUE)

# Fit models
fit <- glmFit(disp, design=mod)

# Perform the tests for the coefficients
time_effect <- glmLRT(fit, coef="Time2 months")
SNL_effect <- glmLRT(fit, coef="protocolL5 SNL")
interaction_effect <- glmLRT(fit, coef="protocolL5 SNL:Time2 months")

## ------------------------------------------------------------------------
plotBCV(disp)

## ------------------------------------------------------------------------
nDE <- lapply(list(time=time_effect, SNL=SNL_effect, interaction=interaction_effect), decideTestsDGE)
nDE <- sapply(nDE, summary)
rownames(nDE) <- c("Down", "Stable", "Up")
nDE

## ------------------------------------------------------------------------
plotSmear(SNL_effect, de.tags=rownames(SNL_effect)[decideTestsDGE(SNL_effect) != 0])

## ------------------------------------------------------------------------
topTags(SNL_effect)

