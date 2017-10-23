## ------------------------------------------------------------------------
# Load the data
library(ABC2017)
library(ggplot2)
theme_set(theme_minimal())
library(edgeR)
data("zebrafish")

## ------------------------------------------------------------------------
# Trim 
EM_zebra <- subset(zebrafish$Expression, rowSums(zebrafish$Expression >= 5) >= 3)

# Calculate normalization factors
dge_zebra <- DGEList(EM_zebra)
dge_zebra <- calcNormFactors(dge_zebra, method="TMM")

## ------------------------------------------------------------------------
mod <- model.matrix(~gallein, data=zebrafish$Design)
mod

## ------------------------------------------------------------------------
disp_zebra <- estimateDisp(dge_zebra, design=mod, robust=TRUE)

## ------------------------------------------------------------------------
plotBCV(y=disp_zebra)

## ------------------------------------------------------------------------
fit_zebra <- glmFit(disp_zebra, design=mod)

## ------------------------------------------------------------------------
# Using the name of the coefficient
gallein <- glmLRT(fit_zebra, coef="galleintreated")

# Using the index of the coefficient
gallein <- glmLRT(fit_zebra, coef=2)

## ------------------------------------------------------------------------
topTags(gallein)

## ------------------------------------------------------------------------
is_de <- decideTestsDGE(gallein, p.value=0.1)
head(is_de)
summary(is_de)

## ------------------------------------------------------------------------
plotSmear(gallein, de.tags=rownames(gallein)[is_de != 0])

## ---- tidy=TRUE----------------------------------------------------------
plot(-log10(PValue)~logFC, data=topTags(gallein, n=Inf, sort.by = "none"), col=factor(decideTestsDGE(gallein) != 0), pch=16)

## ------------------------------------------------------------------------
fitQL <- glmQLFit(y=disp_zebra, design=mod, robust=TRUE)
galleinQL <- glmQLFTest(fitQL, coef="galleintreated")
summary(decideTestsDGE(galleinQL, p.value=0.1))

## ------------------------------------------------------------------------
plotQLDisp(fitQL)

## ------------------------------------------------------------------------
# Load the data
data("pasilla")

# This will makes "treated the intercept!"
mod <- model.matrix(~condition, data=pasilla$Design)
mod

## ------------------------------------------------------------------------
# Relevel to make "untreated" intercept
pasilla$Design$condition <- relevel(pasilla$Design$condition, "untreated")
mod <- model.matrix(~condition, data=pasilla$Design)
mod

## ------------------------------------------------------------------------
# Relevel to make "untreated" intercept
pasilla$Design$condition <- relevel(pasilla$Design$condition, "untreated")
mod_batch <- model.matrix(~condition+type, data=pasilla$Design)
mod_batch

## ------------------------------------------------------------------------
# Trim 
EM <- subset(pasilla$Expression, rowSums(pasilla$Expression >= 5) >= 3)

# Calculate normalization factors
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Calculate dispersion
disp <- estimateDisp(dge, design=mod, robust=TRUE)

# Fit models
fit <- glmFit(disp, design=mod)

# Perform the test of the given coefficient
res <- glmLRT(fit, coef="conditiontreated")

## ------------------------------------------------------------------------
# Calculate dispersion
disp_batch <- estimateDisp(dge, design=mod_batch, robust=TRUE)

# Fit models
fit_batch <- glmFit(disp_batch, design=mod_batch)

# Perform the test of the given coefficient
res_batch <- glmLRT(fit_batch, coef="conditiontreated")

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
plotBCV(disp)
plotBCV(disp_batch)

## ------------------------------------------------------------------------
summary(decideTestsDGE(res))
summary(decideTestsDGE(res_batch))

## ------------------------------------------------------------------------
topTags(res)
topTags(res_batch)

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
plotSmear(res_batch, de.tags=rownames(res)[decideTestsDGE(res) != 0])
plotSmear(res_batch, de.tags=rownames(res_batch)[decideTestsDGE(res_batch) != 0])

## ---- tidy=TRUE----------------------------------------------------------
par(mfrow=c(1,2))
plot(-log10(PValue)~logFC, data=topTags(res, n=Inf, sort.by = "none"), col=factor(decideTestsDGE(res) != 0), pch=16)
plot(-log10(PValue)~logFC, data=topTags(res_batch, n=Inf, sort.by = "none"), col=factor(decideTestsDGE(res_batch) != 0), pch=16)

## ------------------------------------------------------------------------
fitQL_batch <- glmQLFit(y=disp_batch, design=mod_batch, robust=TRUE)
testQL_batch <- glmQLFTest(fitQL_batch, coef=2)
summary(decideTestsDGE(testQL_batch, p.value=0.05))

## ------------------------------------------------------------------------
plotQLDisp(fitQL)

