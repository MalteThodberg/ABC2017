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

## ------------------------------------------------------------------------
res <- as.data.frame(topTags(gallein, n=Inf, sort.by="none"))
ggplot(res, aes(x=logFC, y=-log10(PValue), color=FDR < 0.05)) + geom_point()

## ------------------------------------------------------------------------
fitQL <- glmQLFit(y=disp_zebra, design=mod, robust=TRUE)
galleinQL <- glmQLFTest(fitQL, coef="galleintreated")
summary(decideTestsDGE(galleinQL, p.value=0.1))

## ------------------------------------------------------------------------
plotQLDisp(fitQL)

## ------------------------------------------------------------------------
# Load the data
data("pasilla")

# Set a datasets
dataset <- pasilla

# Experiment with model matrices
mod <- model.matrix(~condition, data=dataset$Design)

# Look how the experiment is designed
dataset$Design

## ------------------------------------------------------------------------
# Trim 
EM <- subset(dataset$Expression, rowSums(dataset$Expression >= 5) >= 3)

# Calculate normalization factors
dge <- calcNormFactors(DGEList(EM), method="TMM")

# Calculate dispersion
disp <- estimateDisp(dge, design=mod, robust=TRUE)

# Fit models
fit <- glmFit(disp, design=mod)

# Perform the test of the given coefficient
res <- glmLRT(fit, coef="conditionuntreated")

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  # Dispersion
#  plotBCV(disp)
#  
#  # topTages
#  topTags(res)
#  
#  # Up and down regulated genes
#  summary(decideTestsDGE(res))
#  
#  # Smear plot
#  plotSmear(res, de.tags=rownames(res)[decideTestsDGE(res) != 0])
#  
#  # Volcano
#  all_genes <- as.data.frame(topTags(res, n=Inf, sort.by="none"))
#  ggplot(all_genes, aes(x=logFC, y=-log10(PValue), color=FDR < 0.05)) + geom_point()

