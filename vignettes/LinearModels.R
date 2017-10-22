## ------------------------------------------------------------------------
library(ABC2017)
library(ggplot2)
theme_set(theme_minimal())
data("lmExamples")
names(lmExamples)

## ------------------------------------------------------------------------
lmExamples$twoGroup1

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$twoGroup1, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
t.test(Expression~Group, data=lmExamples$twoGroup1)

## ------------------------------------------------------------------------
summary(lm(Expression~Group, data=lmExamples$twoGroup1))

## ------------------------------------------------------------------------
mod <- model.matrix(~Group, data=lmExamples$twoGroup1)

## ------------------------------------------------------------------------
mod

## ------------------------------------------------------------------------
solve(t(mod) %*% mod) %*% t(mod) %*% lmExamples$twoGroup1$Expression

## ------------------------------------------------------------------------
summary(lm(Expression~Group, data=lmExamples$twoGroup2))

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$twoGroup2, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
lmExamples$interaction1

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$interaction1, aes(x=paste0(Group, "-", Time), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
summary(lm(Expression~Group*Time, data=lmExamples$interaction1))

## ------------------------------------------------------------------------
model.matrix(Expression~Group*Time, data=lmExamples$interaction1)

## ------------------------------------------------------------------------
model.matrix(Expression~Group+Time+Group:Time, data=lmExamples$interaction1)

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$interaction2, aes(x=paste0(Group, "-", Time), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
summary(lm(Expression~Group*Time, data=lmExamples$interaction2))

## ------------------------------------------------------------------------
lmExamples$batchEffect

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$batchEffect, aes(x=paste0(Group, "-", Batch), y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
summary(lm(Expression~Group, data=lmExamples$batchEffect))

## ------------------------------------------------------------------------
summary(lm(Expression~Group+Batch, data=lmExamples$batchEffect))

## ------------------------------------------------------------------------
model.matrix(~Group+Batch, data=lmExamples$batchEffect)

## ------------------------------------------------------------------------
summary(lm(Expression~Group, data=lmExamples$threeGroup))

## ---- fig.height=4, fig.width=6------------------------------------------
ggplot(lmExamples$threeGroup, aes(x=Group, y=Expression, color=Group)) +
	geom_jitter(alpha=0.75) +
	stat_summary(fun.y="mean", fun.ymax="mean", fun.ymin="mean",
							 geom="crossbar", color="black", alpha=0.5, linetype="dashed")

## ------------------------------------------------------------------------
library(limma)
library(edgeR)
data("zebrafish")

## ------------------------------------------------------------------------
# Trim
above_one <- rowSums(zebrafish$Expression > 1)
trimmed_em <- subset(zebrafish$Expression, above_one > 3)

# Normalize
dge <- DGEList(trimmed_em)
dge <- calcNormFactors(object=dge, method="TMM")
EM <- cpm(x=dge, log=TRUE)

## ------------------------------------------------------------------------
mod <- model.matrix(~gallein, data=zebrafish$Design)
mod

## ------------------------------------------------------------------------
fit <- lmFit(EM, design=mod)
fit

## ------------------------------------------------------------------------
eb <- eBayes(fit, trend=TRUE, robust=TRUE)
plotSA(eb)

## ------------------------------------------------------------------------
tstat <- data.frame(raw=(fit$coef/fit$stdev.unscaled/fit$sigma)[,"galleintreated"],
										shrunken=eb$t[,"galleintreated"])

## ---- fig.height=4, fig.width=4------------------------------------------
ggplot(tstat, aes(x=raw, y=shrunken, color=raw-shrunken)) + 
	geom_point(alpha=0.33) +
	scale_color_distiller(palette = "Spectral") +
	geom_abline(linetype="dashed", alpha=0.75)

