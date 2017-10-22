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

