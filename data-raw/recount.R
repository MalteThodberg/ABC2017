# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
# 	* Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

library(forcats)

#### Rat neuron ####

# Get data
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData")
load(file=con)
close(con)
bot <- hammer.eset

RatNeurons <- list(Design=pData(bot),
									 Annotation=fData(bot),
									 Expression=as.matrix(exprs(bot)))

names(RatNeurons$Expression) <- NULL

save(RatNeurons, file="./data/rat.rda")
