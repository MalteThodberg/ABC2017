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

rat <- list(Design=pData(bot),
									 Annotation=fData(bot),
									 Expression=as.matrix(exprs(bot)))

names(rat$Expression) <- NULL

devtools::use_data(rat, overwrite = TRUE)

#### Asmann Tissues ####

# Get data
con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bot = bodymap.eset

tissues <- list(Design=pData(bot),
											Annotation=fData(bot),
											Expression=as.matrix(exprs(bot)))

names(tissues$Expression) <- NULL

devtools::use_data(tissues, overwrite = TRUE)
