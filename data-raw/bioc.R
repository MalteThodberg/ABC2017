# * Creating `data-raw`.
# * Adding `data-raw` to `.Rbuildignore`.
# Next:
# 	* Add data creation scripts in data-raw
# * Use devtools::use_data() to add data to package

library(forcats)

#### Pasilla ####
library(DESeq2)
library(pasilla)
data("pasillaGenes")

pEM <- counts(pasillaGenes)
pGR <- data.frame(geneName=rownames(pEM),
									row.names=rownames(pEM))

pDM <- pData(pasillaGenes)
pDM <- subset(pDM, select=c(condition, type))

pasilla <- list(Expression=pEM, Design=pDM, Annotation=pGR)

devtools::use_data(pasilla, overwrite=TRUE)

#### Zebrafish ####
library(zebrafishRNASeq)
data(zfGenes)

zEM <- zfGenes
zDM <- data.frame(gallein=factor(x=c("control", "control", "control",
																		 "treated", "treated", "treated"),
																 levels=c("control", "treated")),
									row.names=colnames(zEM))
zGR <- data.frame(id=rownames(zEM), row.names=rownames(zEM))

zebrafish <- list(Expression=zEM, Design=zDM, Annotation=zGR)
devtools::use_data(zebrafish, overwrite=TRUE)

#### Fission yeast ####
library(fission)
data(fission)

# Subset to simpel interaction
small <- subset(fission, select=minute %in% c(0, 30))

fEM <- assay(small, "counts")
fGR <- as.data.frame(rowRanges(small))

fDM <- as.data.frame(colData(small))
fDM$minute <- fct_drop(fDM$minute)

yeast <- list(Expression=fEM, Design=fDM, Annotation=fGR)
devtools::use_data(yeast, overwrite=TRUE)


