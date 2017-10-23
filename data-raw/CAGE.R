library(CAGEtestR)
library(BiocParallel)
register(MulticoreParam(2))

### Mouse data

# Some example data
design <- read.table("~/Desktop/mm9_nanotubes/design_matrix.txt", header=TRUE, stringsAsFactors=FALSE)

# Better formated design
design$Samples <- gsub(x=design$Samples, pattern=".fastq", replacement="")
design$BigWigPlus <- paste0("mm9.", design$Samples, ".plus.bw")
design$BigWigMinus <- paste0("mm9.", design$Samples, ".minus.bw")
design$Class <- factor(design$Class)

# Format and sort
design <- subset(design, Name != "C548_2", select=-c(Samples))
design <- design[order(design$Class, design$Name),]
rownames(design) <- design$Name
design <- DataFrame(design)

# Set up BigWig files
bw_plus <- BigWigFileList(file.path("~/Desktop/mm9_nanotubes", design$BigWigPlus))
bw_minus <- BigWigFileList(file.path("~/Desktop/mm9_nanotubes", design$BigWigMinus))
names(bw_plus) <- rownames(design)
names(bw_minus) <- rownames(design)

### Run

# Test on full
CTSS <- quantifyCTSSs(bw_plus, bw_minus); gc()
colData(CTSS) <- design

# Test on subset
pooled <- subsetBySupport(CTSS, unexpressed=0); gc()
pooled <- calcTPM(pooled); gc()
pooled <- calcPooled(pooled); gc()
pooled <- rowRanges(pooled); gc()

# TCs
TCs <- clusterUnidirectionally(pooled, pooledCutoff=3); gc()
SE <- quantifyClusters(CTSS, TCs); gc()

SE <- calcTPM(SE, outputColumn="TCTags"); gc()
SE <- subsetBySupport(SE, inputAssay="TPM", unexpressed=1, minSamples=2); gc()

# Gene-level
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene
rowRanges(SE)$geneID <- CAGEfightR::assignGeneID(rowRanges(SE), txdb)
GE <- CAGEfightR::sumOverGenes(SE, inputAssay="counts", geneID="geneID")

# Prepare data

mouse <- list(Expression=assay(GE, "counts"),
							Design=as.data.frame(colData(GE)),
							Annotation=data.frame(geneID=rownames(GE), rowData(GE)))

devtools::use_data(mouse, overwrite = TRUE)
