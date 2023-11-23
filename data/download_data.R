library(edgeR)
library(limma)
library(org.Mm.eg.db)


## Read the counts from the downloaded data
seqdata <- read.delim("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
#
# Remove first two columns from seqdata

countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
countdata
colnames(countdata) <- substr(colnames(countdata), 1, 7)
countdata
## Calculate the Counts Per Million measure
myCPM <- cpm(countdata)
## Identify genes with at least 0.5 cpm in at least 2 samples
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
## Convert to an edgeR object
dgeObj <- DGEList(counts.keep)
## Perform TMM normalisation
dgeObj <- calcNormFactors(dgeObj)
## Obtain corrected sample information
sampleinfo <- read.delim("SampleInfo_Corrected.txt")
group <- paste(sampleinfo$CellType,sampleinfo$Status,sep=".")
group

dgeObj <- estimateCommonDisp(dgeObj)
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)


group <- as.character(group)
type <- sapply(strsplit(group, ".", fixed=T), function(x) x[1])
status <- sapply(strsplit(group, ".", fixed=T), function(x) x[2])
# Specify a design matrix with an intercept term
design <- model.matrix(~ type + status)

fit <- glmFit(dgeObj, design)
lrt.BvsL <- glmLRT(fit, coef=2)
results <- as.data.frame(topTags(lrt.BvsL,n = Inf))

results$gene <- rownames(results)

ann <- select(org.Mm.eg.db,keys=rownames(results),columns=c("CHR","ENTREZID","SYMBOL","GENENAME"))

ann$gene <- ann$ENTREZID

library(tidyverse)

results <- results %>% left_join(ann)
expression_matrix <- cpm(dgeObj,log=TRUE) %>% as.data.frame() %>%
  rownames_to_column("gene") %>% pivot_longer(cols = colnames(dgeObj), names_to = "SampleName",values_to = "expression")

save(results,expression_matrix,sampleinfo, file = "expression.RData")
