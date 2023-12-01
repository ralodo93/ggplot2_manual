### Prepare data for plotting ###

library(pasilla)
library(DESeq2)
library(tidyverse)
library(org.Dm.eg.db)

pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package = "pasilla",
                      mustWork = TRUE)
pasAnno <- system.file(
  "extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla",
  mustWork = TRUE
)
cts <- as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
coldata <- read.csv(pasAnno, row.names = 1)
coldata <- coldata[, c("condition", "type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

rownames(coldata) <- sub("fb", "", rownames(coldata))
coldata <- coldata[colnames(cts), ]

dds <-
  DESeqDataSetFromMatrix(countData = cts,
                         colData = coldata,
                         design =  ~ condition)
dds <- DESeq(dds)
res <- results(dds)

#### assign genes to chromosomes ####
x <- org.Dm.egCHR
# Get the entrez gene identifiers that are mapped to a chromosome
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
xx <-
  t(as.data.frame(xx)) %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
  mutate(ENTREZID = gsub("X", "", gene_id),CHR = V1) %>% select(ENTREZID,CHR)

mapping_genes <-
  AnnotationDbi::select(
    org.Dm.eg.db,
    keys = rownames(res),
    columns = "ENTREZID",
    keytype = "ENSEMBL"
  ) %>% na.omit() %>% left_join(xx)


res <- as.data.frame(res) %>% rownames_to_column("ENSEMBL")

matrix_norm <- as.data.frame(counts(dds, normalized = T))
matrix_norm <- matrix_norm %>% rownames_to_column("ENSEMBL") %>% pivot_longer(cols = colnames(matrix_norm))


save(matrix_norm, res, mapping_genes, coldata, file = "course_data.rdata")
