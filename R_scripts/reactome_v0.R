library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ReactomePA)
library(biomaRt)

load("downstream_inputs/DE_DKO_adj_inputs.RData")

# create a named log2foldchange vector for reactome input
log2fc <- res_DKO_adj$log2FoldChange
names(log2fc) <- rownames(res_DKO_adj)
log2fc <- log2fc[!is.na(log2fc)]
log2fc <- sort(log2fc, decreasing = TRUE)

# create lookup table of ENSEMBL to entrez ids cuz reactome friggen doesnt use ENSEMBL
ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=115)
mapping <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = names(log2fc),
  mart = ensembl
)

# dang. 48% of ensembl IDs failed to map to entrez
sum(is.na(mapping)) / dim(mapping)[1]

# rename log2fc vector with entrez ids
entrez_ids <- entrez_ids[!is.na(entrez_ids)]
log2fc_entrez <- log2fc[names(entrez_ids)]
names(log2fc_entrez) <- entrez_ids

# reactomePA ORA analysis
ora <- enrichPathway(
  names(log2fc_entrez),
  organism = "mouse",
)

# save plot to pdf file
pdf(file = "plots/reactome.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

dotplot(ora, title="ReactomePA DKO Case vs Control")

dev.off()
