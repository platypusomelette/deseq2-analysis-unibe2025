library(ggplot2)
library(enrichplot)
library(ggrepel)

load("downstream_inputs/IFNtype_genelist_inputs.RData")

# define volcanoplotter functions -----------------------------------------

volcanoplotter <- function(res, DE, DE_renamed, title, ylim = c(0, 80)) {

  # creating new df out of DE_renamed (gene IDs instead of ENSEMBL) for sig genes labeling
  DE_df <- as.data.frame(DE_renamed)
  volcano_labels <- DE_df[DE_df$padj < 0.05 & abs(DE_df$log2FoldChange) > 1,]
  volcano_labels <- subset(volcano_labels, log2FoldChange >= -10 & log2FoldChange <= 10) # so we stop labeling out-of-bounds dots
  
  # make volcano plot 
  DE_volcano <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = direction)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 2) + 
    geom_text_repel(data = volcano_labels, 
                    aes(label=rownames(volcano_labels)), size = 3, 
                    show.legend = FALSE) +
    scale_color_manual(values = c("cyan3", "grey", "red3")) +
    coord_cartesian(ylim = ylim, xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Direction', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
         title = title, subtitle = 'padj < 0.05 and absolute log2FoldChange > 1') +
    scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis
  
  return(list(
    DE_volcano = DE_volcano,
    tophits = volcano_labels[order(volcano_labels$padj),])
  )

}

IFN_volcanoplotter <- function(res, DE, DE_renamed, title, ylim = c(0, 80)) {
  
  # creating new df out of DE_renamed (gene IDs instead of ENSEMBL) for sig genes labeling
  DE_df <- as.data.frame(DE_renamed)
  volcano_labels <- DE_df[DE_df$padj < 0.05 & abs(DE_df$log2FoldChange) > 1,]
  
  # make new volcano labels for IFN genes
  IFN_genelist <- c(typeI_genelist$Gene, typeII_genelist$Gene)
  volcano_labels_IFN <- DE_df[rownames(DE) %in% IFN_genelist, ]
  volcano_labels_IFN <- subset(volcano_labels_IFN, log2FoldChange >= -10 & log2FoldChange <= 10) # so we stop labeling out-of-bounds dots
  
  # make volcano plot with IFN labels
  DE_volcano <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = direction)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 2) + 
    scale_color_manual(values = c("cyan3", "grey", "red3")) +
    coord_cartesian(ylim = ylim, xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
    labs(color = 'Direction', #legend_title, 
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value"),
         title = title, subtitle = 'IFN genes labeled. padj < 0.05 and absolute log2FoldChange > 1') +
    scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
    geom_point(
      data = volcano_labels_IFN,
      color = "orange2",
      size = 2) +
    geom_text_repel(data = volcano_labels_IFN, 
                    aes(label=rownames(volcano_labels_IFN)), size = 3, 
                    show.legend = FALSE,
                    color = "orange2")
  return(list(
    DE_volcano = DE_volcano,
    tophits = volcano_labels[order(volcano_labels$padj),])
  )
  
}

# make and save da volcano plots ---------------------------------------------------
# save plot to pdf file
pdf(file = "plots/volcanoplot.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6) # The height of the plot in inches

# DE_DKO_adj
load("downstream_inputs/DE_DKO_adj_inputs.RData")
DE_DKO_adj_volcano <- volcanoplotter(res_DKO_adj, DE_DKO_adj, DE_DKO_adj_renamed, 'DKO Diseased vs Healthy Control (adjusted with respect to WT case vs control)')
DE_DKO_adj_volcano_IFN <- IFN_volcanoplotter(res_DKO_adj, DE_DKO_adj, DE_DKO_adj_renamed, 'DKO Diseased vs Healthy Control (adjusted with respect to WT case vs control)')
print(DE_DKO_adj_volcano$DE_volcano)
print(DE_DKO_adj_volcano_IFN$DE_volcano)
# save top hits to tab delimited txt
write.table(DE_DKO_adj_volcano$tophits, file = "DE_DKO_adj_volcano_tophits.txt", 
            sep = "\t", row.names = TRUE, quote = FALSE)

# DE_WT
load("downstream_inputs/DE_WT_inputs.RData")
DE_WT_volcano <- volcanoplotter(res_WT, DE_WT, DE_WT_renamed, 'WT Diseased vs Healthy Control', ylim = c(0,280))
DE_WT_volcano_IFN <- IFN_volcanoplotter(res_WT, DE_WT, DE_WT_renamed, 'WT Diseased vs Healthy Control', ylim = c(0,280))
print(DE_WT_volcano$DE_volcano)
print(DE_WT_volcano_IFN$DE_volcano)
# save top hits to tab delimited txt
write.table(DE_WT_volcano$tophits, file = "DE_WT_volcano_tophits.txt", 
            sep = "\t", row.names = TRUE, quote = FALSE)

# DE_DKO
load("downstream_inputs/DE_DKO_inputs.RData")
DE_DKO_volcano <- volcanoplotter(res_DKO, DE_DKO, DE_DKO_renamed, 'DKO Diseased vs Healthy Control', ylim = c(0,280))
DE_DKO_volcano_IFN <- IFN_volcanoplotter(res_DKO, DE_DKO, DE_DKO_renamed, 'DKO Diseased vs Healthy Control', ylim = c(0,280))
print(DE_DKO_volcano$DE_volcano)
print(DE_DKO_volcano_IFN$DE_volcano)
# save top hits to tab delimited txt
write.table(DE_DKO_volcano$tophits, file = "DE_DKO_volcano_tophits.txt", 
            sep = "\t", row.names = TRUE, quote = FALSE)

# DE_diseased
load("downstream_inputs/DE_diseased_inputs.RData")
DE_diseased_volcano <- volcanoplotter(res_diseased, DE_diseased, DE_diseased_renamed, 'DKO vs WT in disease cases', ylim = c(0,280))
DE_diseased_volcano_IFN <- IFN_volcanoplotter(res_diseased, DE_diseased, DE_diseased_renamed, 'DKO vs WT in disease cases', ylim = c(0,280))
print(DE_diseased_volcano$DE_volcano)
print(DE_diseased_volcano_IFN$DE_volcano)
# save top hits to tab delimited txt
write.table(DE_diseased_volcano$tophits, file = "DE_diseased_volcano_tophits.txt", 
            sep = "\t", row.names = TRUE, quote = FALSE)

# save DE genes in DKO_adj that are not in WT, to tab delimited txt
DKO_adj_vs_WT_tophits <- rownames(DE_DKO_adj_volcano$tophits)[!rownames(DE_DKO_adj_volcano$tophits) %in% rownames(DE_WT_volcano$tophits)]

DKO_adj_vs_WT_tophits <- DE_DKO_adj_volcano$tophits[DKO_adj_vs_WT_tophits,]

write.table(DKO_adj_vs_WT_tophits[6], file = "DKO_adj_vs_WT_volcano_tophits.tsv")

# save plots
dev.off()