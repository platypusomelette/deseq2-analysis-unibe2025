# load supplementary data 3: Genes in the blood modules with average normalised read counts for each group
modules_genelist_raw <- read.csv("raw_inputfiles/41467_2019_10601_MOESM5_ESM_minimal.csv", header=TRUE)
# load supplementary data 5: Annotation of the blood modules
modules_annot_raw <- read.csv("raw_inputfiles/41467_2019_10601_MOESM7_ESM_minimal.csv", header=TRUE)

# create lookup table for custom gene sets using supplementary data ---------------------------------
modules_genelist <- modules_genelist_raw[, c(3,1)] # maps ENSEMBL IDs to module ID
modules_annot_names <- modules_annot_raw[, c(1,4)] # maps module ID to functional annotation

# append functional annotation column to include module ID in parentheses for human-readable labels
modules_annot_names[,2] <- paste0(modules_annot_names[,2], " (", modules_annot_names[,1], ")")
row.names(modules_annot_names) <- modules_annot_raw[,1]

### define function customGO() to run custom gene set enrichment analysis
customGO <- function(DE, DE_upregulated, DE_downregulated, title, subtitle = "") {
  
  # custom gene set enrichment - all genes ----------------------------------------
  ego_custom <- enricher(
    gene = rownames(DE),
    universe = counts_file$Geneid,
    TERM2GENE = modules_genelist,
    pvalueCutoff = 0.05
  )
  
  # swapping out plain module labels for human-readable labels
  ego_custom@result$Description <- modules_annot_names[ego_custom@result$ID, 2]
  
  # create dotplot
  all_dotplot <- dotplot(ego_custom) + 
    ggtitle(paste0("Custom Modules - ", title)) +
    labs(subtitle = subtitle)
  
  # custom gene set enrichment - upregulated genes ----------------------------------------
  ego_custom_up <- enricher(
    gene = rownames(DE_upregulated),
    universe = counts_file$Geneid,
    TERM2GENE = modules_genelist,
    pvalueCutoff = 0.05
  )
  
  # swapping out plain module labels for human-readable labels
  ego_custom_up@result$Description <- modules_annot_names[ego_custom_up@result$ID, 2]
  
  # create dotplot
  up_dotplot <- dotplot(ego_custom_up) + 
    ggtitle(paste0("Custom Modules - Upregulated in ", title)) +
    labs(subtitle = subtitle)
  
  # custom gene set enrichment - downregulated genes ----------------------------------------
  ego_custom_down <- enricher(
    gene = rownames(DE_downregulated),
    universe = counts_file$Geneid,
    TERM2GENE = modules_genelist,
    pvalueCutoff = 0.05
  )
  
  # swapping out plain module labels for human-readable labels
  ego_custom_down@result$Description <- modules_annot_names[ego_custom_down@result$ID, 2]
  
  # create dotplot
  down_dotplot <- dotplot(ego_custom_down) + 
    ggtitle(paste0("Custom Modules - Downregulated in ", title)) +
    labs(subtitle = subtitle)
  
  return(list(
    ego_custom = ego_custom,
    all_dotplot = all_dotplot,
    ego_custom_up = ego_custom_up,
    up_dotplot = up_dotplot,
    ego_custom_down = ego_custom_down,
    down_dotplot = down_dotplot
  ))

}

# do custom gene set analysis for DE_diseased
load("downstream_inputs/DE_diseased_inputs.RData")
DE_diseased_customGO <- customGO(DE_diseased,DE_diseased_upregulated, DE_diseased_downregulated, "DKO vs WT in Disease Case")
all_DE_diseased_dotplot_custom <- DE_diseased_customGO$all_dotplot
up_DE_diseased_dotplot_custom <- DE_diseased_customGO$up_dotplot
down_DE_diseased_dotplot_custom <- DE_diseased_customGO$down_dotplot

# save outputs ------------------------------------------------------------

# save plots to pdf files
pdf(file = "plots/customGO.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

print(all_DE_diseased_dotplot_custom)
print(up_DE_diseased_dotplot_custom)
print(down_DE_diseased_dotplot_custom)

dev.off()

# rank modules by percentage of genes annotated as belonging to type I or type II in supp data 5  --------

modules_genelist_raw$Type.I <- modules_genelist_raw$Type.I != ""
modules_genelist_raw$Type.II <- modules_genelist_raw$Type.II != ""

typeI_hits <- sapply(modules_annot_names[,1],
       function(x) sum(modules_genelist_raw[modules_genelist_raw[,3] == x, 4], na.rm = TRUE)
)
typeII_hits <- sapply(modules_annot_names[,1],
                     function(x) sum(modules_genelist_raw[modules_genelist_raw[,3] == x, 5], na.rm = TRUE)
)
module_totals <- sapply(modules_annot_names[,1],
                      function(x) sum(modules_genelist_raw[,3] == x, na.rm = TRUE)
)

typeI_ranks <- data.frame(sort(typeI_hits / module_totals, decreasing=TRUE))
typeII_ranks <- data.frame(sort(typeII_hits / module_totals, decreasing=TRUE))

# save ranks to tab delimited txt
write.table(typeI_ranks, file = "text_outputs/customGO_typeI_ranks.txt", 
            sep = "\t", row.names = TRUE, col.names = "type_I_percent", quote = FALSE)
write.table(typeII_ranks, file = "text_outputs/customGO_typeII_ranks.txt", 
            sep = "\t", row.names = TRUE, col.names = "type_II_percent", quote = FALSE)


# make dfs for labeling type I or type II genes ----------------------------

typeI_genelist <- modules_genelist_raw[modules_genelist_raw$Type.I,]
typeII_genelist <- modules_genelist_raw[modules_genelist_raw$Type.II,]

save(typeI_genelist, typeII_genelist,
     file = "downstream_inputs/IFNtype_genelist_inputs.RData")


# save DE genes in DKO_adj that are not IFN genes, to tsv --------

# doing search with ENSEMBL IDs then getting the gene symbols back with renamed df, just to maintain consistency with paper's genelist
load("downstream_inputs/DE_DKO_adj_inputs.RData")
DKO_adj_tophits <- DE_DKO_adj[order(DE_DKO_adj$padj),]
DE_DKO_adj_renamed <- data.frame(DE_DKO_adj_renamed[order(DE_DKO_adj_renamed$padj),])
DKO_adj_notIFN_tophits <- DE_DKO_adj_renamed[which(!rownames(DKO_adj_tophits) %in% c(typeI_genelist$Gene, typeII_genelist$Gene)),]

write.table(DKO_adj_notIFN_tophits[6:7], file = "text_outputs/DKO_adj_notIFN_tophits.tsv")
