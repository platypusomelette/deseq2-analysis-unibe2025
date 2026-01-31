# load supplementary data 2: Genes in the lung modules
modules_genelist_raw <- read.csv("raw_inputfiles/41467_2019_10601_MOESM4_ESM_minimal.csv", header=TRUE)
# load supplementary data 4: Annotation of the lung modules
modules_annot_raw <- read.csv("raw_inputfiles/41467_2019_10601_MOESM6_ESM_minimal.csv", header=TRUE)

# create lookup table for custom gene sets using supplementary data ---------------------------------
modules_genelist <- modules_genelist_raw[, c(3,1)]
modules_annot_names <- modules_annot_raw[, c(1,2)]

# append functional annotation column to include module ID in parentheses for human-readable labels
modules_annot_names[,2] <- paste0(modules_annot_names[,2], " (", modules_annot_names[,1], ")")
row.names(modules_annot_names) <- modules_annot_raw[,1]

### define function customGO() to run custom overrepresentation analysis
customGO <- function(DE, DE_upregulated, DE_downregulated, title, subtitle = "") {
  
  # custom overrepresentation - all genes ----------------------------------------
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
  
  # custom overrepresentation - upregulated genes ----------------------------------------
  ego_custom_up <- enricher(
    gene = rownames(DE_upregulated),
    universe = counts_file$Geneid,
    TERM2GENE = modules_genelist,
    pvalueCutoff = 0.05
  )
  
  # swapping out plain module labels for human-readable labels
  ego_custom_up@result$Description <- modules_annot_names[ego_custom_up@result$ID, 2]
  
  # create dotplot
  up_dotplot <- dotplot(ego_custom_up, showCategory = 15) + 
    ggtitle(paste0("Custom Modules - Upregulated in ", title)) +
    labs(subtitle = subtitle)
  
  # custom overrepresentation - downregulated genes ----------------------------------------
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
pdf(file = "plots/customGO_Lmodules.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8) # The height of the plot in inches

print(all_DE_diseased_dotplot_custom)
print(up_DE_diseased_dotplot_custom)
print(down_DE_diseased_dotplot_custom)

dev.off()