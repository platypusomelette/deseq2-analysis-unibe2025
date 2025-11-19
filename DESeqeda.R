library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Step 1: Read counts file and process for DESeq
counts_data <- read.table("counts.txt", 
                          header = TRUE, 
                          row.names = 1,
                          sep = "\t",
                          check.names = FALSE)
# Check dimensions
dim(counts_data)
head(counts_data)

# Clean up column names - extract just SRR numbers
colnames(counts_data) <- gsub(".*/", "", colnames(counts_data))
colnames(counts_data) <- gsub(".sorted.bam", "", colnames(counts_data))

# Check cleaned names
colnames(counts_data)

# Convert to matrix and ensure numeric
counts_data <- as.matrix(counts_data)
mode(counts_data) <- "numeric"

# Check for NAs
sum(is.na(counts_data))

# Remove any rows with NAs if they exist
counts_data <- counts_data[complete.cases(counts_data), ]

# Check again
sum(is.na(counts_data))

# Create coldata
coldata <- data.frame(
  row.names = colnames(counts_data),
  condition = c(
    rep("Disease WT", 5), # SRR7821949-7821953 
    rep("Control WT", 3), # SRR7821968-7821970
    rep("Disease Ifnar-/- x Ifngr-/-", 4), # SRR7821954-7821957
    rep("Control Ifnar-/- x Ifngr-/-", 3)  # SRR7821971-7821973
  )
)
coldata$condition <- factor(coldata$condition)

# Verify order matches
all(rownames(coldata) == colnames(counts_data))

# Step 2: Create DESeq Dataset and run 

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = coldata,
  design = ~ condition
)

# Run DESeq
dds <- DESeq(dds)

# VST with blind=TRUE
vsd <- vst(dds, blind = TRUE)

# Optional try rlog()
#rld <- rlog(dds, blind = TRUE)


# Step 3: PCA Plot 

# PCA plot
plotPCA(vsd, intgroup = "condition")


#Step 4: Heatmap 


# Select top 20 genes by mean expression
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]

# Create annotation dataframe 
df <- as.data.frame(colData(dds)[, "condition", drop = FALSE])


# Step 4a: Basic heatmap run this 

# Create heatmap using vsd (variance-stabilized data)
pheatmap(assay(vsd)[select, ], 
         cluster_rows = FALSE, 
         show_rownames = FALSE,
         cluster_cols = FALSE, 
         annotation_col = df)

# Step 4b: Run this too for clustered heatmap 

# Select top 20 genes by mean expression
select <- order(rowMeans(counts(dds, normalized = TRUE)),
                decreasing = TRUE)[1:20]

# Create annotation dataframe
df <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

sampleDists <- dist(t(assay(vsd)))


pheatmap(assay(vsd)[select, ], 
         cluster_rows = TRUE, 
         show_rownames = TRUE,
         cluster_cols = TRUE, 
         annotation_col = df,
         main = "Top 20 Expressed Genes")


