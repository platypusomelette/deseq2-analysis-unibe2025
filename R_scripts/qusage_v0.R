library(qusage)
library(dplyr)

# Create design vector: 0 = reference/control, 1 = treatment/experiment
group_factor <- factor()
group_numeric <- as.numeric(group_factor) - 1  # Control = 0, Treatment = 1

# qusage expects a matrix: genes x samples
# Ensure rownames are gene symbols or IDs
expr_matrix <- as.matrix(expr_log2)

# ---------------------------
# Run qusage
# ---------------------------

qusage_res <- qusage(expr_matrix, 
                     group_numeric, 
                     gene_sets, 
                     n.iter = 2000)  # increase iterations if needed