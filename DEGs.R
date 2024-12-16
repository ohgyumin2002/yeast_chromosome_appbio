
#install.packages("BiocManager")
library(BiocManager)
library(dplyr)
library(ggplot2)
#BiocManager::install("DESeq2")
library(DESeq2)
setwd("../deseq/")
# Load the counts matrix (change path to wherever you put allCounts to)
counts <- read.csv("allCounts.csv", sep=",", row.names = NULL)

# Display the loaded counts data
View(counts)

# Ensure unique row names
uniq_name <- make.names(counts$X, unique = TRUE)
row.names(counts) <- uniq_name

View(counts)

counts <- counts[, c(-1)]

# Remove the first column (assumed to be gene names) to retain only the count data
mat <- counts

# Display the processed counts matrix
View(counts)
dim(counts)

final_counts <- counts

View(final_counts)
mat <- final_counts

# Create a metadata dataframe specifying conditions for each sample
samples <- data.frame(theSample = rep(c("wt", "sy14"), each=6))
View(samples)


# Create DESeq2 dataset
ds <- DESeqDataSetFromMatrix(countData=mat, colData=samples, design=~theSample)

# Assigning sample names in the data set object may be useful in downstream analysis
colnames(ds) <- colnames(mat)

# Run DESeq2 analysis
ds <- DESeq(ds)

# Extract results
res <- results(ds)
View(res)



# Filter significant results based on adjusted p-value (e.g., p < 0.05)
sig <- res[which(res$padj < 0.001), ]
sig

# Further filter results based on log2 fold-change thresholds
sigLf <- sig[ which(sig$log2FoldChange < -1 | sig$log2FoldChange > 1), ]
sigLf

# Save results to csv file
write.csv(as.data.frame(res), file = "deseq2_results.csv")

# MA plot
plotMA(ds)

# PCA plot
data_vst <- vst(ds, blind = FALSE)
plotPCA(data_vst, intgroup = "Condition")

