# Load required libraries
library(DESeq2)
library(ggplot2)
library(data.table)

#################################
## Defining output directories ##
#################################

## Main output directory:
OUTPUT_DIR <- file.path("results")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

## Generated data:
DATA_DIR <- file.path(OUTPUT_DIR, "data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

#########################
## Simulating the data ##
#########################

## DATA 1
## Create example RNA-seq count matrix and export it

COUNT_MAT_PATH <- file.path(DATA_DIR, "count_matrix.RData")

set.seed(123)
num_genes <- 1000
num_samples <- 10

count_matrix <- matrix(data = round(runif(num_genes * num_samples, min = 0, max = 100)), nrow = num_genes)
colnames(count_matrix) <- paste0("Sample", 1:num_samples)
rownames(count_matrix) <- paste0("Gene", 1:num_genes)

save(
    count_matrix,
    file = COUNT_MAT_PATH)

## DATA 2
## Create metadata and export it
metadata <- data.frame(
    Sample = colnames(count_matrix),
    Condition = sample(c("Control", "Treatment"), size = num_samples, replace = TRUE),
    stringsAsFactors = FALSE)



####################################
## Reading and checking the input ##
####################################

# Step 2: Read the exported matrix file and check for numerical values and missing values
count_matrix <- read.csv("count_matrix.csv", header = TRUE, row.names = 1)
if (!all(sapply(count_matrix, is.numeric)) || any(is.na(count_matrix))) {
  stop("Matrix contains non-numeric values or missing values")
}



# Step 3: Impute missing values using average
count_matrix[is.na(count_matrix)] <- rowMeans(count_matrix, na.rm = TRUE)




# Step 4: Perform differential gene expression using DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ Condition)
dds <- DESeq(dds)


## Pre-filter lowly expressed genes:
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]



# Extract differential expression results
res <- results(dds)

# Step 5: Plot differential gene expression results using ggplot2
# For example, let's plot a volcano plot
volcano_plot <- ggplot(data = res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1.5, color = ifelse(res$padj < 0.05, "red", "black")) +
  labs(x = "Log2 Fold Change", y = "-Log10 p-value", title = "Volcano Plot") +
  theme_minimal()

print(volcano_plot)
