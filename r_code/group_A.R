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

## Produced results:
RESULTS_DIR <- file.path(OUTPUT_DIR, "results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

#########################
## Simulating the data ##
#########################

COUNT_MAT_PATH <- file.path(DATA_DIR, "count_matrix.RData")
METADATA_PATH  <- file.path(DATA_DIR, "metadata.tsv")

## DATA 1
## Create example RNA-seq count matrix and export it

set.seed(567)
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

data.table::fwrite(
    metadata,
    file = METADATA_PATH,
    sep = "\t",
    col.names = TRUE)

####################################
## Reading and checking the input ##
####################################

## Loading the count matrix:
count_matrix <- get(load(COUNT_MAT_PATH))

## Reading the metadata:
metadata <- data.table::fread(
    METADATA_PATH,
    header = TRUE,
    sep = "\t")

## Checking the input:
if (!all(sapply(count_matrix, is.numeric))) {
    stop("Matrix contains non-numeric values.")
}

if (any(is.na(count_matrix))) {
    warning("Matrix contains missing values. Consider adding pseudocounts or data imputation.")
}

##########################
## Processing the input ##
##########################

## Impute missing values using average:
count_matrix[is.na(count_matrix)] <- rowMeans(count_matrix, na.rm = TRUE)

## Add pseudocount:
pseudocount <- 1
count_matrix <- count_matrix + pseudocount

## Filtering lowly expressed genes:
min_counts <- 10
keep_genes <- rowSums(count_matrix >= min_counts) > 0
count_matrix <- count_matrix[keep_genes, ]

#######################################################
## Performming differential gene expression analysis ##
#######################################################

## Perform differential gene expression using DESeq2
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ Condition)
dds <- DESeq(dds)

## Extract differential expression results
res <- results(dds)

###################################################
## Plotting differential gene expression results ##
###################################################

results_for_plotting <- as.data.frame(res)

volcano_plot <- ggplot(
        data = results_for_plotting,
        aes(
            x = log2FoldChange,
            y = -log10(pvalue))) +
    geom_point(
        size = 1.5,
        color = ifelse(res$padj < 0.05, "red", "black")) +
    labs(
        x = "Log2 Fold Change",
        y = "-Log10 p-value",
        title = "DGE results") +
    theme_minimal()

###########################
## Exporting the results ##
###########################

PATH_TO_DGE_RES <- file.path(RESULTS_DIR, "DGE_results.tsv")
PATH_TO_DGE_PLOT <- file.path(RESULTS_DIR, "DGE_results.pdf")

data.table::fwrite(
    results_for_plotting,
    file = PATH_TO_DGE_RES,
    sep = "\t",
    col.names = TRUE)

ggsave(
    volcano_plot,
    file = PATH_TO_DGE_PLOT,
    height = 6,
    width  = 6,
    limitsize = FALSE)


