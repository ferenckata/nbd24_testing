#############################
## Load required libraries ##
#############################

library(DESeq2)
library(ggplot2)
library(data.table)

############################
## Sourcing the functions ##
############################

source("r_code/group_B_functions.R")

#################################
## Defining output directories ##
#################################

OUTPUT_DIR <- "results"
create_directory(OUTPUT_DIR)

DATA_DIR <- file.path(OUTPUT_DIR, "data")
create_directory(DATA_DIR)

RESULTS_DIR <- file.path(OUTPUT_DIR, "results")
create_directory(RESULTS_DIR)

#########################
## Simulating the data ##
#########################

COUNT_MAT_PATH <- file.path(DATA_DIR, "count_matrix.RData")
METADATA_PATH  <- file.path(DATA_DIR, "metadata.tsv")

counts_metadata_list <- generate_counts_with_metadata(
    counts_output_path = COUNT_MAT_PATH,
    number_of_genes = 1000,
    number_of_samples = 10,
    seed = 567,
    condition_vector = c("Control", "Treatment"),
    metadata_output_path = METADATA_PATH
)

## DATA 1
## Create example RNA-seq count matrix and export it

count_matrix <- counts_metadata_list$counts

## DATA 2
## Create metadata and export it

metadata <- counts_metadata_list$metadata

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