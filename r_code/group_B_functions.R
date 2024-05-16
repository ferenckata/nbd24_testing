#' Creating a directory
#' 
#' @param dir_path Path to the directory (mandatory)
#' @returns None
#' @examples
#' create_directory("output/results")
create_directory <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
}


#' Creating a simulated count matrix
#' 
#' @param output_path Path to RData file, where count matrix is saved (mandatory)
#' @param number_of_genes Number of genes to have in the matrix (default is 1000)
#' @param number_of_samples Number of samples to have in the matrix (default is 10)
#' @param seed Random seed for reproducibility (default is no seed)
#' @returns count_matrix
#' @examples
#' generate_counts(
#'  output_path = "output/results/count_matrix.RData",
#'  number_of_genes = 1000,
#'  number_of_samples = 10,
#'  seed = 567)
generate_counts <- function(
    output_path = path_to_output_rdata,
    number_of_genes = 1000,
    number_of_samples = 10,
    seed = NA) {
    
    if (!is.na(seed)) {
        set.seed(seed)
    }

    count_matrix <- matrix(
        data = round(runif(number_of_genes * number_of_samples, min = 0, max = 100)),
        nrow = number_of_genes)
    colnames(count_matrix) <- paste0("Sample", 1:number_of_samples)
    rownames(count_matrix) <- paste0("Gene", 1:number_of_genes)

    save(
        count_matrix,
        file = output_path)
    return(count_matrix)
}


#' Creating a simulated count matrix and associated metadata
#' 
#' @param output_path Path to RData file, where count matrix is saved (mandatory)
#' @param number_of_genes Number of genes to have in the matrix (default is 1000)
#' @param number_of_samples Number of samples to have in the matrix (default is 10)
#' @param seed Random seed for reproducibility (default is no seed)
#' @param condition_vector Names of different conditions (default is "Treatment" and "Control")
#' @param metadata_output_path Path to the .tsv file, where metadata is saved (mandatory)
#' @returns list(counts = count_matrix, metadata = metadata)
#' @examples
#' generate_counts_with_metadata(
#'  counts_output_path = COUNT_MAT_PATH,
#'  number_of_genes = 1000,
#'  number_of_samples = 10,
#'  seed = 567,
#'  condition_vector = c("Control", "Treatment"),
#'  metadata_output_path = METADATA_PATH)
generate_counts_with_metadata <- function(
    counts_output_path = path_to_output_rdata,
    number_of_genes = 1000,
    number_of_samples = 10,
    seed = NA,
    condition_vector = c("Control", "Treatment"),
    metadata_output_path = path_to_metadata) {

    ## Generating counts:
    count_matrix <- generate_counts(
        output_path = counts_output_path,
        number_of_genes = number_of_genes,
        number_of_samples = number_of_samples,
        seed = seed)

    ## Generating metadata:
    metadata <- data.frame(
        Sample = colnames(count_matrix),
        Condition = sample(
            condition_vector, 
            size = number_of_samples,
            replace = TRUE),
        stringsAsFactors = FALSE)
    
    data.table::fwrite(
        metadata,
        file = metadata_output_path,
        sep = "\t",
        col.names = TRUE)

    return(list(
        counts = count_matrix,
        metadata = metadata))

}
