## To create directory:

create_directory <- function(dir_path) {
    if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
}

## To create count matrix:

generate_counts <- function(
    output_path = path_to_output_rdata,
    number_of_genes = 1000,
    number_of_samples = 10,
    seed = 567) {
    
    set.seed(seed)

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


## Generate counts with metadata:

generate_counts_with_metadata <- function(
    counts_output_path = path_to_output_rdata,
    number_of_genes = 1000,
    number_of_samples = 10,
    seed = 567,
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
