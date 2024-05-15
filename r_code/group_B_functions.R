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
}
