library(orthogene)

convert_orthologous_genes <- function(genes,
                          input_species = "human", 
                          output_species = "mouse") {
    
    genes <- data.frame(genes)
    colnames(genes) <- "gene"
    
    converted_genes <- convert_orthologs(
        gene_df = genes,
        gene_input = "gene",
        input_species = input_species,
        output_species = output_species
    )
    return(rownames(converted_genes))
}
