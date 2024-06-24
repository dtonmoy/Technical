# ##################################################
# LOAD LIBRARY
library(Seurat)
library(ggplot2)
library(tibble)
library(dplyr)
library(openxlsx)
library(SeuratDisk)
library(rhdf5)
library(grDevices)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON = "/Users/tonmoy/Research/Spatial_proximity_project/env/decoupleR/bin/python")
use_condaenv(condaenv = "decoupleR", conda = "auto", required = NULL)

dc <- import("decoupler")
pd <- import("pandas")
np <- import("numpy")
scipy <- import("scipy")
sklearn <- import("sklearn")


# ##################################################
# FUNCTIONS
# ----------
# General
# ----------
convert_rownames_to_column <- function(matrix, colname) {
    
    df <- as.data.frame(matrix)
    df[[colname]] <- rownames(df)
    rownames(df) <- NULL
    shuffled_df <- df[, c(ncol(df), 1:(ncol(df)-1))]
    
    return(shuffled_df)
    
}
convert_column_to_rownames <- function(matrix, colname) {
    
    countMat <- as.data.frame(matrix)
    countMat <- na.omit(countMat)
    rownames(countMat) <- countMat[[colname]] 
    countMat <- countMat[, -which(colnames(countMat) == colname)]
    
    return(countMat)
    
}
generate_random_palette <- function(n, seed=1234) {
    
    set.seed(seed)
    
    # ----------
    # Get all available colors
    # ----------
    all_colors <- colors()
    
    # ----------
    # Randomly sample n colors
    # ----------
    palette <- sample(all_colors, n)
    
    return(palette)
}

# ----------
# Seurat specific
# ----------
seurat_fix_dimnames <- function(seurat_object, other_layers=TRUE) {
    
    if (class(seurat_object@assays[[1]]) == "Assay5") {
        
        other_layers_names <- names(seurat_object@assays)
        
        for (i in other_layers_names) {
            
            if (typeof(seurat_object@assays[[i]]@layers[["data"]]) == "NULL") {
                
                # -------------------------
                # fix the stupid dimension names of counts
                # -------------------------
                seurat_object@assays[[i]]@layers[["counts"]]@Dimnames[[1]] <- rownames(seurat_object@assays[[i]])
                seurat_object@assays[[i]]@layers[["counts"]]@Dimnames[[2]] <- colnames(seurat_object@assays[[i]])
                
            } else if (typeof(seurat_object@assays[[i]]@layers[["data"]]) != "NULL") {
                # browser()
                # -------------------------
                # fix the stupid dimension names of both counts and data
                # -------------------------
                seurat_object@assays[[i]]@layers[["counts"]]@Dimnames[[1]] <- rownames(seurat_object@assays[[i]])
                seurat_object@assays[[i]]@layers[["counts"]]@Dimnames[[2]] <- colnames(seurat_object@assays[[i]])
                seurat_object@assays[[i]]@layers[["data"]]@Dimnames[[1]] <- rownames(seurat_object@assays[[i]])
                seurat_object@assays[[i]]@layers[["data"]]@Dimnames[[2]] <- colnames(seurat_object@assays[[i]])
                
            }
        }
        
    } else if (class(seurat_object@assays[[1]]) == "Assay") {
        cat("The version of this seurat object is 4.\nNo need to fix the dimnames!!!")
    }
    
    return(seurat_object)
    
}
get_hkg_list <- function() {
    
    # ----------
    # Load and clean the housekeeping genes
    # ----------
    hkg <- read.csv("/Users/tonmoy/Research/GBM_3/data/genesets/human_housekeeping_genesets.csv")
    hkg <- subset(hkg, select = -references)
    hkg <- hkg[-c(2, 3, 4), ]
    hkg[1, 1] <- NA
    hkg <- lapply(colnames(hkg), function(x) {na.omit(hkg[[x]])})
    names(hkg) <- c("gtex", "hpa", "hrt", "cellminer", "klijn")
    return(hkg)
}
seurat_QC <- function(seurat_object) {
    
    # -------------------------
    # calculate MT content
    # -------------------------
    seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent_mt")
    
    # -------------------------
    # calculate ribosomal content
    # -------------------------
    seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^(RP[1|2]|RPS\\d|RPSA|RPL\\d|RPLP\\d|UBA52|DAP3|FAU)", col.name = "percent_ribo")
    
    # -------------------------
    # calculate hemoglobin content
    # -------------------------
    seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^HBA1|^HBA2|^HBB|^HBD|^HBE1|^HBG1|^HBG2|^HBM|^HBQ1|^HBZ", col.name = "percent_hb")
    
    # -------------------------
    # calculate pseudo gene content
    # -------------------------
    seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^(AC\\d+\\.\\d+|^AL\\d+\\.\\d+|^AP\\d+\\.\\d+)", col.name = "percent_pseudo")
    
    # -------------------------
    # calculate housekeeping gene content
    # -------------------------
    hkg <- get_hkg_list()
    hkg <- lapply(hkg, function(x) {intersect(x, rownames(seurat_object))})
    seurat_object <- PercentageFeatureSet(seurat_object, features = hkg[["gtex"]], col.name = "gtex_hkg")
    seurat_object <- PercentageFeatureSet(seurat_object, features = hkg[["hpa"]], col.name = "hpa_hkg")
    seurat_object <- PercentageFeatureSet(seurat_object, features = hkg[["hrt"]], col.name = "hrt_hkg")
    seurat_object <- PercentageFeatureSet(seurat_object, features = hkg[["cellminer"]], col.name = "cellminer_hkg")
    seurat_object <- PercentageFeatureSet(seurat_object, features = hkg[["klijn"]], col.name = "klijn_hkg")
    
    return(seurat_object)
    
}
remove_noninformative_genes <- function(seurat_object) {
    
    # -------------------------
    # Get the gene list to remove (MALAT1/Mitochondrial/Ribosomal/Hb/Pseudo)
    # -------------------------
    remove_me <- rownames(seurat_object)[grepl("MALAT1|^MT-|^(RP[1|2]|RPS\\d|RPSA|RPL\\d|RPLP\\d|UBA52|DAP3|FAU)|^HBA1|^HBA2|^HBB|^HBD|^HBE1|^HBG1|^HBG2|^HBM|^HBQ1|^HBZ|^(AC\\d+\\.\\d+|^AL\\d+\\.\\d+|^AP\\d+\\.\\d+)", rownames(seurat_object))]
    
    # -------------------------
    # Remove non-informative genes
    # -------------------------
    seurat_object <- subset(seurat_object, features = rownames(seurat_object)[!rownames(seurat_object) %in% remove_me])
    
    return(seurat_object)
    
}
calculate_n_unique <- function (seurat_object) {
    
    # -----
    # Get the count matrix and calculate the genes with low spot reads and remove them
    # -----
    if (class(seurat_object@assays[[1]]) == "Assay5") {
        sparse_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
    } else if (class(seurat_object@assays[[1]]) == "Assay") {
        sparse_mat <- seurat_object@assays[["RNA"]]@counts
    }
    
    # -------------------------
    # Function to calculate unique values in a col
    # -------------------------
    calculate_unique_values <- function(col) {
        sort(unique(col[col != 0]))
    }
    
    # -----
    # Calculate unique values per row
    # -----
    unique_values_per_col <- apply(sparse_mat, 2, calculate_unique_values)
    names(unique_values_per_col) <- colnames(seurat_object)
    spot_read_length <- as.data.frame(sapply(names(unique_values_per_col), function(i) {length(unique_values_per_col[[i]])}))
    names(spot_read_length) <- "nUnique_RNA"
    
    # -----
    # Add the n_unique information in the seurat object
    # -----
    seurat_object <- AddMetaData(seurat_object, spot_read_length)
    
    return(seurat_object)
    
}
filter_low_quality_genes <- function (seurat_object, num_reads = 10, top_unique_reads = 10) {
    
    # -----
    # Get the count matrix and calculate the genes with low spot reads and remove them
    # -----
    if (class(seurat_object@assays[[1]]) == "Assay5") {
        sparse_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
    } else if (class(seurat_object@assays[[1]]) == "Assay") {
        sparse_mat <- seurat_object@assays[["RNA"]]@counts
    }
    
    print(paste0("Removing ", sum(Matrix::rowSums(sparse_mat != 0) < num_reads), " low quality genes"))
    
    seurat_object <- subset(seurat_object, features = names(Matrix::rowSums(sparse_mat != 0) > num_reads)[Matrix::rowSums(sparse_mat != 0) > num_reads])
    
    # -------------------------
    # Function to calculate unique values in a row
    # -------------------------
    calculate_unique_values <- function(row) {
        sort(unique(row[row != 0]))
    }
    
    # -----
    # Calculate unique values per row
    # -----
    unique_values_per_row <- apply(sparse_mat, 1, calculate_unique_values)
    names(unique_values_per_row) <- rownames(seurat_object)
    
    # -----
    # Convert the unique values into string
    # -----
    my_string <- lapply(names(unique_values_per_row), function(i) {paste(unique_values_per_row[[i]], collapse = ",")})
    
    # -----
    # Print the top unique values
    # -----
    all_str <- do.call(c, my_string)
    print(paste0("Printing the number of top ", top_unique_reads, " unique reads"))
    print(head(sort(table(all_str, deparse.level = 0), decreasing = TRUE), n = top_unique_reads))
    
    # -----
    # Set the names of my_string
    # -----
    names(my_string) = names(unique_values_per_row)
    
    # -----
    # Remove the genes that have only "1" or "2" or "1,2", "1,2,3" reads 
    # -----
    print("Removing genes that have only '1' or '2' or '1,2' or '1,2,3' reads")
    
    # -----
    # get the combinations
    # -----
    elements <- c("1", "2", "3")
    
    combinations_of_elements_list <- do.call("c", lapply(seq_along(elements), function(i) {utils::combn(elements, i, FUN=list)}))
    combinations_of_elements <- sapply(combinations_of_elements_list, function(i) {paste0(i, collapse = ",")})
    
    seurat_object <- subset(seurat_object, features = names(my_string[!my_string %in% combinations_of_elements]))
    print(paste0("The total number of genes after filtering : ", length(rownames(seurat_object))))
    
    return(seurat_object)
    
}
filter_low_quality_cells <- function (seurat_object, num_reads = 200, top_unique_reads = 10) {
    
    # -----
    # Get the count matrix and calculate the genes with low spot reads and remove them
    # -----
    if (class(seurat_object@assays[[1]]) == "Assay5") {
        sparse_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
    } else if (class(seurat_object@assays[[1]]) == "Assay") {
        sparse_mat <- seurat_object@assays[["RNA"]]@counts
    }
    
    # sparse_mat <- seurat_object@assays[["RNA"]]@layers[["counts"]]
    
    print(paste0("Removing ", sum(Matrix::colSums(sparse_mat != 0) < num_reads), " low quality cells"))
    
    # -----
    # subset the cells
    # -----
    seurat_object <- subset(seurat_object, cells = colnames(sparse_mat)[Matrix::colSums(sparse_mat != 0) > num_reads])
    
    # -------------------------
    # Function to calculate unique values in a row
    # -------------------------
    calculate_unique_values <- function(row) {
        sort(unique(row[row != 0]))
    }
    
    # -----
    # Calculate unique values per row
    # -----
    unique_values_per_col <- apply(sparse_mat, 2, calculate_unique_values)
    names(unique_values_per_col) <- colnames(seurat_object)
    
    # -----
    # Convert the unique values into string
    # -----
    my_string <- lapply(names(unique_values_per_col), function(i) {paste(unique_values_per_col[[i]], collapse = ",")})
    
    # -----
    # Print the top unique values
    # -----
    all_str <- do.call(c, my_string)
    print(paste0("Printing the number of top ", top_unique_reads, " unique reads"))
    print(head(sort(table(all_str, deparse.level = 0), decreasing = TRUE), n = top_unique_reads))
    
    # -----
    # Set the names of my_string
    # -----
    names(my_string) = names(unique_values_per_col)
    
    # -----
    # Remove the cells that have only "1" or "2" or "1,2" or "1,2,3" or "1,2,3,4" reads 
    # -----
    print("Removing cells that have only '1' or '2' or '1,2', '1,2,3', '1,2,3,4' reads")
    
    # -----
    # get the combinations
    # -----
    elements <- c("1", "2", "3", "4")
    
    combinations_of_elements_list <- do.call("c", lapply(seq_along(elements), function(i) {utils::combn(elements, i, FUN=list)}))
    combinations_of_elements <- sapply(combinations_of_elements_list, function(i) {paste0(i, collapse = ",")})
    
    seurat_object <- subset(seurat_object, cells = names(my_string[!my_string %in% combinations_of_elements]))
    print(paste0("The total number of cells after filtering : ", length(colnames(seurat_object))))
    
    return(seurat_object)
    
}
perform_QC <- function(seurat_object, remove_zero_genes=FALSE,
                       nFeature_RNA_threshold_low=300,
                       nFeature_RNA_threshold_high = 3000,
                       nCount_RNA_threshold_low=300, 
                       nCount_RNA_threshold_high=15000, 
                       percent_hb_threshold=3,
                       percent_mt_threshold=15, percent_pseudo_threshold=3,
                       nUnique_RNA_threshold=3, gtex_hkg_threshold=25,
                       hpa_hkg_threshold=25, hrt_hkg_threshold=25, 
                       kljin_hkg_threshold=25,
                       filter_low_quality_genes=FALSE,
                       filter_low_quality_cells=FALSE) {
    
    # ##################################################
    # QC and select cells for further analysis
    # -----
    # Remove genes that have zero counts (THIS PART IS REDUNDENT as the next functions have already these steps implemented)
    # Also if the data slot is not present in the seurat object, then this does not work (super weird)
    # -----
    if (remove_zero_genes) {
        
        # Find cells with non zero expression ---------------------------
        keep_cells <- colnames(seurat_object[, colSums(seurat_object) != 0])
        
        # Filter big matrix ---------------------------------------------
        seurat_object <- seurat_object[, keep_cells]
    }
    
    # -----
    # remove low quality genes
    # -----
    if (filter_low_quality_genes) {
        seurat_object <- filter_low_quality_genes(seurat_object)
    }
    
    # -----
    # remove low quality spots
    # -----
    if (filter_low_quality_cells) {
        seurat_object <- filter_low_quality_cells(seurat_object)
    }
    
    # -----
    # subset seurat objects based on out pre-selected threshold
    # -----
    seurat_object <- subset(seurat_object, 
                            subset = nFeature_RNA > nFeature_RNA_threshold_low &
                                nFeature_RNA < nFeature_RNA_threshold_high &
                                nCount_RNA > nCount_RNA_threshold_low &
                                nCount_RNA < nCount_RNA_threshold_high &
                                percent_hb < percent_hb_threshold &
                                percent_mt < percent_mt_threshold &
                                percent_pseudo < percent_pseudo_threshold &
                                gtex_hkg > gtex_hkg_threshold &
                                hpa_hkg > hpa_hkg_threshold &
                                hrt_hkg > hrt_hkg_threshold &
                                klijn_hkg > kljin_hkg_threshold &
                                nUnique_RNA > nUnique_RNA_threshold)
    
    return(seurat_object)
    
}
get_largest_gene <- function(seurat_object) {
    
    seurat_object$largest_count <- apply(seurat_object@assays$RNA@layers$counts, 2, max)
    seurat_object$largest_index <- apply(seurat_object@assays$RNA@layers$counts, 2, which.max)
    
    seurat_object$largest_gene <- rownames(seurat_object)[seurat_object$largest_index]
    seurat_object$percent_largest_gene <- 100 * seurat_object$largest_count / seurat_object$nCount_RNA    
    
    return(seurat_object)
    
}

# ----------
# decoupleR specific
# ----------
convert_geneset_to_list <- function(geneset_of_interest, seurat_object) {
    
    # ----------
    # Get the unique gene lists
    # ----------
    genesetList <- lapply(1:ncol(geneset_of_interest), function(i) {geneset_of_interest[, i]})
    names(genesetList) <- colnames(geneset_of_interest)
    genesetList <- lapply(genesetList, function(i) {unique(i[!is.na(i)])})
    
    # ----------
    # Remove the genes from the list if they are not present in the seurat object
    # ----------
    genesetList <- sapply(names(genesetList), function(i) {rownames(seurat_object)[rownames(seurat_object) %in% genesetList[[i]]]})
    
    return(genesetList)
    
}
get_net_df <- function(genesetList) {
    
    # ##################################################
    # Create the net dataframe
    # -----
    # append the genes in a new dataframe
    # -----
    demo_df <- data.frame(geneset = rep(names(genesetList), sapply(genesetList, length)),
                          genesymbol = unlist(genesetList))
    
    # -----
    # add collection and weight columns
    # -----
    # demo_df$collection <- "Query_geneset" # ignore this
    demo_df$weight <- 1
    
    # ##################################################
    # Remove the duplicated genes from the net dataframe
    blademo_df <- data.frame()
    
    for (i in unique(demo_df$geneset)) {
        
        # ----------
        # Subset the dataframe based on genesets
        # ----------
        genesetSub <- demo_df[demo_df$geneset == i, ]
        genesetSub <- genesetSub[!duplicated(genesetSub$genesymbol), ]
        # ----------
        # Append or assign the subsetted DataFrame to demo_df
        # ----------
        if (nrow(demo_df) == 0) {
            blademo_df <- genesetSub
        } else {
            blademo_df <- rbind(blademo_df, genesetSub)
        }
    }
    
    # ##################################################
    # Add a category column defining the type of genesets
    blademo_df$collection <- sapply(blademo_df$geneset, function(i) {
        if (grepl("random", i)) {
            strsplit(i, "_\\d+$")[[1]][1]
        } else if (grepl("rnd_per_gsType", i)) {
            "rnd_per_gsType"
        } else {
            "reference"
        }})
    
    # -----
    # set the column name as NULL
    # -----
    rownames(blademo_df) <- NULL
    
    return(blademo_df)
    
}
remove_duplicates_goi <- function(goi) {
    
    demo_df <- data.frame()
    
    for (i in unique(goi$geneset)) {
        
        # ----------
        # Subset the dataframe based on genesets
        # ----------
        genesetSub <- goi[goi$geneset == i, ]
        genesetSub <- genesetSub[!duplicated(genesetSub$genesymbol), ]
        # ----------
        # Append or assign the subsetted DataFrame to demo_df
        # ----------
        if (nrow(demo_df) == 0) {
            demo_df <- genesetSub
        } else {
            demo_df <- rbind(demo_df, genesetSub)
        }
    }
    
    return(goi)
}
run_dc_method_py <- function(seurat_object, net, source="geneset", run_type="decouple", target="genesymbol", mor="weight", 
                             methods=list("viper", "ulm", "ora"), min_n=0, data_to_misc=TRUE, include_consensus=TRUE) {
    
    # ----------
    # Get the python libraries [It's important to have these libraries 
    # installed in the background conda environment that is running in R]
    # ----------
    dc <- import("decoupler")
    pd <- import("pandas")
    
    # ----------
    # Get the matrix from the seurat object
    # ----------
    mat <- log1p(t(as.data.frame(seurat_object@assays[["RNA"]]@layers[["counts"]])))
    
    # ----------
    # Fix the dataframe formatting to run using python decoupleR
    # ----------
    mat <- pd$DataFrame(mat)
    rownames(mat) <- colnames(seurat_object)
    colnames(mat) <- rownames(seurat_object)    
    
    # ##################################################
    # Run the decoupleR method - INDEPENDENT
    if (run_type == "ulm") {
        acts <- dc$run_ulm(mat=mat, net=net, source=source, target=target,
                           min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "ora") {
        acts <- dc$run_ora(mat=mat, net=net, source=source, target=target,
                           min_n=0, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "mlm") {
        acts <- dc$run_mlm(mat=mat, net=net, source=source, target=target,
                           min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "aucell") {
        acts <- dc$run_aucell(mat=mat, net=net, source=source, target=target, 
                              min_n=0, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "gsva") {
        acts <- dc$run_gsva(mat=mat, net=net, source=source, target=target,
                            min_n=min_n, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "mdt") {
        acts <- dc$run_mdt(mat=mat, net=net, source=source, target=target,
                           min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "udt") {
        acts <- dc$run_udt(mat=mat, net=net, source=source, target=target,
                           min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "viper") {
        acts <- dc$run_viper(mat=mat, net=net, source=source, target=target,
                             min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "wmean") {
        acts <- dc$run_wmean(mat=mat, net=net, source=source, target=target,
                             min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    } else if(run_type == "wsum"){
        acts <- dc$run_wsum(mat=mat, net=net, source=source, target=target,
                            min_n=min_n, weight=mor, verbose=TRUE, use_raw=FALSE)
    }
    
    # ##################################################
    # Run the decoupleR method - CONSENSUS
    else if(run_type == "decouple"){
        acts <- dc$decouple(mat=mat, net=net, source=source, 
                            target=target, 
                            weight=mor,
                            methods=methods,
                            consensus=TRUE, min_n=min_n)
    } else {
        stop("No Valid decoupler run type! Try anoter one.")
    }
    
    if (data_to_misc) {
        
        # ----------
        # Append the data in seurat object
        # ----------
        if (run_type %in% c("ulm", "mlm", "ora", "viper", "wmean", "wsum")) {
            seurat_object@misc[["decoupleR"]][[run_type]] <- t(acts[[1]])
        } else if (run_type %in% c("aucell", "gsva", "mdt", "udt")) {
            seurat_object@misc[["decoupleR"]][[run_type]] <- t(acts)    
        } else if (run_type == "decouple") {
            if (include_consensus) {
                for (i in c(methods, "consensus")) {
                    seurat_object@misc[["decoupleR"]][[i]] <- t(acts[[paste0(i, "_estimate")]])
                }
            } else {
                for (i in c(methods)) {
                    seurat_object@misc[["decoupleR"]][[i]] <- t(acts[[paste0(i, "_estimate")]])
                }
            }
        }
    }
    
    return(seurat_object)
    
}

# ----------
# infer cnv specific
# ----------
run_iCNV <- function(countMat, annFile, genePos, refNames) {
    
    # ----------
    # load package
    # ----------
    library(infercnv)
    
    # ----------
    # create the infer_cnv object
    # ----------
    infercnv_obj <- CreateInfercnvObject(
        
        raw_counts_matrix = countMat,
        annotations_file = annFile,
        delim = "\t",
        gene_order_file = genePos,
        ref_group_names = refNames
        
    )
    
    # ----------
    # set the run iCNV function
    # ----------
    infercnv_obj_run <- infercnv::run(
        infercnv_obj,
        cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
        out_dir=out_dir,
        cluster_by_groups=TRUE, 
        plot_steps=FALSE,
        denoise=TRUE,
        HMM=FALSE,
        no_prelim_plot=TRUE,
        output_format = "pdf", # set to 'png' if you want the 'png' output
        # plot_chr_scale = TRUE,
        useRaster = TRUE
    )
    
    return(infercnv_obj_run)
    
}