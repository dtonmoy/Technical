# ----------
# Load library
# ----------
library(Seurat)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/tonmoy/Research/Spatial_proximity_project/env/decoupleR/bin/python")
use_condaenv(condaenv = "decoupleR", conda = "auto", required = NULL)


# ##################################################
# Get the image information from the h5ad
get_image_h5ad <- function(adata, assay="Spatial", slice="slice1") {

    # ----------
    # Fix the image slot
    # ----------
    image <- adata$uns$spatial[[names(adata$uns$spatial)]]$image$lowres
    
    # ----------
    # Fix the scale factors
    # ----------
    scale.factors <- list(spot = adata$uns$spatial[[names(adata$uns$spatial)]]$scalefactors$spot_diameter_fullres,
                          fiducial = adata$uns$spatial[[names(adata$uns$spatial)]]$scalefactors$fiducial_diameter_fullres,
                          hires = adata$uns$spatial[[names(adata$uns$spatial)]]$scalefactors$tissue_hires_scalef,
                          lowres = adata$uns$spatial[[names(adata$uns$spatial)]]$scalefactors$tissue_lowres_scalef)
    
    class(scale.factors) <- "scalefactors"

    # ----------
    # Get the tissue coordinates
    # ----------
    coordinates <- data.frame(tissue = adata$obs$in_tissue,
                              row = adata$obs$array_row,
                              col = adata$obs$array_col,
                              imagerow = adata$obsm$spatial[, 2],
                              imagecol = adata$obsm$spatial[, 1])
    rownames(coordinates) <- adata$obs_names
    
    # ----------
    # Create an `sp` compatible `FOV` instance
    # ----------
    fov <- CreateFOV(
        coordinates[, c("imagerow", "imagecol")],
        type = "centroids",
        radius = scale.factors[["spot"]],
        assay = assay,
        key = Key(slice, quiet = TRUE))
    
    # ----------
    # Build the final `VisiumV2` - essentially just adding `image` and `scale.factors` to the object
    # ----------
    visium.fov <- new(
        Class = "VisiumV2",
        boundaries = fov@boundaries,
        molecules = fov@molecules,
        assay = fov@assay,
        key = fov@key,
        image = image,
        scale.factors = scale.factors
    )
        
    return(visium.fov)
    
}


# ##################################################
# Convert the h5ad to seurat
h5ad_to_seurat <- function(adata) {

    # ----------
    # Check for the assay type and create seurat object
    # ----------
    if ("spatial" %in% names(adata$uns) & "spatial" %in% names(adata$obsm)) {
        seurat_object <- CreateSeuratObject(counts = t(adata$X), assay="Spatial")
        # === Add the image slot # ---------------------------------------------
        seurat_object@images <- list(slice1 = get_image_h5ad(adata))
    } else {
        seurat_object <- CreateSeuratObject(counts = t(adata$X), assay="RNA")
    }
    
    # ----------
    # Add the metadata information 
    # ----------
    seurat_object <- AddMetaData(seurat_object, metadata = adata$obs)
    
    # ----------
    # Add other layers into the misc slot
    # ----------
    seurat_object@misc <- list(
        obsm = adata$obsm,
        obsm = adata$obsp,
        var = adata$var,
        varm = adata$varm,
        varp = adata$varp,
        uns = adata$uns
    )
    
    # =========================
    # EMBEDDINGS
    # =========================
    # === Check if the umap layer is present or not # --------------------------
    if (!is.null(adata$obsm[["X_umap"]])) {
        
        # ----------
        # Add cell embedding
        # ----------
        embedding <- as.data.frame(adata$obsm["X_umap"])
        rownames(embedding) <- adata$obs_names
        colnames(embedding) <- c("umap_1", "umap_2")
        embedding <- as.matrix(embedding)
        seurat_object@reductions[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")
    }
    
    # === Check if the PCA layer is present or not # ---------------------------
    if (!is.null(adata$obsm[["X_pca"]])) {
        
        # ----------
        # Add cell embedding
        # ----------
        embedding <- as.data.frame(adata$obsm["X_pca"])
        rownames(embedding) <- adata$obs_names
        colnames(embedding) <- paste0("PC_", 1:ncol(embedding))
        embedding <- as.matrix(embedding)
        seurat_object@reductions[["pca"]] <- CreateDimReducObject(embedding, key = "pca_")
        
        # ----------
        # Add the feature loading
        # ----------
        seurat_object@reductions[["pca"]]@feature.loadings <- adata$varm[["PCs"]]
    }

    return(seurat_object)
    
}


# ##################################################
# Check if the conversion is working or not (h5ad to seurat)
# ##################################################
spatial_adata <- read_h5ad("/Users/tonmoy/Downloads/UKF334_T_ST.h5ad")
sc_bla <- read_h5ad("/Users/tonmoy/Downloads/Dummy_object.h5ad")
seurat_object_converted <- h5ad_to_seurat(sc_bla) # spatial_adata/sc_bla --> works FINE


# ##################################################
# Convert the seurat to h5ad
seurat_to_h5ad <- function(seurat_object) {

    if (all(c("obsm", "var", "varm", "varp") %in% names(seurat_object@misc))) {
        
        # ----------
        # Check for the assay type and create adata object
        # ----------
        if (names(seurat_object@assays) == "Spatial") {
            # == Spatial assay # -----------------------------------------------
            key <- "Spatial"
        } else if (names(seurat_object@assays) == "RNA") {
            # == RNA assay # ---------------------------------------------------
            key <- "RNA"
        }
        
        # ----------
        # Check for the assay type and create adata object
        # ----------
        adata <- AnnData(X = t(seurat_object@assays[[key]]@layers[["counts"]]),
                         obs = seurat_object@meta.data,
                         obsm = seurat_object@misc[["obsm"]],
                         var = seurat_object@misc[["var"]],
                         varm = seurat_object@misc[["varm"]],
                         varp = seurat_object@misc[["varp"]])

    } else {
        
        # ----------
        # Create var dataframe to pass it on the adata object
        # ----------
        var <- data.frame(rownames(seurat_object))
        colnames(var) <- c("gene")
        rownames(var) <- var$gene

        # ----------
        # Prepare the dimentionality reduction slots
        # ----------
        # == PCA -- cell embedding # -------------------------------------------
        if ("pca" %in% names(seurat_object@reductions)) {
            pca <- seurat_object@reductions[["pca"]]@cell.embeddings
            rownames(pca) <- NULL
            colnames(pca) <- NULL
        } else {
            pca <- NULL
        }

        # == PCA -- feature loading # ------------------------------------------
        if ("pca" %in% names(seurat_object@reductions)) {

            # =========================
            # Because seurat can store the informative sub setted features rather
            # than the full set of features, we need to fix this issue for the scanpy
            # =========================
            if (!is.null(seurat_object@reductions[["pca"]]@feature.loadings)) {
                
                # == Create a blank matrix with the same dimensions # ----------
                pca_fl <- matrix(0, nrow = nrow(seurat_object), ncol = ncol(seurat_object@reductions[["pca"]]@feature.loadings))
                
                # == Fix the row names at first # ------------------------------
                rownames(pca_fl) <- rownames(seurat_object)
                
                # == Replace the rows with the same PCA loading # --------------
                pca_fl[rownames(seurat_object@reductions[["pca"]]@feature.loadings), ] <- seurat_object@reductions[["pca"]]@feature.loadings
                
                rownames(pca_fl) <- NULL
                colnames(pca_fl) <- NULL
            
            }
            
        } else {
            pca_fl <- NULL
        }
        
        # == UMAP -- cell embedding # ------------------------------------------
        if ("umap" %in% names(seurat_object@reductions)) {
            umap <- seurat_object@reductions[["umap"]]@cell.embeddings
            rownames(umap) <- NULL
            colnames(umap) <- NULL
        } else {
            umap <- NULL
        }

        # ----------
        # Prepare the assay names and version ("v1" or "v2")
        # ----------
        if (names(seurat_object@assays) == "Spatial") {
            # == Spatial assay # -----------------------------------------------
            key <- "Spatial"
            # == Check the version of the spatial assay # ----------------------
            if (class(seurat_object@images[["slice1"]])[1] == "VisiumV1") {
                assay_version <- "v1"
            } else if (class(seurat_object@images[["slice1"]])[1] == "VisiumV2") {
                assay_version <- "v2"
            } else {
                stop("The spatial assay version is not supported for now!!")
            }
            
        } else if (names(seurat_object@assays) == "RNA") {
            # == RNA assay # ---------------------------------------------------
            key <- "RNA"
        } else {
            stop("The seurat object does not have the required slots to convert it into h5ad")
        }
            
        
        # =========================
        # Create the image slot structure for adata if the assay is spatial
        # =========================
        if (key == "Spatial") {
            
            image_slot_structure <- list(
                spatial_sample = list(
                    images = list(
                        lowres = seurat_object@images[["slice1"]]@image
                    ),
                    scalefactors = list(
                        fiducial_diameter_fullres = seurat_object@images[["slice1"]]@scale.factors[["fiducial"]],
                        spot_diameter_fullres = seurat_object@images[["slice1"]]@scale.factors[["spot"]],
                        tissue_hires_scalef = seurat_object@images[["slice1"]]@scale.factors[["hires"]],
                        tissue_lowres_scalef = seurat_object@images[["slice1"]]@scale.factors[["lowres"]]
                    )
                )
            )
            
        }
        
        # =========================
        # Create the adata object
        # =========================
        adata <- AnnData(X = t(seurat_object@assays[[key]]@layers[["counts"]]),
                         obs = seurat_object@meta.data,
                         var = var)

        # ----------
        # Add spatial information to the adata object
        # ----------
        if (key == "Spatial") {
            adata$uns$spatial <- image_slot_structure
            
            # == Add the spatial coordinates # ---------------------------------
            if (assay_version == "v1") {
                
                adata$obs <- cbind(seurat_object@images[["slice1"]]@coordinates[["tissue"]],
                                   seurat_object@images[["slice1"]]@coordinates[["row"]],
                                   seurat_object@images[["slice1"]]@coordinates[["col"]],
                                   adata$obs)
                # ==========
                # Fix the column names
                # ==========
                colnames(adata$obs)[colnames(adata$obs) == 'seurat_object@images[["slice1"]]@coordinates[["tissue"]]'] <- "in_tissue"
                colnames(adata$obs)[colnames(adata$obs) == 'seurat_object@images[["slice1"]]@coordinates[["row"]]'] <- "array_row"
                colnames(adata$obs)[colnames(adata$obs) == 'seurat_object@images[["slice1"]]@coordinates[["col"]]'] <- "array_col"
                
                # =========================
                # Add spatial info in the obsm
                # =========================
                adata$obsm$spatial <- cbind(seurat_object@images[["slice1"]]@coordinates[["imagecol"]],
                                            seurat_object@images[["slice1"]]@coordinates[["imagerow"]])
                colnames(adata$obsm$spatial) <- NULL
                
            } else if (assay_version == "v2") {
                
                # ==========
                # Append the in_tissue column
                # ==========
                adata$obs$in_tissue <- 1
                
                adata$obs <- cbind(data.frame(seurat_object@images[["slice1"]]@boundaries[["centroids"]]@coords)$x,
                                   data.frame(seurat_object@images[["slice1"]]@boundaries[["centroids"]]@coords)$y,
                                   adata$obs)
                
                # =========================
                # Add spatial info in the obsm
                # =========================
                adata$obsm$spatial <- seurat_object@images[["slice1"]]@boundaries[["centroids"]]@coords
                colnames(adata$obsm$spatial) <- NULL
                
                # ==========
                # Fix the column names
                # ==========
                colnames(adata$obs)[colnames(adata$obs) == 'data.frame(seurat_object@images[["slice1"]]@boundaries[["centroids"]]@coords)$x'] <- "array_row"
                colnames(adata$obs)[colnames(adata$obs) == 'data.frame(seurat_object@images[["slice1"]]@boundaries[["centroids"]]@coords)$y'] <- "array_col"
                
                
            }
        }
        
        # =========================
        # Add layers
        # =========================
        # == Count data # ------------------------------------------------------
        if (class(seurat_object@assays[[key]]@layers[["counts"]]) == "dgCMatrix") {
            adata$layers[["counts"]] <- t(seurat_object@assays[[key]]@layers[["counts"]])
        }
        
        # == log2 normalized data # --------------------------------------------
        if (class(seurat_object@assays[[key]]@layers[["data"]]) == "dgCMatrix") {
            adata$layers[["log2norm_counts"]] <- t(seurat_object@assays[[key]]@layers[["counts"]])
        }
        
        # == Scaled data # -----------------------------------------------------
        if ((class(seurat_object@assays[[key]]@layers[["scale.data"]]) == "matrix")[1]) {
            scaled_data <- t(seurat_object@assays[[key]]@layers[["scale.data"]])
            colnames(scaled_data) <- NULL
            adata$layers[["scaled"]] <- scaled_data
        }
            
        # ----------
        # Check if the following slots are present or not
        # If present, then add them to the adata object
        # ----------
        # == PCA cell embedding # ----------------------------------------------
        if (!is.null(pca)) {
            adata$obsm$X_pca <- pca
        }
        
        # == UMAP cell embedding # ---------------------------------------------
        if (!is.null(umap)) {
            adata$obsm$X_umap <- umap
        }
        
        # == PCA feature loading # ---------------------------------------------
        if (!is.null(pca_fl)) {
            adata$varm$PCs <- pca_fl
        }
        
    } 
    
    return(adata)

}


# ##################################################
# Check if the conversion is working or not (seurat to h5ad)
# ##################################################
lalala <- readRDS("/Users/tonmoy/Research/Spatial_proximity_project/data/all_seurat/seurat_list_processed.rds")
lololo <- seurat_to_h5ad(lalala[[1]]) # --> works fine for v1 spatial

lipili <- Load10X_Spatial("/Users/tonmoy/Research/Spatial_proximity_project/data/LGG/LGG/IDHm_BWH23_oligo/")
lopolo <- seurat_to_h5ad(lipili) # --> works fine for v2 spatial

lamala <- readRDS("/Users/tonmoy/Research/GBM_3/data/ALA_count/ala_seurat_filtered_final.rds")
lomolo <- seurat_to_h5ad(lamala) # --> works fine for sc

# ----------
# Save the results
# ----------
write_h5ad(
    lomolo,
    "/Users/tonmoy/Research/GBM_3/test/seurat_adata_conversion/sc.h5ad",
    compression = NULL,
    compression_opts = NULL,
    as_dense = list()
)







