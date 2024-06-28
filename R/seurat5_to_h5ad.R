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

spatial_adata <- read_h5ad("/Users/tonmoy/Downloads/UKF334_T_ST.h5ad")
sc_bla <- read_h5ad("/Users/tonmoy/Downloads/Dummy_object.h5ad")
seurat_object_converted <- h5ad_to_seurat(sc_bla) # spatial_adata/sc_bla --> works FINE


# ##################################################
# Get the image information from the seurat
sc_seurat <- readRDS("/Users/tonmoy/Research/GBM_3/data/ALA_count/ala_seurat_filtered_final.rds")


# ##################################################
# Convert the seurat to h5ad
seurat_to_h5ad <- function(seurat_object) {

    if (all(c("obsm", "var", "varm", "varp") %in% names(seurat_object@misc))) {
        
        # ----------
        # Check for the assay type and create adata object
        # ----------
        if (names(seurat_object@assays) == "Spatial") {
            
            # == Spatial assay # -----------------------------------------------
            adata <- AnnData(X = t(seurat_object@assays[["Spatial"]]@layers[["counts"]]),
                             obs = seurat_object@meta.data,
                             obsm = seurat_object@misc[["obsm"]],
                             var = seurat_object@misc[["var"]],
                             varm = seurat_object@misc[["varm"]],
                             varp = seurat_object@misc[["varp"]])
        } else if (names(seurat_object@assays) == "RNA") {
            
            # == RNA assay # ---------------------------------------------------
            adata <- AnnData(X = t(seurat_object@assays[["RNA"]]@layers[["counts"]]),
                             obs = seurat_object@meta.data,
                             obsm = seurat_object@misc[["obsm"]],
                             var = seurat_object@misc[["var"]],
                             varm = seurat_object@misc[["varm"]],
                             varp = seurat_object@misc[["varp"]])
        }

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
            if (!is.null(pca_fl)) {
                
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

        if (names(seurat_object@assays) == "Spatial") {
            
            # =========================
            # Create the image slot structure for adata
            # =========================
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
            
            # =========================
            # Spatial assay
            # =========================
            adata <- AnnData(X = t(seurat_object@assays[["Spatial"]]@layers[["counts"]]),
                             obs = seurat_object@meta.data,
                             var = var)
            
            # =========================
            # Add layers
            # =========================
            # == Count data # -----------------------------------------------------
            if (class(seurat_object@assays[["Spatial"]]@layers[["counts"]]) == "dgCMatrix") {
                adata$layers[["counts"]] <- t(seurat_object@assays[["Spatial"]]@layers[["counts"]])
            }
            # == log2 normalized data # --------------------------------------------
            if (class(seurat_object@assays[["Spatial"]]@layers[["data"]]) == "dgCMatrix") {
                adata$layers[["log2norm_counts"]] <- t(seurat_object@assays[["Spatial"]]@layers[["counts"]])
            }
            # == scaled data # -----------------------------------------------------
            if ((class(seurat_object@assays[["Spatial"]]@layers[["scale.data"]]) == "matrix")[1]) {
                scaled_data <- t(seurat_object@assays[["Spatial"]]@layers[["scale.data"]])
                colnames(scaled_data) <- NULL
                adata$layers[["scaled"]] <- scaled_data
            }
            
        } else if (names(seurat_object@assays) == "RNA") {
            
            # =========================
            # RNA assay
            # =========================
            adata <- AnnData(X = t(seurat_object@assays[["RNA"]]@layers[["counts"]]),
                             obs = seurat_object@meta.data,
                             var = var)

            # =========================
            # Add layers
            # =========================
            # == Count data # -----------------------------------------------------
            if (class(seurat_object@assays[["RNA"]]@layers[["counts"]]) == "dgCMatrix") {
                adata$layers[["counts"]] <- t(seurat_object@assays[["RNA"]]@layers[["counts"]])
            }
            # == log2 normalized data # --------------------------------------------
            if (class(seurat_object@assays[["RNA"]]@layers[["data"]]) == "dgCMatrix") {
                adata$layers[["log2norm_counts"]] <- t(seurat_object@assays[["RNA"]]@layers[["counts"]])
            }
            # == scaled data # -----------------------------------------------------
            if ((class(seurat_object@assays[["RNA"]]@layers[["scale.data"]]) == "matrix")[1]) {
                scaled_data <- t(seurat_object@assays[["RNA"]]@layers[["scale.data"]])
                colnames(scaled_data) <- NULL
                adata$layers[["scaled"]] <- scaled_data
                
        } else {
            stop("The seurat object does not have the required slots to convert it into h5ad")
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
        
        # == PCA feature loading # ---------------------------------------------
        if (!is.null(pca_fl)) {
            adata$varm$PCs <- pca_fl
        }
            
        }
        
    } 
    
    return(adata)

    }
}

dim(sc_bla$varm$PCs)

# 
# labla <- seurat_list[[1]]
balala <- seurat_to_h5ad(sc_seurat)

