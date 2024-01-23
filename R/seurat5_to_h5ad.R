# ----------
# Load library
# ----------
library(Seurat)
library(SeuratDisk)


seurat_to_h5ad <- function(seurat_object, output_directory) {
    
    # -----
    # Remove layers to convert it into the old Seurat 4 format
    # -----
    object <- seurat_object
    object[["RNA3"]] <- as(object = object[["RNA"]], Class = "Assay")
    DefaultAssay(object) <- "RNA3"
    object[["RNA"]] <- NULL
    object <- RenameAssays(object = object, RNA3 = 'RNA')
    
    # -----
    # Save the file
    # -----
    SaveH5Seurat(object, filename = file.path(output_directory, "seurat.h5Seurat"))
    Convert(file.path(output_directory, "seurat.h5Seurat"), dest = "h5ad")
    
}


# ----------
# Set the output directory
# ----------
output_directory <- "your/output/directory/path"

# ----------
# Convert and save the h5ad
# ----------
seurat_to_h5ad(seurat_object, output_directory)
