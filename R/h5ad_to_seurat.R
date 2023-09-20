##########
# Load libraries
library(Seurat) 
library(anndata)

##########
# Convert adata to seuart
adata <- read_h5ad("xxx.h5ad")
seurat.object <- CreateSeuratObject(counts = t(as.matrix(adata$layers["counts"])), meta.data = adata$obs)