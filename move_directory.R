##########
# Load libraries
library(ff)

##########
# Set parameters
from <- "/Directory/containing/the/folders/that/are/to/be/moved"
to <- "/Directory/where/the/folders/will/be/moved"

##########
# Get the directory lists
dirLists <- list.dirs(from, recursive = FALSE)

for (i in dirLists) {
    
    #-----
    # get the new directory paths
    #-----
    newDir <- file.path(to, basename(i))
    
    #-----
    # create new directory to move files
    #-----
    dir.create(newDir, recursive = TRUE, showWarnings = TRUE)
    
    #-----
    # move the files
    #-----
    file.move(i, newDir)
    
}