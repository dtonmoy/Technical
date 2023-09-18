##########
# Load libraries
library(ff)

##########
# Set parameters
from <- "/Users/tonmoy/Research/CELLICATE_project/test/directory_move/test/first"
to <- "/Users/tonmoy/Research/CELLICATE_project/test/directory_move/test/second"

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