################################################################################
##### Get path and information for Integrating Sphere files
################################################################################

# Look and extract information of .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param root_path select the root path of 'Integrating Sphere' folder

# example
root_path <- "/media/antonio/antonio_ssd/ligth_quality/Integrating Sphere"

#-------------------------------------------------------------------------------
#Function
sig_files <- function(root_path) {
  
  #Search for paths
  files <- list.files(path = paste0(root_path), 
                      pattern = ".sig", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "file")
  
  #Get site
  frame[, site := strsplit(file, "[.]")[[1]][1], by = seq_along(1:nrow(frame))]
  frame[, plot := strsplit(file, "[.]")[[1]][2], by = seq_along(1:nrow(frame))]
  frame[, sp := strsplit(file, "[.]")[[1]][3], by = seq_along(1:nrow(frame))]
  frame[, abr1 := strsplit(file, "[.]")[[1]][4], by = seq_along(1:nrow(frame))]
  frame[, abr2 := strsplit(file, "[.]")[[1]][5], by = seq_along(1:nrow(frame))]
  frame$abr <- paste0(frame$abr1, ".", frame$abr2) 
  
  #Clean
  remove <- is.na(frame$abr2)
  
  frame[remove == TRUE, site := NA]
  frame[remove == TRUE, plot := NA]
  frame[remove == TRUE, sp := NA]
  frame[remove == TRUE, abr1 := NA]
  frame[remove == TRUE, abr2 := NA]
  frame[remove == TRUE, abr := NA]
  
  
  #Final to export
  frame <- frame[, c(1, 2, 3, 4, 5, 8)]
  
  return(frame)
  
}


