################################################################################
##### Get path and information for ASD files
################################################################################

# Look and extract information of .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param root_path select the root path of 'Integrating Sphere' folder

# example
root_path <- "/media/antonio/antonio_ssd/ligth_quality/FAB1_2_Irradiance_Measurements/ASD"

#-------------------------------------------------------------------------------
#Function
asd_files <- function(root_path) {
  
  #Search for paths
  files <- list.files(path = paste0(root_path), 
                      pattern = ".asd", 
                      all.files = TRUE,
                      full.names = FALSE, 
                      recursive = TRUE)
  
  #Arrange path in frame
  frame <- data.table(matrix(unlist(strsplit(files, "/")), 
                             nrow= length(files), 
                             byrow=TRUE), stringsAsFactors=FALSE)
  colnames(frame) <- c("folder", "file")
  
  return(frame)
  
}


