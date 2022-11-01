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
root_path <- "F:/ligth_quality/New_ASD/Cloquet"

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


stack <- data.table()

for(i in 1:nrow(frame)) {
  
  data <- as.data.frame(fread(paste0(root_path, "/", files[i])))
  data <- c(frame$folder[i], frame$file[i], data[, 2])
  data <- matrix(data, nrow = 1)
  stack <- rbind(stack, data)
  
}

colnames(stack) <- c("date", "file", data$Wavelength)

fwrite(stack, "ASD_IDENT-Cloquet.csv")

