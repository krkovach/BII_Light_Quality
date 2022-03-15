################################################################################
##### Function to read SSM files and get metadata
################################################################################

# Look and extract information of .SSM files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param path paths of .SSM files

#-------------------------------------------------------------------------------
#Function
read_ssm <- function(paths) {
  
  #To fill
  to_complete <- data.table()
  to_meta <- data.frame(Wave = NA,
                        Pix = NA,
                        Val = NA,
                        Time = NA)
  
  #Wavelength
  wv <- NA
  
  #loop over paths
  for(i in 1:length(paths)) {
    
    #Metadata
    info <- fread(paths[i], skip= 1, nrows = 1, header = FALSE, sep = " ")
    
    #Fill
    to_meta[i, 1] <- as.numeric(substr(info$V3, 6, 11))
    to_meta[i, 2] <- as.numeric(substr(info$V4, 5, 7))
    
    if(ncol(info) == 15) {
    
      to_meta[i, 3] <- info$V6
      to_meta[i, 4] <- as.numeric(substr(info$V7, 6, 8))
    
    } else {
      
      to_meta[i, 3] <- as.numeric(substr(info$V5, 5, 11))
      to_meta[i, 4] <- as.numeric(substr(info$V6, 6, 8))
    
    }
    
    #Spectra
    spectra <- fread(paths[i], skip= 2)
    
    if(i == 1) {
      wv <- spectra$V1
    }
    
    to_complete <- cbind(to_complete, spectra[, .SD , .SDcols= c(2)])
    
  }
  
  #Arrange
  transp <- t(to_complete)
  transp <- as.data.table(transp)
  colnames(transp) <- as.character(wv)
  
  #Merge meta with spectra 
  frame <- cbind(to_meta, transp)
  
  return(frame)
  
}