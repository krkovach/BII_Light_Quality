################################################################################
##### Function to read PSM files and get metadata
################################################################################

# Look and extract information of .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param path paths of .sed files
#' @param column column of interest to extract information. 1, reference irradiance; 
#' 2, target irradiance; 3, targ/refe.

#-------------------------------------------------------------------------------
#Function
psm_files <- function(paths, column = 2) {
  
  #To fill
  to_complete <- data.table()
  to_meta <- data.frame(date = NA,
                        time = NA,
                        type = NA)
  
  #Wavelength
  wv <- NA
  
  #loop over paths
  for(i in 1:length(paths)) {
    
    #Metadata
    time <- fread(paths[i], skip= 7, nrows = 1, header = FALSE)
    time <- time$V2
    
    date <- fread(paths[i], skip= 6, nrows = 1, header = FALSE)
    date <- date$V2
    
    type <- fread(paths[i], skip= 14, nrows = 1, header = FALSE)
    type <- type$V3
    
    #Fill
    to_meta[i, 1] <- date
    to_meta[i, 2] <- time
    to_meta[i, 3] <- type
    
    #Spectra
    spectra <- fread(paths[i], skip= 26)
    
    if(i == 1) {
      wv <- spectra$Wvl
    }
    
    to_complete <- cbind(to_complete, spectra[, .SD , .SDcols= c(column+1)])
  
  }
  
  #Arrange
  transp <- t(to_complete)
  transp <- as.data.table(transp)
  colnames(transp) <- as.character(wv)
  
  #Merge meta with spectra 
  frame <- cbind(to_meta, transp)
  
  return(frame)
  
}


