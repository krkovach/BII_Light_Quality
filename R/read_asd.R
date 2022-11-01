################################################################################
##### Function to read ASD txt files and get metadata
################################################################################

# Look and extract information of .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param path paths of .asd files

#-------------------------------------------------------------------------------
#Function
read_asd <- function(paths) {
  
  #To fill
  to_complete <- data.table()
  to_meta <- data.frame(date = NA,
                        time = NA)
  
  #Wavelength
  wv <- NA
  
  #loop over paths
  for(i in 1:length(paths)) {
    
    #Metadata
    date <- fread(paths[i], skip= 6, nrows = 1, header = FALSE)
    date <- as.IDate(paste0(substr(date$V3[1], 7, 10), "-",
                            substr(date$V3[1], 4, 5), "-",
                            substr(date$V3[1], 1, 2)))
    
    time <- fread(paths[i], skip= 6, nrows = 1, header = FALSE)
    time <- time$V5
    
    #Fill
    to_meta[i, 1] <- date
    to_meta[i, 2] <- time
    
    #Spectra
    spectra <- fread(paths[i], skip= 39, header = TRUE)
    
    if(i == 1) {
      wv <- spectra$Wavelength
    }
    
    to_complete <- cbind(to_complete, spectra[, .SD , .SDcols= c(2)])
    
  }
  
  #Arrange
  transp <- t(to_complete)
  transp <- as.data.table(transp)
  colnames(transp) <- as.character(wv)
  
  #As date
  to_meta$date <- as.IDate(to_meta$date)
  
  #Merge meta with spectra 
  frame <- cbind(to_meta, transp)
  
  return(frame)
  
}


