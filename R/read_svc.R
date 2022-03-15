################################################################################
##### Function to read SVC files and get metadata
################################################################################

# Look and extract information of .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
#Arguments

#' @param path paths of .sed files
#' @param column column of interest to extract information. 1, reference irradiance; 
#' @param new_bands an interger to describe the new bands to resample
#' 2, target irradiance; 3, targ/refe.

#-------------------------------------------------------------------------------
#Function
read_svc <- function(paths, column = 2, new_bands = 338:1990) {
  
  #To fill
  to_complete <- data.table()
  to_meta <- data.frame(date = NA,
                        time = NA,
                        type = NA)
  
  #loop over paths
  for(i in 1:length(paths)) {
    
    #Metadata
    datetime <- fread(paths[i], skip= 17, nrows = 1, header = FALSE)$V2
    datetime <- as.POSIXct(datetime, tz = "", "%m/%d/%Y %H:%M:%S")
    
    date <- as.Date(datetime)
    time <- format(datetime, tz = "", "%H:%M:%S")
    type <- fread(paths[i], skip= 16, nrows = 1, header = FALSE)$V2
    
    #Fill
    to_meta[i, 1] <- as.character(date)
    to_meta[i, 2] <- time
    to_meta[i, 3] <- type
    
    #Spectra
    spectra <- fread(paths[i], skip= 30)
    spectra <- spectra[, .SD , .SDcols= c(1, column+1)]
    colnames(spectra) <- c("wv", "value")
    
    #Create and predict a spline function
    resample_function <- smooth.spline(x = spectra$wv,
                                       y = spectra$value,
                                       spar= 0.05)
    
    spectra_smoothed <- predict(resample_function, new_bands)
    
    #Fill
    to_complete <- cbind(to_complete, spectra_smoothed$y)
    
  }
  
  #Arrange
  transp <- t(to_complete)
  transp <- as.data.table(transp)
  colnames(transp) <- as.character(new_bands)
  
  #Merge meta with spectra 
  frame <- cbind(to_meta, transp)
  
  return(frame)
  
}


