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
  to_meta <- data.frame(date = NA,
                        time = NA,
                        Wave = NA,
                        Pix = NA,
                        Val = NA,
                        iTime = NA)
  to_meta$date <- as.Date(to_meta$date)
  
  #Wavelength
  wv <- NA
  
  #loop over paths
  for(i in 1:length(paths)) {
    
    #File metadata
    meta_file <- file.mtime(paths[i])
    to_meta[i, 1] <- as.Date(meta_file) #get date
    to_meta[i, 2] <- strftime(meta_file, format="%H:%M:%S") #get time
    
    #Metadata
    info <- fread(paths[i], skip= 1, nrows = 1, header = FALSE, sep = ":")
    
    #Fill
    to_meta[i, 3] <- strsplit(info$V2, " ")[[1]][1] #get wave
    to_meta[i, 4] <- strsplit(info$V3, " ")[[1]][1] #get pix
    to_meta[i, 5] <- strsplit(info$V4, " ")[[1]][1] #get val
    to_meta[i, 6] <- strsplit(info$V5, " ")[[1]][1] #get itime
    
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