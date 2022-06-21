################################################################################
#### Get coefficients of difference between sensors
################################################################################

# Compute the coefficients of difference between sensors

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
#' Arguments
#' @param nr_file A frame from the 02-nearest_record
#' @param target A data.table product from 01-matching_files.
#' @param target_col A data.table product from 01-matching_files.
#' @param reference A data.table product from 01-matching_files.
#' @param reference_col A data.table product from 01-matching_files.
#' @param metadata A frame with the target metadata.
#' 
#-------------------------------------------------------------------------------
#' @example

file <- fread("/media/antonio/antonio_ssd/ligth_quality/Data processing/file-index/FAB_file-ind.txt")
index_target <- 7
index_reference <- 6

nr_file <- nearest_record(file, index_target, index_reference)
head(nr_file)
# Subset by threshold
nr_file <- subset(nr_file, abs(time_difference) <= 40)
head(nr_file)

target <- fread("/media/antonio/antonio_ssd/ligth_quality/Data processing/SVC/FAB_svc.txt")
target_col <- 6:2223
reference <- fread("/media/antonio/antonio_ssd/ligth_quality/Data processing/PSM/FAB_psm.txt")
reference_col <- 6:2156

metadata <- fread("/media/antonio/antonio_ssd/ligth_quality/Data processing/Metadata/FAB_SVC-locations.txt")

# Test function with out calibration file
coeff <- get_coefficients(nr_file, target, target_col, reference, reference_col, metadata)
head(coeff[, 1:10])
fwrite(coeff, 
       "/media/antonio/antonio_ssd/ligth_quality/Data processing/Coefficients/FAB_coefficients_svc-psm.txt", 
       sep = "\t")

#-------------------------------------------------------------------------------
#' Function

get_coefficients <- function(nr_file, 
                             target, 
                             target_col, 
                             reference, 
                             reference_col, 
                             metadata) {
  
  #Get just the open observations
  open <- merge(target, metadata, by = "file", all.x = TRUE, all.y = FALSE)
  open <- open$ID == "open" & is.na(open$ID) == FALSE
  observations <- 1:length(open)
  observations <- observations[open == TRUE]
  
  #Select spectra
  target_spect <- target[, .SD , .SDcols = target_col] 
  reference_spect <- reference[, .SD , .SDcols = reference_col] 
  
  #Get spectral ranges
  target_wv <- as.numeric(colnames(target_spect)) #283 988
  reference_wv <- as.numeric(colnames(reference_spect)) #350 1523
  
  #Min
  twv_match <- which(target_wv%in%reference_wv)
  rwv_match <- which(reference_wv%in%target_wv)
  
  #Match spectral range
  target_match <- target_spect[, .SD , .SDcols = twv_match]
  reference_match <- reference_spect[, .SD , .SDcols = rwv_match]
  
  #Get wv
  twv <- as.numeric(colnames(target_match))
  rwv <- as.numeric(colnames(reference_match))
  
  if(all.equal(twv, rwv) != TRUE) {
    stop("It is likely that the spectral resolution between sensors does not match")
  }
  
  #Empty frame to fill
  frame <- as.data.frame(target_match[0, ])
  meta <- data.table()
  
  for(i in 1:length(observations)) {
    
    #Select rows
    files <- subset(nr_file, target == observations[i])
    
    #Get spectra by index
    inside <- as.numeric(target_match[files$target[1], ])
    outise <- as.numeric(reference_match[files$reference[1], ])
    
    #Estimate coefficients
    coefficients <- outise/inside
    
    #Fill frame
    frame <- rbind(frame, matrix(coefficients, nrow = 1)) 
    
    #Fill files
    files$target[1] <- target$file[files$target[1]]
    files$reference[1] <- reference$file[files$reference[1]]
    meta <- rbind(meta, files)
    
  }
  
  #Correct frame features
  frame <- as.data.table(frame)
  colnames(frame) <- colnames(target_match)
  
  #Add index
  frame <- cbind(meta, frame)
  
  return(frame)
  
}
