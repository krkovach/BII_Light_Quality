################################################################################
#### Transmittance estimation
################################################################################

# From the indexes of nearest record between sensors estimates the transmittance 
# spectra from target and reference spectrometers.

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
# Load source

source("R/02-nearest_record.R")

#-------------------------------------------------------------------------------
#' Arguments
#' @param nr_file A frame from the 02-nearest_record
#' @param target A data.table product from 01-matching_files.
#' @param target_col A data.table product from 01-matching_files.
#' @param reference A data.table product from 01-matching_files.
#' @param reference_col A data.table product from 01-matching_files.
#' @param intercalibration Intercalibration coefficients between the sensors. If NULL,
#' it assumes to be 1.

#-------------------------------------------------------------------------------
#' @example

#Get index file
 file <- fread("F:/ligth_quality/Data processing/file-index/FAB_file-ind.txt")
 index_target <- 7
 index_reference <- 6

#Get nearest record
 nr_file <- nearest_record(file, index_target, index_reference)
 head(nr_file)
 nr_file <- subset(nr_file, abs(time_difference) <= 40)
 head(nr_file)

#Get target sensor
 target <- fread("F:/ligth_quality/Data processing/SVC/FAB_svc.txt")
 target_col <- 6:2223
 
#Get reference sensor
 reference <- fread("F:/ligth_quality/Data processing/PSM/FAB_psm.txt")
 reference_col <- 6:2156

#Get intercalibration file
 intercalibration <- fread("F:/ligth_quality/Data processing/Coefficients/FAB_coefficients.txt",
                           header = TRUE)

# Test function with out calibration file
 trans <- get_transmittance(nr_file, target, target_col, reference, reference_col, intercalibration)
 head(trans[, 1:10])
 fwrite(trans, 
        "F:/ligth_quality/Data processing/Transmittance/FAB_transmittance_svc-psm.txt", 
        sep = "\t")

#-------------------------------------------------------------------------------
#' Function

get_transmittance <- function(nr_file, 
                              target, 
                              target_col, 
                              reference, 
                              reference_col, 
                              intercalibration = NULL) {
 
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
  
  #Add target and reference files
  nr_file$target_file <- NA
  nr_file$reference_file <- NA
  
  #Empty frame to fill
  frame <- as.data.frame(target_match[0, ])
  
  if(is.null(intercalibration) != TRUE) {
    
    intercal_vec <- as.numeric(intercalibration[1, ])
    
  }
  
  for(i in 1:nrow(nr_file)) {
    
    #Get spectra by index
    inside <- as.numeric(target_match[nr_file$target[i], ])
    outise <- as.numeric(reference_match[nr_file$reference[i], ])
    
    #Estimate transmittance 
    transm <- inside/outise
    
    if(is.null(intercalibration) != TRUE) {
      
      transm <- transm*intercal_vec
    
    }
    
    #Fill frame
    frame <- rbind(frame, matrix(transm, nrow = 1)) 
    
    #Fill files
    nr_file$target_file[i] <- target$file[nr_file$target[i]]
    nr_file$reference_file[i] <- reference$file[nr_file$reference[i]]
    
  }
  
  #Correct frame features
  frame <- as.data.table(frame)
  colnames(frame) <- colnames(target_match)
  
  #Add index
  frame <- cbind(nr_file, frame)
  
  return(frame)
  
}
