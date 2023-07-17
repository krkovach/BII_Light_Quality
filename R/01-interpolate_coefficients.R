################################################################################
#### Interpolate coefficients
################################################################################

# This function interpolate the coefficients thought the time between two sets 
# of measurements.

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(lubridate)
library(imputeTS)

#-------------------------------------------------------------------------------
#' Arguments
#' @param dates A vector of dates to consider.
#' @param time_range A range of times to consider.
#' @param time_span A time span in seconds of the records to consider.
#' @param coefficients A data.table of coefficients.

#-------------------------------------------------------------------------------
#' @example

root_path <- "/media/antonio/antonio_ssd/ligth_quality/Data processing/Coefficients"
coefficients <- fread(paste0(root_path, "/IDENT-Freiburg_coefficients_svc-asd_subday_averaged.txt"), header = TRUE)
time_span <- 1
out_path <- paste0(root_path, "/IDENT-Freiburg_coefficients__svc-asd_interpolation.txt")

interpolate_coefficients(coefficients, time_span, out_path)

#-------------------------------------------------------------------------------
#' Function

interpolate_coefficients <- function(coefficients, time_span, out_path) {
  
  #Dates and times
  coefficients$date <- as.IDate(coefficients$date)
  dates <- unique(coefficients$date)
  coefficients$time <- as.POSIXct(coefficients$time, format="%H:%M:%OS", tz="GMT")
  coefficients$time <- as.ITime(coefficients$time)
  
  #Create sequence of time
  time_range <- as.ITime(range(as.ITime(coefficients$time)))
  time_sequence <- seq(time_range[1], time_range[2], by = 1)
  
  #Complete
  complete <- data.table()
  
  #Loop over dates
  for(i in 1:length(dates)) {
    
    #Complete frame to fill
    frame <- CJ(date = dates[i], time = time_sequence)
    
    #Merge by date
    date_coefficients <- subset(coefficients, date == dates[i])
    empty <- merge(frame, date_coefficients, by = c("date", "time"), all.x = TRUE, all.y = TRUE)
    time_range <- as.ITime(range(date_coefficients$time))
    empty <- subset(empty, time >= time_range[1])
    empty <- subset(empty, time <= time_range[2])
    
    #Loop over columns
    cols <- colnames(empty)
    fill <- as.data.frame(empty)
    
    for(j in 3:length(cols)) {
      
      fill[, j] <- na_interpolation(fill[, j])
      
    }
    
    complete <- rbind(complete, fill)
    
  }
  
  fwrite(complete, out_path, sep = "\t")
  
  return(complete)
  
}


