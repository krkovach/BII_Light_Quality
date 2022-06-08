################################################################################
##### Match frames by date and time
################################################################################

# Function to frames files by date and time

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
#' Arguments
#' @param file1 A data.table file as master.
#' @param file2 A data.table file as slave to anchor with the master.
#' @param threshold A integer describing the threshold in seconds to match with 
#' the master file.

#-------------------------------------------------------------------------------
#' Function

matching <- function(file1, file2, threshold = 1) {
  
  #Matching function
  date_time_1 <- ymd_hms(paste0(file1$date, " ", file1$time))
  date_time_2 <- ymd_hms(paste0(file2$date, " ", file2$time))
  
  #Anchor vector
  complete <- data.table(rows_1 = 1:length(date_time_1), 
                         rows_2 = rep(NA, length(date_time_1)))
  
  #Loop over file1 rows
  for(i in 1:length(date_time_1)) {
    
    #Estimate difference
    difference <- abs(date_time_1[i] - date_time_2)
    
    #Minimum value index
    min_value <- which.min(difference)
    
    #Value
    value <- as.numeric(difference[min_value], units = "secs")
    
    #Threshold logic
    if(value <= threshold) {
      
      complete$rows_2[i] <- min_value
      
    }
  }
  
  #Clean return
  complete <- na.exclude(complete)
  
  #Return frame
  return(complete)
  
}
