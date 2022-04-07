################################################################################
##### Match files by date and time
################################################################################

# Function to match files by date and time

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
#' Arguments
#' @param file1 A data.table file as master.
#' @param file2 A data.table file as slave to anchor with the master.
#' @param threshold A interger describing the threshold in seconds to match with 
#' the master file.
 
file1 <- fread("J:/ligth_quality/Data processing/Pyranometer/FAB_pyranometer.txt")
file2 <- fread("J:/ligth_quality/Data processing/BLK-C/FAB_BLK-C.txt")


#-------------------------------------------------------------------------------
#' Function

matching_files <- function(file1, file2, threshold = 300) {
  
  #Matching function
  date_time_1 <- ymd_hms(paste0(file1$date, " ", file1$time))
  date_time_2 <- ymd_hms(paste0(file2$date, " ", file2$time))
  
  #Anchor vector
  return <- data.table(master = 1:length(date_time_1), 
                       anchor = rep(NA, length(date_time_1)))
  
  #Loop over file1 rows
  for(i in 1:length(anchor)) {
    
    #Estimate difference
    difference <- abs(date_time_1[i] - date_time_2)
    
    #Minimum value index
    min_value <- which.min(difference)
    
    #Value
    value <- as.numeric(difference[min_value], units = "secs")
    
    
    if(value <= threshold) {
      
      return$anchor[i] <- min_value
      
    }
  }
  
  
}


