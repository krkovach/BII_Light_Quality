################################################################################
##### Nearest record between sensors
################################################################################

# From the matching files search for the nearest record between sensors using a 
# target sensor.

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
#' Arguments
#' @param file A data.table product from 01-matching_files.
#' @param index_target A column index of the target sensor.
#' @param index_auxiliary A column index of the auxiliary sensor.

#-------------------------------------------------------------------------------
#' @example

# file <- fread("F:/ligth_quality/Data processing/file-index/FAB_file-ind.txt")
# index_target <- 7
# index_auxiliary <- 6

# test <- nearest_record(file, index_target, index_auxiliary)
# head(test)
# Subset by threshold
# test <- subset(test, abs(time_difference) <= 10)
# head(test)

#-------------------------------------------------------------------------------
#' Function

nearest_record <- function(file, index_target, index_auxiliary) {
  
  #Check for date and time
  if(class(file$dates)[1] != "IDate") {
    file$dates <- as.IDate(file$dates)
  }
  
  if(class(file$time_sequence) != "ITime") {
    file$time_sequence <- as.ITime(file$time_sequence)
  }
  
  #Copy frame
  file_copy <- file

  #Matching function
  file_copy$date_time <- ymd_hms(paste0(file_copy$date, " ", 
                                        file_copy$time_sequence))
  
  #Select sensors
  target <- file_copy[,.SD,.SDcols= c(ncol(file_copy), index_target)]
  target <- na.exclude(target)
  auxiliary <- file_copy[,.SD,.SDcols= c(ncol(file_copy), index_auxiliary)]
  auxiliary <- na.exclude(auxiliary)
  
  #to complete
  complete <- data.frame(date = as.Date("1999-01-01"),
                         time = as.ITime("00:00:00 UTC"),
                         target = 0,
                         auxiliary  = 0,
                         time_difference = 0)
  
  #Loop over file1 rows
  for(i in 1:nrow(target)) {
    
    #Estimate difference
    difference <- difftime(target$date_time[i], auxiliary$date_time, units = "secs")
    
    #Minimum value index
    min_value <- which.min(abs(difference))
    
    #Value
    value <- difference[min_value]
    
    #Fill
    complete[i, 1] <- as.Date(target$date_time[i])
    complete[i, 2] <- as.ITime(target$date_time[i])
    complete[i, 3] <- target[i , .SD, .SDcols= 2]
    complete[i, 4] <- auxiliary[min_value , .SD, .SDcols= 2]
    complete[i, 5] <- value
    
  }
  
  complete <- as.data.table(complete)
  colnames(complete) <- c("date", "time", "target", "auxiliary", "time_difference")
  
  return(complete)
  
}
