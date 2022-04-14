################################################################################
##### Match files by date and time
################################################################################

# Function to several frames by date and time

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(lubridate)

#-------------------------------------------------------------------------------
# Load source code

source("R/matching.R")

#-------------------------------------------------------------------------------
#' Arguments
#' @param dates A vector of dates to consider.
#' @param time_range A range of times to consider.
#' @param time_span A time span in seconds of the records to consider.
#' @param list_files A data.table file as slave to anchor with the master.
#' @param threshold A interger describing the threshold in seconds to match with 
#' the master file.

dates <- c(as.Date("2021-07-19"), as.Date("2021-07-20"))
time_range <- as.POSIXct(c("08:00:00", "17:00:00"), format="%H:%M:%OS", tz="GMT")
time_span <- 1
list_files <- 
  
#-------------------------------------------------------------------------------
#' Function

matching_files <- function(dates, time_range, time_span, list_files, threshold = 0.99) {
  
  #Create sequence of time
  as_time <- time_range <- as.POSIXct(time_range, format="%H:%M:%OS", tz="GMT")
  time_sequence <- seq(time_range[1], time_range[2], by = "1 secs")
  time_sequence <- as.ITime(time_sequence)
  
  #Complete frame to fill
  complete <- CJ(dates, time_sequence)
  
  #N of files to match
  n_files <- length(list_files)
  
  for(i in 1:n_files) {
    
    file <- matching(complete, list_files[[i]], threshold)
    
    new <- complete[file$rows_1, 1:2]
    new$index <- file$rows_2
    
    complete <- merge(complete, new, all.x = TRUE)
    colnames(complete)[(2+i)] <- names(list_files[i])
    
  }
  
  #Check for all NA
  na_logict <- apply(complete[, .SD, .SDcols = c(3:ncol(complete))], 
                     1, 
                     function(x) all(is.na(x)))
  
  #Just indices
  complete <- complete[na_logict == FALSE]
  
  return(complete)
  
}