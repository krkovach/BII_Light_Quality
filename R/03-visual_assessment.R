################################################################################
#### Transmittance visual quality estimation
################################################################################

# Produce several plots of transmittance to visual assess the quality of the 
# measurement.

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(ggplot2)

#-------------------------------------------------------------------------------
#' Arguments
#' @param transmittance A data.table with the transmittance estimations. Output
#' from '03-get_transmittance.R'.
#' @param site_name A character vector to label the outputs.
#' @param output_path A path to storage the plots.

#-------------------------------------------------------------------------------
#' @example

transmittance <- fread("F:/ligth_quality/Data processing/Transmittance/IDENT-Cloquet_transmittance_svc-psm.txt")
meta <- fread("F:/ligth_quality/Data processing/Metadata/IDENT-Cloquet_SVC-locations.txt")
transmittance <- merge(meta, transmittance, all.x = FALSE, all.y = TRUE, by.x = "file", by.y = "target_file")

site_name <- "IDENT-Cloquet"
output_path <- "F:/ligth_quality/Data processing/visual_assessment/IDENT-Cloquet"

visual_assessment(transmittance, site_name, output_path)

#-------------------------------------------------------------------------------
#' Function
visual_assessment <- function(transmittance, site_name, output_path) {
  
  #Read wavelength
  Wavelength <- as.numeric(colnames(transmittance)[10:ncol(transmittance)])
  
  #Files
  files <- nrow(transmittance)
  
  #Create label names
  name_time <- paste0("Site: ",
                      site_name,
                      "  -  Date: ", 
                      transmittance$date,
                      "  -  ", 
                      "Time: ",
                      transmittance$time)
  name_files <- paste0("Target: ",
                       transmittance$target_file,
                       "  -  ",
                       "Reference: ",
                       transmittance$reference_file)
  
  #Loop over rows
  for(i in 1:files) {
    
    Transmittance <- as.numeric(transmittance[i, 10:2158])
    
    aveg <- mean(Transmittance, na.rm = TRUE)
    
    if(is.nan(aveg) == FALSE) {
      
      if(aveg >= 1) {
        label <- "Average transmittance is higher than 1"
      } else {
        label <- ""
      }
      
      plot <- ggplot() +
        geom_line(aes(x = Wavelength, y = Transmittance)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(ylim= c(0, 1.5), xlim = c(min(Wavelength), max(Wavelength))) +
        geom_hline(yintercept= 1, linetype = "dashed", color = "red") +
        theme_light() + 
        labs(
          title = name_time[i],
          subtitle = name_files[i],
          caption = label,
        )
      
      #Provide name
      filename <- paste0(output_path, "/row_", i, ".png")
      
      #Export
      ggsave(plot, filename = filename, device = "png")
      
    }
    
  }
}


