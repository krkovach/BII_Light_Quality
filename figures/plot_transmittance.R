################################################################################
#### Plot transmittance
################################################################################

# Plot the coefficients derived from comparisons between sensors.

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(lubridate)
library(ggplot2)

#-------------------------------------------------------------------------------
#' Arguments
#' @param transmittance A frame with the coefficients from 03-get_coefficients.R
#' 

#-------------------------------------------------------------------------------
#' @usage 

#Read coefficients
trans <- fread("F:/ligth_quality/Data processing/Transmittance/IDENT-Cloquet_transmittance_svc-psm.txt")
trans <- na.exclude(trans)
head(trans[, 1:10])
trans <- trans[, c(1, 2, 5, 6, 7, 8:2158)]

#Melt wavelength

coeff_melt <- melt(trans,
                   id.vars = c("date", "time", "target_file", "reference_file", "time_difference"),
                   variable.name = "wavelength",
                   value.name = "transmittance")

coeff_melt$wavelength <- as.numeric(as.character(coeff_melt$wavelength))
coeff_melt$time <- as.ITime(coeff_melt$time)

#Get median and 95% intervals
ggplot(coeff_melt) +
  geom_line(aes(x = wavelength, y = transmittance, group = target_file, colour = date)) +
  coord_cartesian(ylim=c(0.0, 5)) +
  ylab("Transmittance") +
  xlab("Wavelength (nm)") +
  theme_bw()

