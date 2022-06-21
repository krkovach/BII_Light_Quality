################################################################################
#### Plot coefficients
################################################################################

# Plot the coefficients derived from comparisons between sensors.

#-------------------------------------------------------------------------------
#' Libraries

library(data.table)
library(lubridate)
library(ggplot2)


#-------------------------------------------------------------------------------
#' Arguments
#' @param coefficients A frame with the coefficients from 03-get_coefficients.R
#' 

#-------------------------------------------------------------------------------
#' @usage 

#Read coefficients
coeff <- fread("F:/ligth_quality/Data processing/Coefficients/IDENT-Cloquet_coefficients_svc-psm.txt")
coeff <- na.exclude(coeff)
head(coeff[, 1:10])

#Melt wavelength

coeff_melt <- melt(coeff,
                   id.vars = c("date", "time", "target", "reference", "time_difference"),
                   variable.name = "wavelength",
                   value.name = "coefficient")

coeff_melt$wavelength <- as.numeric(as.character(coeff_melt$wavelength))
coeff_melt$time <- as.ITime(coeff_melt$time)

#Get median and 95% intervals
coeff_median <- coeff_melt[, .(median = median(coefficient), sd = sd(coefficient))
                           , by = c("wavelength")]

fwrite(coeff_median, "IDENT-Cloquet_coefficients.csv")

#Plot

ggplot(coeff_median) +
  geom_ribbon(aes(x = wavelength, ymin= (median-sd), ymax = (median+sd)), fill="grey") +
  geom_line(aes(x = wavelength, y = median)) +
  coord_cartesian(ylim=c(0.0, 5)) +
  ylab("Reference / Target") +
  xlab("Wavelength (nm)") +
  theme_bw()

ggplot(coeff_melt) +
  geom_line(aes(x = wavelength, y = coefficient, group = target, colour = date)) +
  coord_cartesian(ylim=c(0.0, 5)) +
  ylab("Reference / Target") +
  xlab("Wavelength (nm)") +
  theme_bw()
  
