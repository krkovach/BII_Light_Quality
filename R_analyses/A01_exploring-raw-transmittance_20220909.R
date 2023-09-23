# Script to explore canopy transmittance
# Plot canopy transmittance values
# Identify outliers/readings to remove

#_______________________________________________________________________________
# Install packages:
library(data.table)
library(stringr)

#_______________________________________________________________________________
# Input files:
  # transmittance SVC-ASD
trans_fab_sa_ <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-01-05_version/FAB_transmittance_svc-asd.txt")
trans_cl_sa_ <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-01-05_version/IDENT-Cloquet_transmittance_svc-asd.txt")
trans_fr_sa_ <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-01-05_version/IDENT-Freiburg_transmittance_svc-asd.txt")

cl_meta <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/IDENT-Cloquet_SVC-locations.txt")
cl_meta2 <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Cloquet_plots-composition-diversity.csv")

fab_meta <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/FAB_SVC-locations.txt")
fab_meta2 <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/FAB1&2_Composition of plots.csv")

fr_meta <- fread("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/IDENT-Freiburg_SVC-locations.txt")
fr_meta2 <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Freiburg_plot-lookup.csv")

colnames(cl_meta)[1] <- c("target_file") # renaming column to match data
colnames(fab_meta)[1] <- c("target_file") # renaming column to match data
colnames(fr_meta)[1] <- c("target_file") # renaming column to match data

cl_meta$ID <- substring(cl_meta$ID,9) # removing site name from ID and keeping rest
# fab_meta$ID <- substring(fab_meta$ID,6) # removing site name from ID and keeping rest
fr_meta$ID <- substring(fr_meta$ID,4) # removing site name from ID and keeping rest


trans_cl_sa <- merge(x=cl_meta, y=trans_cl_sa_, by="target_file", all.y=T, all.x=F) # merging -- CHECK IF NEED TO SET SORT TO FALSE
trans_fab_sa <- merge(x=fab_meta, y=trans_fab_sa_, by="target_file", all.y=T, all.x=F)
trans_fr_sa <- merge(x=fr_meta, y=trans_fr_sa_, by="target_file", all.y=T, all.x=F)

trans_cl_sa <- trans_cl_sa[!is.na(trans_cl_sa$ID),] # removing rows with NAs for plot ID
trans_cl_sa <- trans_cl_sa[which(trans_cl_sa$ID!="NA"),]
trans_fab_sa <- trans_fab_sa[!is.na(trans_fab_sa$ID),] 
trans_fab_sa <- trans_fab_sa[which(trans_fab_sa$ID!="FAB1_NA"&trans_fab_sa$ID!="FAB2_NA"),] 
trans_fr_sa <- trans_fr_sa[!is.na(trans_fr_sa$ID),] 
trans_fr_sa <- trans_fr_sa[which(trans_fr_sa$ID!="NA"),] 

table(trans_cl_sa$ID)
table(trans_fab_sa$ID)
table(trans_fr_sa$ID)


#_______________________________________________________________________________
# Unit vector normalise transmittance values
# First exclude noisy and unreliable bands

col_index <- colnames(trans_cl_sa)%in%c(391:1350, 1451:1800, 1981:2200) # c(350:390, 1350:1450, 1800:1980, 2200:2500))   # columns that want to keep (following Hovi and Rautiainen 2020)
trans_cl_sa_trimmed <- subset(trans_cl_sa,, col_index)
vn <- sqrt(apply(trans_cl_sa_trimmed^2, 1, sum))
trans_cl_sa_vn_ <- trans_cl_sa_trimmed/vn
trans_cl_sa_vn <- cbind(trans_cl_sa[,c(1,2,6,7,10)], trans_cl_sa_vn_) # adding back key metadata

col_index <- colnames(trans_fr_sa)%in%c(391:1350, 1451:1800, 1981:2200) # c(350:390, 1350:1450, 1800:1980, 2200:2500))   # columns that want to keep 
trans_fr_sa_trimmed <- subset(trans_fr_sa, , col_index)
vn <- sqrt(apply(trans_fr_sa_trimmed^2, 1, sum))
trans_fr_sa_vn_ <- trans_fr_sa_trimmed/vn
trans_fr_sa_vn <- cbind(trans_fr_sa[,c(1,2,6,7,10)], trans_fr_sa_vn_) # adding back key metadata

col_index <- colnames(trans_fab_sa)%in%c(391:1350, 1451:1800, 1981:2200) # c(350:390, 1350:1450, 1800:1980, 2200:2500))   # columns that want to keep 
#trans_fab_sa_trimmed <- trans_fab_sa[,..col_index] # another option for subsetting, both work
trans_fab_sa_trimmed <- subset(trans_fab_sa, subset=T, select=col_index)
vn <- sqrt(apply(trans_fab_sa_trimmed^2, 1, sum))
trans_fab_sa_vn_ <- trans_fab_sa_trimmed/vn
trans_fab_sa_vn <- cbind(trans_fab_sa[,c(1,2,6,7,10)], trans_fab_sa_vn_) # adding back key metadata


#_______________________________________________________________________________
# Plotting to check and identify outliers - vector normalized (to check for difference in shape)

# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Cloquet_plots_VN_2023-09-19.pdf")
for(i in 1:length(unique(trans_cl_sa_vn$ID))){
  rns <- which(trans_cl_sa_vn$ID==unique(trans_cl_sa_vn$ID)[i])
  
  plot(x=as.numeric(colnames(trans_cl_sa_vn)[6:1535]), y=trans_cl_sa_vn[rns[1],6:1535], ylim=c(0,0.1), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=as.numeric(colnames(trans_cl_sa_vn)[6:1535]), y=trans_cl_sa_vn[rns[j],6:1535], col=color_palette[j], lwd=2)
  }
  
  title(paste(trans_cl_sa_vn$ID[rns[1]], " - vector normalized", "\n(", trans_cl_sa_vn$date[rns[1]], trans_cl_sa_vn$time[rns[1]],"max.diff =",
              max(abs(trans_cl_sa_vn$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  
  rect(xleft = 1351, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1801, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  box(bty="l", lwd=2)

  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_cl_sa_vn$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n")
  }
# dev.off()


# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/FAB_plots_VN_2023-09-19.pdf")
for(i in 1:length(unique(trans_fab_sa_vn$ID))){
  rns <- which(trans_fab_sa_vn$ID==unique(trans_fab_sa_vn$ID)[i])
  
  plot(x=as.numeric(colnames(trans_fab_sa_vn)[6:1535]), y=trans_fab_sa_vn[rns[1],6:1535], ylim=c(0,0.1), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=as.numeric(colnames(trans_fab_sa_vn)[6:1535]), y=trans_fab_sa_vn[rns[j],6:1535], col=color_palette[j], lwd=2)
  }
  
  title(paste(trans_fab_sa_vn$ID[rns[1]], " - vector normalized", "\n(", trans_fab_sa_vn$date[rns[1]], trans_fab_sa_vn$time[rns[1]],"max.diff =",
              max(abs(trans_fab_sa_vn$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  
  rect(xleft = 1351, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1801, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  box(bty="l", lwd=2)
  
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fab_sa_vn$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n")
}
# dev.off()



# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Freiburg_plots_VN_2023-09-19.pdf")
for(i in 1:length(unique(trans_fr_sa_vn$ID))){
  rns <- which(trans_fr_sa_vn$ID==unique(trans_fr_sa_vn$ID)[i])
  
  plot(x=as.numeric(colnames(trans_fr_sa_vn)[6:1535]), y=trans_fr_sa_vn[rns[1],6:1535], ylim=c(0,0.1), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=as.numeric(colnames(trans_fr_sa_vn)[6:1535]), y=trans_fr_sa_vn[rns[j],6:1535], col=color_palette[j], lwd=2)
  }
  
  title(paste(trans_fr_sa_vn$ID[rns[1]], " - vector normalized", "\n(", trans_fr_sa_vn$date[rns[1]], trans_fr_sa_vn$time[rns[1]],"max.diff =",
              max(abs(trans_fr_sa_vn$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  
  rect(xleft = 1351, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1801, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  box(bty="l", lwd=2)
  
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fr_sa_vn$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n")
}
# dev.off()


#_______________________________________________________________________________
# Plotting to check and identify outliers (not vector normalized - to check for differences in irradiance)

# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Cloquet_plots_2023-09-19.pdf")
for(i in 1:length(unique(trans_cl_sa$ID))){
  rns <- which(trans_cl_sa$ID==unique(trans_cl_sa$ID)[i])
  
  plot(x=seq(350,2500,1), y=trans_cl_sa[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=seq(350,2500,1), y=trans_cl_sa[rns[j],12:2162], col=color_palette[j], lwd=2)
  }
  
  title(paste(trans_cl_sa$ID[rns[1]], "\n(", trans_cl_sa$date[rns[1]], trans_cl_sa$time[rns[1]],"max.diff =",
              max(abs(trans_cl_sa$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  
  abline(h=1, col="gray", lty=2)
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_cl_sa$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n")
}
# dev.off()


# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/FAB_plots_2023-09-19.pdf")
for(i in 1:length(unique(trans_fab_sa$ID))){
  rns <- which(trans_fab_sa$ID==unique(trans_fab_sa$ID)[i])
  
  plot(x=seq(350,2500,1), y=trans_fab_sa[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){
    lines(x=seq(350,2500,1), y=trans_fab_sa[rns[j],12:2162], col=color_palette[j], lwd=2)
  }
  
  title(paste(trans_fab_sa$ID[rns[1]], "\n(", trans_fab_sa$date[rns[1]], trans_fab_sa$time[rns[1]],"max.diff =",
              max(abs(trans_fab_sa$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fab_sa$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n")
  # legend("topleft", legend=rns, col=color_palette, lty=1, lwd=2, bty="n")
  abline(h=1, col="gray", lty=2)
}
# dev.off()


# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Freiburg_plots_2023-09-19.pdf")
for(i in 1:length(unique(trans_fr_sa$ID))){
  rns <- which(trans_fr_sa$ID==unique(trans_fr_sa$ID)[i])
  
  plot(x=seq(350,2500,1), y=trans_fr_sa[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  
  for(j in 1:length(rns)){
    lines(x=seq(350,2500,1), y=trans_fr_sa[rns[j],12:2162], col=color_palette[j], lwd=2)
  }
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fr_sa$time_difference[rns]), " s)",sep=""), col=color_palette, lty=1, lwd=2, bty="n", cex=0.5)
  # legend("topleft", legend=rns, col=color_palette, lty=1, lwd=2, bty="n")
  title(paste(trans_fr_sa$ID[rns[1]], "\n(", trans_fr_sa$date[rns[1]], trans_fr_sa$time[rns[1]],"max.diff =",max(abs(trans_fr_sa$time_difference[rns])),")"), 
        adj=0, font.main=1, cex.main=1)
  
  # rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  # rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  # rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  # box(bty="l")
  abline(h=1, col="gray", lty=2)
}
# dev.off()


#_______________________________________________________________________________
# ID outliers and plot with outliers identified:
# Cloquet:
trans_cl_sa$Exclude <- 0
trans_cl_sa$MaybeExclude <- 0
trans_cl_sa$MaybeExclude_set2 <- 0
trans_cl_sa$Exclude_vn <- 0
trans_cl_sa$MaybeExclude_vn <- 0
rns <- c(399, 25, 76, 890, 907, 858, 357, 386, 388, 978, 399, 967, 407, 417, 1051, 740) # definitely exclude
rns2 <- c(16, 21, 159, 222, 249,338, 360,361,362, 363,364,365, 994, 1039, 1001, 811, 854 ) # maybe exclude
rns3 <- c(118, 171, 200, 863,864, 870, 295,296, 909, 958, 989, 1038,1049) # maybe, maybe exclude
rns_vn <- c(16, 17, 20, 21, 25, 26, 27, 76, 222, 248, 249, 250, 907, 357, 386, 384, 388, 387, 407, 417, 818 ) # exclude based on vector normalized data (also considering time difference in ASD match; and if evidence of a lot of noise)
rns_vn2 <- c(41, 43, 58, 61, 72, 83, 90, 118, 120, 123, 124, 125, 126, 159, 158, 171, 172, 175, 200, 
             230, 231, 234, 239, 240, 241, 242, 244, 589, 247, 251, 252, 253, 586, 257, 258, 259, 262,
             261, 260, 265, 566, 864, 870, 295, 296, 304, 888, 890, 892, 314, 318, 319, 320, 958, 332, 
             934, 338, 944, 359, 360, 361, 362, 363, 364, 365, 991, 990, 989, 389, 977, 978, 396, 397, 
             399, 401, 967, 402, 966, 1039, 1046, 1048, 1051, 1081, 459, 995, 478, 479, 1069, 715, 740, 823,
             854, 856) # maybe exclude
trans_cl_sa$Exclude[rns] <- 1
trans_cl_sa$MaybeExclude[rns2] <- 1
trans_cl_sa$MaybeExclude_set2[rns3] <- 1
trans_cl_sa$Exclude_vn[rns_vn] <- 1
trans_cl_sa$MaybeExclude_vn[rns_vn2] <- 1


trans_fr_sa$Exclude <- 0
trans_fr_sa$MaybeExclude <- 0
trans_fr_sa$MaybeExclude_set2 <- 0
trans_fr_sa$Exclude_vn <- 0
trans_fr_sa$MaybeExclude_vn <- 0
rns <- c(146, 425, 434, 431, 437, 482 ) # definitely exclude
rns2 <- c(148:150, 218:220, 224:226, 240, 306, 330, 349, 381, 411, 416, 483) # maybe exclude; 224-226 mislabelled(?)
rns3 <- c(278,279, 340, 374, 384:386, 228:229) # maybe, maybe exclude (plot 175_6 angio_b, 194_BEPE_m)
rns_vn <- c(507, 508, 509, 95, 96, 97, 98, 99, 100, 146, 145, 174, 175, 240, 330, 411, 416, 431, 434, 482) # exclude based on vector normalised data
rns_vn2 <- c(443, 444, 445, 147, 144, 192, 228, 229, 230, 306, 582, 349, 354, 374, 372, 373, 381, 387,
             388, 425, 483) 
trans_fr_sa$Exclude[rns] <- 1
trans_fr_sa$MaybeExclude[rns2] <- 1
trans_fr_sa$MaybeExclude_set2[rns3] <- 1
trans_fr_sa$Exclude_vn[rns_vn] <- 1
trans_fr_sa$MaybeExclude_vn[rns_vn2] <- 1


trans_fab_sa$Exclude <- 0
trans_fab_sa$MaybeExclude <- 0
trans_fab_sa$MaybeExclude_set2 <- 0
trans_fab_sa$Exclude_vn <- 0
trans_fab_sa$MaybeExclude_vn <- 0
rns <- c(66, 343:345) # definitely exclude (what is going on with FAB2_9_BEPA_m, b, FAB2_10_m, b and on -- wavelengths cut out at ~1800 nm 343:416)
rns2 <- c(183, 339, 423:425, 426:428, 429:431, 432:434, 437, 441:443, 444:446, 453:455, 456, 453:465, 586) # maybe exclude
rns3 <- c(95, 132:134, 135:137, 186:187, # no longer excluding 180:182, 162:164, 120:122, 210:212 (less sure about these), 114:116 (keeping all and averaging but could exclude 116?), 87:89, 
          234:236, 240:242, 243:245, 246:248, 249:275, 285:287, 288:290, 533:535, 291:293, 536:538,
          294:296, 517:519, 297:299, 520:522, 300:302, 505:507, 303:305, 508:510, 417:419, 447:449, 
          466:471, 472:474, 475:477, 478:480, 481:483, 484:486, 487:489, 490:492, 493:495, 496:498, 
          499:501, 502:504, 511:513, 514:516, 527:532, 545:553, 563:583, 587:589, 593:595, 596:601) # maybe, maybe exclude (FAB1_62_QUAL_m and other measurements around this time -- are calibrations off?)
rns_vn <- c(66, 186, 187) # Exclude based on vector normalised data
rns_vn2 <- c(80, 85, 102, 116, 192, 200, 255, 282:284, 285:287, 288:290, 291:293, 294:296, 297:299, 300:302, 303:305, 328) # could keep 116 but clearly different to other readings on plot; could exclude 89 

trans_fab_sa$Exclude[rns] <- 1
trans_fab_sa$MaybeExclude[rns2] <- 1
trans_fab_sa$MaybeExclude_set2[rns3] <- 1
trans_fab_sa$Exclude_vn[rns_vn] <- 1
trans_fab_sa$MaybeExclude_vn[rns_vn2] <- 1

str_split(trans_cl_sa$ID, "_", simplify=T)
ids_cl <- data.frame(str_split(trans_cl_sa$ID, "_", simplify=T))
colnames(ids_cl) <- c("Block", "Plot", "Mixture", "Height")
trans_cl_sa2 <- cbind(trans_cl_sa, ids_cl)
trans_cl_sa2$Mixture_height <- c(paste(trans_cl_sa2$Mixture, trans_cl_sa2$Height, sep="_"))

##______________________________________________________________________________
## Saving files with outliers identified

write.csv(trans_cl_sa, "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-09-19_version/IDENT-Cloquet_transmittance_svc-asd_outliers-identified.csv")
write.csv(trans_fr_sa, "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-09-19_version/IDENT-Freiburg_transmittance_svc-asd_outliers-identified.csv")
write.csv(trans_fab_sa, "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance/2023-09-19_version/FAB_transmittance_svc-asd_outliers-identified.csv")

## Note: make sure that Freiburg plot 47 and 226 have correct composition assigned (ACPL and LADE-LALA, respectively)

##______________________________________________________________________________
## Plotting to show outliers:
## Cloquet:
# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Cloquet_plots-outliers_2023-09-19.pdf")
for(i in 1:length(unique(trans_cl_sa2$ID))){
  rns <- which(trans_cl_sa2$ID==unique(trans_cl_sa2$ID)[i])
  plot(x=seq(350,2500,1), y=trans_cl_sa2[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
    color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=seq(350,2500,1), y=trans_cl_sa2[rns[j],12:2162], col=color_palette[j], lwd=2, lty=(trans_cl_sa2[rns[j]]$Exclude+1))
  }
  for(k in 1:length(rns)){ # potential outliers highlighted with gray dashes
    if(trans_cl_sa2[rns[k]]$MaybeExclude==1){
      lines(x=seq(350,2500,1), y=trans_cl_sa2[rns[k],12:2162], col="lightgray", lwd=2, lty="5B")}
  }
    for(k in 1:length(rns)){ # potential outliers highlighted with black dots
      if(trans_cl_sa2[rns[k]]$MaybeExclude_set2==1){
        lines(x=seq(350,2500,1), y=trans_cl_sa2[rns[k],12:2162], col="black", lwd=2, lty=3)}
    }
    
  title(paste(trans_cl_sa2$ID[rns[1]], "\n(", trans_cl_sa2$date[rns[1]], trans_cl_sa2$time[rns[1]],"max.diff =",
              max(abs(trans_cl_sa2$time_difference[rns])),")"), adj=0, font.main=1, cex.main=1)
  abline(h=1, col="gray", lty=2)
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_cl_sa2$time_difference[rns]), "s, ", 
                                 trans_cl_sa2$date[rns], ", ", trans_cl_sa2$time[rns], ")", 
                                 sep=""), col=color_palette, lty=1, lwd=2, bty="n", cex=0.7)
  legend("bottomleft", legend=c("White dash = exclude; gray dash = maybe exclude; black dash = maybe exclude set 2"), cex=0.6, bty="n")
}
# dev.off()




## FAB:
# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/FAB_plots-outliers_2023-09-19.pdf")
for(i in 1:length(unique(trans_fab_sa$ID))){
  rns <- which(trans_fab_sa$ID==unique(trans_fab_sa$ID)[i])
  plot(x=seq(350,2500,1), y=trans_fab_sa[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=seq(350,2500,1), y=trans_fab_sa[rns[j],12:2162], col=color_palette[j], lwd=2, lty=(trans_fab_sa[rns[j]]$Exclude+1))
  }
  for(k in 1:length(rns)){ # potential outliers highlighted with gray dashes
    if(trans_fab_sa[rns[k]]$MaybeExclude==1){
      lines(x=seq(350,2500,1), y=trans_fab_sa[rns[k],12:2162], col="lightgray", lwd=2, lty="5B")}
  }
  for(k in 1:length(rns)){ # potential outliers highlighted with black dots
    if(trans_fab_sa[rns[k]]$MaybeExclude_set2==1){
      lines(x=seq(350,2500,1), y=trans_fab_sa[rns[k],12:2162], col="black", lwd=2, lty=3)}
  }
  
  title(paste(trans_fab_sa$ID[rns[1]], "\n(", trans_fab_sa$date[rns[1]], trans_fab_sa$time[rns[1]],"max.diff =",
              max(abs(trans_fab_sa$time_difference[rns])),")"), adj=0, font.main=1, cex.main=1)
  abline(h=1, col="gray", lty=2)
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fab_sa$time_difference[rns]), "s, ", 
                                 trans_fab_sa$date[rns], ", ", trans_fab_sa$time[rns], ")", 
                                 sep=""), col=color_palette, lty=1, lwd=2, bty="n", cex=0.7)
  legend("bottomleft", legend=c("White dash = exclude; gray dash = maybe exclude; black dash = maybe exclude set 2"), cex=0.6, bty="n")
}
# dev.off()



## Freiburg:
# pdf(file="/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Freiburg_plots-outliers_2023-09-19.pdf")
for(i in 1:length(unique(trans_fr_sa$ID))){
  rns <- which(trans_fr_sa$ID==unique(trans_fr_sa$ID)[i])
  plot(x=seq(350,2500,1), y=trans_fr_sa[rns[1],12:2162], ylim=c(0,1.5), bty="l", type="n", 
       ylab="Transmittance", xlab="Wavelength (nm)")
  color_palette <- colorRampPalette(c("#0D47A1", "#00BCD4", "#FBC02D"))(length(rns))
  for(j in 1:length(rns)){ # outliers highlighted with white dashes
    lines(x=seq(350,2500,1), y=trans_fr_sa[rns[j],12:2162], col=color_palette[j], lwd=2, lty=(trans_fr_sa[rns[j]]$Exclude+1))
  }
  for(k in 1:length(rns)){ # potential outliers highlighted with gray dashes
    if(trans_fr_sa[rns[k]]$MaybeExclude==1){
      lines(x=seq(350,2500,1), y=trans_fr_sa[rns[k],12:2162], col="lightgray", lwd=2, lty="5B")}
  }
  for(k in 1:length(rns)){ # potential outliers highlighted with black dots
    if(trans_fr_sa[rns[k]]$MaybeExclude_set2==1){
      lines(x=seq(350,2500,1), y=trans_fr_sa[rns[k],12:2162], col="black", lwd=2, lty=3)}
  }
  
  title(paste(trans_fr_sa$ID[rns[1]], "\n(", trans_fr_sa$date[rns[1]], trans_fr_sa$time[rns[1]],"max.diff =",
              max(abs(trans_fr_sa$time_difference[rns])),")"), adj=0, font.main=1, cex.main=1)
  abline(h=1, col="gray", lty=2)
  legend("topleft", legend=paste(rns, " (time diff: ",abs(trans_fr_sa$time_difference[rns]), "s, ", 
                                 trans_fr_sa$date[rns], ", ", trans_fr_sa$time[rns], ")", 
                                 sep=""), col=color_palette, lty=1, lwd=2, bty="n", cex=0.7)
  legend("bottomleft", legend=c("White dash = exclude; gray dash = maybe exclude; black dash = maybe exclude set 2"), cex=0.6, bty="n")
}
# dev.off()

