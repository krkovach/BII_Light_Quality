##______________________________________________________________________________
## Loading data:
## These are canopy transmittance (sa = SVC vs ASD) with outliers identified
trans_fab_sa <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/Transmittance/2023-09-19_version/FAB_transmittance_svc-asd_outliers-identified.csv")
trans_cl_sa <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/Transmittance/2023-09-19_version/IDENT-Cloquet_transmittance_svc-asd_outliers-identified.csv")
trans_fr_sa <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/Transmittance/2023-09-19_version/IDENT-Freiburg_transmittance_svc-asd_outliers-identified.csv")


## These are LAI values:
lai_cl <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/LAI/IDENTCloquet2021_LAI2000_withCorrection.csv")
lai_fab <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/LAI/FAB1_data_LAI2000_withCorrection.csv") # this excludes all FAB2 data
lai_fr <- read.csv("/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing/LAI/IDENTFreiburg2021_LAI2000_withCorrection.csv")

##______________________________________________________________________________
## Adding mixture, block, plot info:
cl_ids <- data.frame(str_split(trans_cl_sa$ID, "_", simplify=T))
colnames(cl_ids) <- c("Block", "Plot", "Mixture", "CanopyHt")
trans_cl_sa2 <- cbind(trans_cl_sa, cl_ids)
trans_cl_sa2$BlMix <- paste(trans_cl_sa2$Block, trans_cl_sa2$Mixture, sep="_")

fab_ids <- data.frame(str_split(trans_fab_sa$ID, "_", simplify=T))
colnames(fab_ids) <- c("Site", "Plot", "Mixture", "CanopyHt")
trans_fab_sa2 <- cbind(trans_fab_sa, fab_ids)
trans_fab_sa2$PlMix <- paste(trans_fab_sa2$Plot, trans_fab_sa2$Mixture, sep="_")

fr_ids <- data.frame(str_split(trans_fr_sa$ID, "_", simplify=T))
colnames(fr_ids) <- c("Plot", "Mixture", "CanopyHt")
trans_fr_sa2 <- cbind(trans_fr_sa, fr_ids)
trans_fr_sa2$PlMix <- paste(trans_fr_sa2$Plot, trans_fr_sa2$Mixture, sep="_")
trans_fr_sa
# Correcting mislabelled mixtures
trans_fr_sa2[which(trans_fr_sa2$Plot==47),]$Mixture <- "ACPL" # ACPL-PIAB plots should be plot #45
trans_fr_sa2[which(trans_fr_sa2$Plot==47),]$PlMix <- "47_ACPL"
trans_fr_sa2[which(trans_fr_sa2$Plot==47&trans_fr_sa2$CanopyHt=="m"),]$ID <- "47_ACPL_m"
trans_fr_sa2[which(trans_fr_sa2$Plot==47&trans_fr_sa2$CanopyHt=="b"),]$ID <- "47_ACPL_b"

trans_fr_sa2[which(trans_fr_sa2$Plot==226),]$Mixture <- "LADE-LALA" # Some measurements labelled LALA-LADE rather than LADE-LALA
trans_fr_sa2[which(trans_fr_sa2$Plot==226),]$PlMix <- "226_LADE-LALA"
trans_fr_sa2[which(trans_fr_sa2$Plot==226&trans_fr_sa2$CanopyHt=="m"),]$ID <- "226_LADE-LALA_m"
trans_fr_sa2[which(trans_fr_sa2$Plot==226&trans_fr_sa2$CanopyHt=="b"),]$ID <- "226_LADE-LALA_b"


##______________________________________________________________________________
## Removing outliers
trans_fab_sa2$Exclude_all <- rowSums(trans_fab_sa2[,which(colnames(trans_fab_sa2)%in%c("Exclude", "MaybeExclude", "MaybeExclude_set2", "Exclude_vn", "MaybeExclude_vn"))])
trans_fab_sa3 <- trans_fab_sa2[which(trans_fab_sa2$Exclude_all==0),]

trans_cl_sa2$Exclude_all <- rowSums(trans_cl_sa2[,which(colnames(trans_cl_sa2)%in%c("Exclude", "MaybeExclude", "MaybeExclude_set2", "Exclude_vn", "MaybeExclude_vn"))])
trans_cl_sa3 <- trans_cl_sa2[which(trans_cl_sa2$Exclude_all==0),]

trans_fr_sa2$Exclude_all <- rowSums(trans_fr_sa2[,which(colnames(trans_fr_sa2)%in%c("Exclude", "MaybeExclude", "MaybeExclude_set2", "Exclude_vn", "MaybeExclude_vn"))])
trans_fr_sa3 <- trans_fr_sa2[which(trans_fr_sa2$Exclude_all==0),]


##______________________________________________________________________________
## Average spectra for plots
trans_fab_ag <- aggregate(trans_fab_sa3[,c(13:2163)], by=list(ID=trans_fab_sa3$ID, 
                                                              Site=trans_fab_sa3$Site, 
                                                              Plot=trans_fab_sa3$Plot, 
                                                              Mix=trans_fab_sa3$Mixture, 
                                                              CanopyHt=trans_fab_sa3$CanopyHt, 
                                                              PlMix=trans_fab_sa3$PlMix), FUN="mean")
trans_cl_ag <- aggregate(trans_cl_sa3[,c(13:2163)], by=list(ID=trans_cl_sa3$ID, 
                                                            Block=trans_cl_sa3$Block, 
                                                            Plot=trans_cl_sa3$Plot, 
                                                            Mix=trans_cl_sa3$Mixture, 
                                                            CanopyHt=trans_cl_sa3$CanopyHt, 
                                                            BlMix=trans_cl_sa3$BlMix), FUN="mean")
trans_fr_ag <- aggregate(trans_fr_sa3[,c(13:2163)], by=list(ID=trans_fr_sa3$ID,
                                                            Plot=trans_fr_sa3$Plot, 
                                                            Mix=trans_fr_sa3$Mixture, 
                                                            CanopyHt=trans_fr_sa3$CanopyHt, 
                                                            PlMix=trans_fr_sa3$PlMix), FUN="mean")

trans_fab1_ag <- trans_fab_ag[which(substring(trans_fab_ag$ID,1,4)=="FAB1"),] # separating FAB1 from FAB2
trans_fab2_ag <- trans_fab_ag[which(substring(trans_fab_ag$ID,1,4)=="FAB2"),]

trans_fab1_ag$ID2 <- trans_fab1_ag$ID # adding standardised form IDs for FAB1 and FAB2
trans_fab2_ag$ID2 <- trans_fab2_ag$ID
trans_fab1_ag$ID <- with(trans_fab1_ag, paste(Site,Plot,CanopyHt, sep="_"))
trans_fab2_ag$ID <- with(trans_fab2_ag, paste(Site,Plot,CanopyHt, sep="_")) 

# Fr: 175_6 angio_b is missing

##______________________________________________________________________________
## LAI: checking, averaging, adding to spectral data frames
# whether to clean?
#...Cloquet:
lai_cl$ID <- with(lai_cl, paste(BlPl,Mixture,Height, sep="_"))
lai_cl$LAI <- as.numeric(lai_cl$LAI)
lai_cl$MeanGapFraction <- as.numeric(lai_cl$MeanGapFraction)

lai_cl_av_all <- aggregate(cbind(LAI_allloc=LAI, LAI_ExclOuter_allloc=LAI_ExclOutermostRing)~ID, lai_cl, FUN=mean) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
lai_cl_av_specloc <- aggregate(cbind(LAI_specloc=LAI, LAI_ExclOuter_specloc=LAI_ExclOutermostRing)~ID, lai_cl[lai_cl$Location=="NE",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured
gf_cl_av_all <- aggregate(MeanGapFraction~ID, lai_cl, FUN=mean, na.rm=T) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
gf_cl_av_specloc <- aggregate(cbind(MeanGapFraction_specloc=MeanGapFraction)~ID, lai_cl[lai_cl$Location=="NE",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured
lai_cl_av_all_corrected <- aggregate(cbind(LAI_allloc_corrected=LAI_corrected, LAI_ExclOuter_allloc_corrected=LAI_ExclOuter_corrected)~ID, lai_cl, FUN=mean) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
lai_cl_av_specloc_corrected <- aggregate(cbind(LAI_specloc_corrected=LAI_corrected, LAI_ExclOuter_specloc_corrected=LAI_ExclOuter_corrected)~ID, lai_cl[lai_cl$Location=="NE",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured

lai_cl_av1 <- merge(lai_cl_av_all, lai_cl_av_specloc)
lai_cl_av2 <- merge(lai_cl_av1, gf_cl_av_all)
lai_cl_av3 <- merge(lai_cl_av2, gf_cl_av_specloc)
lai_cl_av4 <- merge(lai_cl_av3, lai_cl_av_all_corrected)
lai_cl_av <- merge(lai_cl_av4, lai_cl_av_specloc_corrected)

#...Freiburg:
lai_fr2 <- merge(lai_fr, unique(trans_fr_ag[,2:3]), by.x="Plot_number", by.y="Plot", all.x=T, all.y=F) # adding Mix
lai_fr2$ID <- with(lai_fr2, paste(Plot_number, Mix, Height, sep="_"))
# (this was incorrect -- repeating rows because of mislabelled mixtures, now fixed)

lai_fr_av_all <- aggregate(cbind(LAI_allloc=LAI, LAI_ExclOuter_allloc=LAI_ExclOutermostRing)~ID, lai_fr2, FUN=mean, na.rm=T) # aggregating all
lai_fr_av_specloc <- aggregate(cbind(LAI_specloc=LAI, LAI_ExclOuter_specloc=LAI_ExclOutermostRing)~ID, lai_fr2[lai_fr2$Location=="NE",], FUN=mean) # selecting just location where spectra measured
gf_fr_av_all <- aggregate(MeanGapFraction~ID, lai_fr2, FUN=mean, na.rm=T) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
gf_fr_av_specloc <- aggregate(cbind(MeanGapFraction_specloc=MeanGapFraction)~ID, lai_fr2[lai_fr2$Location=="NE",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured
lai_fr_av_all_corrected <- aggregate(cbind(LAI_allloc_corrected=LAI_corrected, LAI_ExclOuter_allloc_corrected=LAI_ExclOuter_corrected)~ID, lai_fr2, FUN=mean) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
lai_fr_av_specloc_corrected <- aggregate(cbind(LAI_specloc_corrected=LAI_corrected, LAI_ExclOuter_specloc_corrected=LAI_ExclOuter_corrected)~ID, lai_fr2[lai_fr2$Location=="NE",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured

lai_fr_av1 <- merge(lai_fr_av_all, lai_fr_av_specloc, all.x=T)
lai_fr_av2 <- merge(lai_fr_av1, gf_fr_av_all, all.x=T)
lai_fr_av3 <- merge(lai_fr_av2, gf_fr_av_specloc, all.x=T)
lai_fr_av4 <- merge(lai_fr_av3, lai_fr_av_all_corrected, all.x=T)
lai_fr_av <- merge(lai_fr_av4, lai_fr_av_specloc_corrected, all.x=T)


### NB: NO LAI value for spectral location (NE) for Freiburg 2_PIGL_m --- used average of the other locations to infill (i.e., same LAI as all locations)
lai_fr_av[rowSums(is.na(lai_fr_av)) > 0,]
lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"), which(colnames(lai_fr_av)=="LAI_specloc")] <- 
  lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"),which(colnames(lai_fr_av)=="LAI_allloc")]
lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"), which(colnames(lai_fr_av)=="LAI_ExclOuter_specloc")] <- 
  lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"),which(colnames(lai_fr_av)=="LAI_ExclOuter_allloc")]
lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"), which(colnames(lai_fr_av)=="MeanGapFraction_specloc")] <- 
  lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"),which(colnames(lai_fr_av)=="MeanGapFraction")]
lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"), which(colnames(lai_fr_av)=="LAI_specloc_corrected")] <- 
  lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"),which(colnames(lai_fr_av)=="LAI_allloc_corrected")]
lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"), which(colnames(lai_fr_av)=="LAI_ExclOuter_specloc_corrected")] <- 
  lai_fr_av[which(lai_fr_av$ID=="2_PIGL_m"),which(colnames(lai_fr_av)=="LAI_ExclOuter_allloc_corrected")]


#...FAB:
lai_fab$ID <- with(lai_fab, paste(Expt,Plot,Height, sep="_"))

lai_fab_av_all <- aggregate(cbind(LAI_allloc=LAI, LAI_ExclOuter_allloc=LAI_ExclOuter)~ID, lai_fab, FUN=mean, na.rm=T) # aggregating all
lai_fab_av_specloc <- aggregate(cbind(LAI_specloc=LAI, LAI_ExclOuter_specloc=LAI_ExclOuter)~ID, lai_fab[lai_fab$Location=="C",], FUN=mean, na.rm=T) # selecting just location where spectra measured
gf_fab_av_all <- aggregate(MeanGapFraction~ID, lai_fab, FUN=mean, na.rm=T) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
gf_fab_av_specloc <- aggregate(cbind(MeanGapFraction_specloc=MeanGapFraction)~ID, lai_fab[lai_fab$Location=="C",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured
lai_fab_av_all_corrected <- aggregate(cbind(LAI_allloc_corrected=LAI_corrected, LAI_ExclOuter_allloc_corrected=LAI_ExclOuter_corrected)~ID, lai_fab, FUN=mean) # aggregating all measured locations (at the exact same spot as the spectra as well as three other locations)
lai_fab_av_specloc_corrected <- aggregate(cbind(LAI_specloc_corrected=LAI_corrected, LAI_ExclOuter_specloc_corrected=LAI_ExclOuter_corrected)~ID, lai_fab[lai_fab$Location=="C",], FUN=mean, na.rm=T) # selecting just the location on the plot where spectra were measured

lai_fab_av1 <- merge(lai_fab_av_all, lai_fab_av_specloc)
lai_fab_av2 <- merge(lai_fab_av1, gf_fab_av_all)
lai_fab_av3 <- merge(lai_fab_av2, gf_fab_av_specloc)
lai_fab_av4 <- merge(lai_fab_av3, lai_fab_av_all_corrected)
lai_fab_av <- merge(lai_fab_av4, lai_fab_av_specloc_corrected)


##______________________________________________________________________________
## Creating dataframes to save:
# removing noisy parts of spectrum
trans_cl_ag_outrm <- trans_cl_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)]
trans_fr_ag_outrm <- trans_fr_ag[,c(1:5,6:1006, 1106:1456, 1636:1856)]
trans_fab1_ag_outrm <- trans_fab1_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)]
trans_fab2_ag_outrm <- trans_fab2_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)]

trans_cl_ag_lai <- merge(trans_cl_ag_outrm, lai_cl_av, by="ID", all.x=T, all.y=F)
trans_fr_ag_lai <- merge(trans_fr_ag_outrm, lai_fr_av, by="ID", all.x=T, all.y=F)

trans_fab1_ag_lai <- merge(trans_fab1_ag_outrm, lai_fab_av, by="ID", all.x=T, all.y=F)
trans_fab2_ag_lai <- merge(trans_fab2_ag_outrm, lai_fab_av, by="ID", all.x=T, all.y=F)

trans_cl_ag_lai$Site <- "Cloquet"
trans_fr_ag_lai$Site <- "Freiburg"

# trans_cl_ag_lai$Plot <- with(trans_cl_ag_lai, paste(Block, Plot, sep="_")) # only run this once

trans_cl_ag_lai$ID
trans_fr_ag_lai$ID
trans_fab1_ag_lai$ID
trans_fab2_ag_lai$ID

trans_cl_ag_lai2 <- cbind(trans_cl_ag_lai[1:228,c(1, 1590, 3:6, 1580:1589, 7:1579)]) # excluding rows that contain reference readings (dark, open, spectralon)
trans_fr_ag_lai2 <- cbind(trans_fr_ag_lai[2:112,c(1, 1589, 2:5, 1579:1588, 6:1578)]) # previously had 2:116 (175 6 angio no B)
trans_fab1_ag_lai2 <- cbind(trans_fab1_ag_lai[1:68,c(1, 2:6, 1580:1589, 7:1579)])
trans_fab2_ag_lai2 <- cbind(trans_fab2_ag_lai[1:36,c(1, 2:6, 1580:1589, 7:1579)]) # not adding LAI here

# updating column names to be consistent:
colnames(trans_cl_ag_lai2)[1:10]
colnames(trans_fr_ag_lai2)[6] <- "BlMix"
colnames(trans_fab1_ag_lai2)[6]  <- "BlMix"
colnames(trans_fab2_ag_lai2)[6]  <- "BlMix"

# checking columns consistently named:
trans_cl_ag_lai2[1:3,1:12]
trans_fr_ag_lai2[1:3,1:12]
trans_fab1_ag_lai2[1:3,1:12]
trans_fab2_ag_lai2[1:3,1:12]

trans_all_lai <- rbind(trans_cl_ag_lai2, trans_fr_ag_lai2, trans_fab1_ag_lai2, trans_fab2_ag_lai2)

file_dir <- "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance"
fwrite(trans_all_lai, file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230919.txt"), sep = "\t")  

##______________________________________________________________________________
## Plotting mean spectra alongside individual spectra (flagging included and excluded):

## Cloquet:
# remove noisy (outlier) sections of the spectrum
# keep the following: 350-1350, 1450-1800, 1980-2200
# exclude the following: 1351-1449, 1801-1979, 2201-2500
trans_cl_ag_outrm <- trans_cl_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)] # average spectra per plot and canopy ht (without bad readings)
trans_cl_sa2_outrm <- trans_cl_sa2[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)] # individual spectra (including those identified as bad readings)
trans_cl_sa3_outrm <- trans_cl_sa3[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)] # individual spectra excluding those identified as bad readings
  


# pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Cloquet_plots_mean-transmittance.pdf")
for(i in 1:length(unique(trans_cl_ag_outrm$BlMix))){
  BlMix <- unique(trans_cl_ag_outrm$BlMix)[i]
  a_m <- trans_cl_ag_outrm[which(trans_cl_ag_outrm$BlMix==BlMix&trans_cl_ag_outrm$CanopyHt=="m"),]
  a_b <- trans_cl_ag_outrm[which(trans_cl_ag_outrm$BlMix==BlMix&trans_cl_ag_outrm$CanopyHt=="b"),]
  
  aa_m <- trans_cl_sa2_outrm[which(trans_cl_sa2_outrm$BlMix==BlMix&trans_cl_sa2_outrm$CanopyHt=="m"),]
  aa_b <- trans_cl_sa2_outrm[which(trans_cl_sa2_outrm$BlMix==BlMix&trans_cl_sa2_outrm$CanopyHt=="b"),]
  
  aaa_m <- trans_cl_sa3_outrm[which(trans_cl_sa3_outrm$BlMix==BlMix&trans_cl_sa3_outrm$CanopyHt=="m"),]
  aaa_b <- trans_cl_sa3_outrm[which(trans_cl_sa3_outrm$BlMix==BlMix&trans_cl_sa3_outrm$CanopyHt=="b"),]

      # mean spectra:
  x <- c(350:1350,1450:1800,1980:2200)
  plot(x=x, y=a_m[i,7:ncol(a_m)], ylim=c(0,1.3), type="n", bty="l", 
     xlab="wavelength (nm)", ylab="Transmittance")

    # adding lines for mean spectrum:
  lines(x=x, y=a_m[1,7:ncol(a_m)], lwd=3, col="#FBC02D")
  lines(x=x, y=a_b[1,7:ncol(a_m)], lwd=3, col="#0D47A1") # "#FBC02D", "##00BCD4"
  title(BlMix)
  
  # adding lines for spectral outliers
  for(i in 1:nrow(aa_m)){lines(x=x, y=aa_m[i,13:(ncol(aa_m)-2)], lwd=1, lty=2, col="#FF990050")} #"#FBC02D99"
  for(i in 1:nrow(aa_b)){lines(x=x, y=aa_b[i,13:(ncol(aa_b)-2)], lwd=1, lty=2, col="#3399FF50")} #"#0D47A199"

  # adding lines for individual spectra:
  for(i in 1:nrow(aaa_m)){lines(x=x, y=aaa_m[i,13:(ncol(aaa_m)-2)], lwd=1, lty=1,  col="#FF9900")}
  for(i in 1:nrow(aaa_b)){lines(x=x, y=aaa_b[i,13:(ncol(aaa_b)-2)], lwd=1, lty=1, col="#3399FF")}
  
  # blanking out noisy sections of the spectrum (straight lines connecting points)
  rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  
  box(bty="l")
  legend("topleft", legend=c("middle canopy - mean", "included spectrum", "excluded spectrum", "bottom canopy - mean", "included spectrum", "excluded spectrum"), 
         col=c("#FBC02D","#FF9900", "#FF990050", "#0D47A1","#3399FF", "#3399FF50"), lwd=c(3,1,1, 3, 1, 1), lty=c(1,1,2,1,1,2), bty="n", cex=0.7)
}
# dev.off()


## Freiburg:
# remove noisy (outlier) sections of the spectrum
# keep the following: 350-1350, 1450-1800, 1980-2200
# exclude the following: 1351-1449, 1801-1979, 2201-2500
trans_fr_ag_outrm <- trans_fr_ag[,c(1:5,6:1006, 1106:1456, 1636:1856)]
trans_fr_sa2_outrm <- trans_fr_sa2[,c(1:12,13:1013, 1113:1463, 1643:1863, 2171:2172)]
trans_fr_sa3_outrm <- trans_fr_sa3[,c(1:12,13:1013, 1113:1463, 1643:1863, 2171:2172)]


# pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/Freiburg_plots_mean-transmittance.pdf")
for(i in 1:length(unique(trans_fr_ag_outrm$PlMix))){
  PlMix <- unique(trans_fr_ag_outrm$PlMix)[i]
  a_m <- trans_fr_ag_outrm[which(trans_fr_ag_outrm$PlMix==PlMix&trans_fr_ag_outrm$CanopyHt=="m"),]
  a_b <- trans_fr_ag_outrm[which(trans_fr_ag_outrm$PlMix==PlMix&trans_fr_ag_outrm$CanopyHt=="b"),]
  
  aa_m <- trans_fr_sa2_outrm[which(trans_fr_sa2_outrm$PlMix==PlMix&trans_fr_sa2_outrm$CanopyHt=="m"),]
  aa_b <- trans_fr_sa2_outrm[which(trans_fr_sa2_outrm$PlMix==PlMix&trans_fr_sa2_outrm$CanopyHt=="b"),]
  
  aaa_m <- trans_fr_sa3_outrm[which(trans_fr_sa3_outrm$PlMix==PlMix&trans_fr_sa3_outrm$CanopyHt=="m"),]
  aaa_b <- trans_fr_sa3_outrm[which(trans_fr_sa3_outrm$PlMix==PlMix&trans_fr_sa3_outrm$CanopyHt=="b"),]
  
  # mean spectra:
  x <- c(350:1350,1450:1800,1980:2200)
  plot(x=x, y=a_m[i,6:ncol(a_m)], ylim=c(0,1.3), type="n", bty="l", 
       xlab="wavelength (nm)", ylab="Transmittance")
  
  # adding lines for mean spectrum:
  lines(x=x, y=a_m[1,6:ncol(a_m)], lwd=3, col="#FBC02D")
  lines(x=x, y=a_b[1,6:ncol(a_m)], lwd=3, col="#0D47A1") # "#FBC02D", "##00BCD4"
  title(PlMix)
  
  # adding lines for spectral outliers
  for(i in 1:nrow(aa_m)){lines(x=x, y=aa_m[i,13:(ncol(aa_m)-2)], lwd=1, lty=2, col="#FF990050")} #"#FBC02D99"
  for(i in 1:nrow(aa_b)){lines(x=x, y=aa_b[i,13:(ncol(aa_b)-2)], lwd=1, lty=2, col="#3399FF50")} #"#0D47A199"
  
  # adding lines for individual spectra:
  for(i in 1:nrow(aaa_m)){lines(x=x, y=aaa_m[i,13:(ncol(aaa_m)-2)], lwd=1, lty=1,  col="#FF9900")}
  for(i in 1:nrow(aaa_b)){lines(x=x, y=aaa_b[i,13:(ncol(aaa_b)-2)], lwd=1, lty=1, col="#3399FF")}
  
  # blanking out noisy sections of the spectrum (straight lines connecting points)
  rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  
  box(bty="l")
  legend("topleft", legend=c("middle canopy - mean", "included spectrum", "excluded spectrum", "bottom canopy - mean", "included spectrum", "excluded spectrum"), 
         col=c("#FBC02D","#FF9900", "#FF990050", "#0D47A1","#3399FF", "#3399FF50"), lwd=c(3,1,1, 3, 1, 1), lty=c(1,1,2,1,1,2), bty="n", cex=0.7)
}
# dev.off()


## FAB:
# remove noisy (outlier) sections of the spectrum
# keep the following: 350-1350, 1450-1800, 1980-2200
# exclude the following: 1351-1449, 1801-1979, 2201-2500

trans_fab1_ag_outrm <- trans_fab1_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)]
trans_fab2_ag_outrm <- trans_fab2_ag[,c(1:6,7:1007, 1107:1457, 1637:1857)]

trans_fab1_sa2 <- trans_fab_sa2[which(substring(trans_fab_sa2$ID,1,4)=="FAB1"),]
trans_fab2_sa2 <- trans_fab_sa2[which(substring(trans_fab_sa2$ID,1,4)=="FAB2"),]

trans_fab1_sa3 <- trans_fab_sa3[which(substring(trans_fab_sa3$ID,1,4)=="FAB1"),]
trans_fab2_sa3 <- trans_fab_sa3[which(substring(trans_fab_sa3$ID,1,4)=="FAB2"),]

trans_fab1_sa2_outrm <- trans_fab1_sa2[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)]
trans_fab2_sa2_outrm <- trans_fab2_sa2[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)]

trans_fab1_sa3_outrm <- trans_fab1_sa3[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)]
trans_fab2_sa3_outrm <- trans_fab2_sa3[,c(1:12,13:1013, 1113:1463, 1643:1863, 2172:2173)]

   ##...FAB1
# pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/FAB1_plots_mean-transmittance.pdf")
for(i in 1:length(unique(trans_fab1_ag_outrm$PlMix))){
  PlMix <- unique(trans_fab1_ag_outrm$PlMix)[i]
  a_m <- trans_fab1_ag_outrm[which(trans_fab1_ag_outrm$PlMix==PlMix&trans_fab1_ag_outrm$CanopyHt=="m"),]
  a_b <- trans_fab1_ag_outrm[which(trans_fab1_ag_outrm$PlMix==PlMix&trans_fab1_ag_outrm$CanopyHt=="b"),]
  
  aa_m <- trans_fab1_sa2_outrm[which(trans_fab1_sa2_outrm$PlMix==PlMix&trans_fab1_sa2_outrm$CanopyHt=="m"),]
  aa_b <- trans_fab1_sa2_outrm[which(trans_fab1_sa2_outrm$PlMix==PlMix&trans_fab1_sa2_outrm$CanopyHt=="b"),]
  
  aaa_m <- trans_fab1_sa3_outrm[which(trans_fab1_sa3_outrm$PlMix==PlMix&trans_fab1_sa3_outrm$CanopyHt=="m"),]
  aaa_b <- trans_fab1_sa3_outrm[which(trans_fab1_sa3_outrm$PlMix==PlMix&trans_fab1_sa3_outrm$CanopyHt=="b"),]
  
  # mean spectra:
  x <- c(350:1350,1450:1800,1980:2200)
  plot(x=x, y=a_m[i,7:ncol(a_m)], ylim=c(0,1.3), type="n", bty="l", 
       xlab="wavelength (nm)", ylab="Transmittance")
  
  # adding lines for mean spectrum:
  lines(x=x, y=a_m[1,7:ncol(a_m)], lwd=3, col="#FBC02D")
  lines(x=x, y=a_b[1,7:ncol(a_m)], lwd=3, col="#0D47A1") # "#FBC02D", "##00BCD4"
  title(PlMix)
  
  # adding lines for spectral outliers
  for(i in 1:nrow(aa_m)){lines(x=x, y=aa_m[i,13:(ncol(aa_m)-2)], lwd=1, lty=2, col="#FF990050")} #"#FBC02D99"
  for(i in 1:nrow(aa_b)){lines(x=x, y=aa_b[i,13:(ncol(aa_b)-2)], lwd=1, lty=2, col="#3399FF50")} #"#0D47A199"
  
  # adding lines for individual spectra:
  for(i in 1:nrow(aaa_m)){lines(x=x, y=aaa_m[i,13:(ncol(aaa_m)-2)], lwd=1, lty=1,  col="#FF9900")}
  for(i in 1:nrow(aaa_b)){lines(x=x, y=aaa_b[i,13:(ncol(aaa_b)-2)], lwd=1, lty=1, col="#3399FF")}
  
  # blanking out noisy sections of the spectrum (straight lines connecting points)
  rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  
  box(bty="l")
  legend("topleft", legend=c("middle canopy - mean", "included spectrum", "excluded spectrum", "bottom canopy - mean", "included spectrum", "excluded spectrum"), 
         col=c("#FBC02D","#FF9900", "#FF990050", "#0D47A1","#3399FF", "#3399FF50"), lwd=c(3,1,1, 3, 1, 1), lty=c(1,1,2,1,1,2), bty="n", cex=0.7)
}
# dev.off()


   ##...FAB2
# pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/visual_assessment/2023-09-19_version/FAB2_plots_mean-transmittance.pdf")
for(i in 1:length(unique(trans_fab2_ag_outrm$PlMix))){
  PlMix <- unique(trans_fab2_ag_outrm$PlMix)[i]
  a_m <- trans_fab2_ag_outrm[which(trans_fab2_ag_outrm$PlMix==PlMix&trans_fab2_ag_outrm$CanopyHt=="m"),]
  a_b <- trans_fab2_ag_outrm[which(trans_fab2_ag_outrm$PlMix==PlMix&trans_fab2_ag_outrm$CanopyHt=="b"),]

  aa_m <- trans_fab2_sa2_outrm[which(trans_fab2_sa2_outrm$PlMix==PlMix&trans_fab2_sa2_outrm$CanopyHt=="m"),]
  aa_b <- trans_fab2_sa2_outrm[which(trans_fab2_sa2_outrm$PlMix==PlMix&trans_fab2_sa2_outrm$CanopyHt=="b"),]

  aaa_m <- trans_fab2_sa3_outrm[which(trans_fab2_sa3_outrm$PlMix==PlMix&trans_fab2_sa3_outrm$CanopyHt=="m"),]
  aaa_b <- trans_fab2_sa3_outrm[which(trans_fab2_sa3_outrm$PlMix==PlMix&trans_fab2_sa3_outrm$CanopyHt=="b"),]

  # mean spectra:
  x <- c(350:1350,1450:1800,1980:2200)
  plot(x=x, y=a_m[i,7:ncol(a_m)], ylim=c(0,1.6), type="n", bty="l",
       xlab="wavelength (nm)", ylab="Transmittance")

  # adding lines for mean spectrum:
  lines(x=x, y=a_m[1,7:ncol(a_m)], lwd=3, col="#FBC02D")
  lines(x=x, y=a_b[1,7:ncol(a_m)], lwd=3, col="#0D47A1") # "#FBC02D", "##00BCD4"
  title(PlMix)

  # adding lines for spectral outliers
  for(i in 1:nrow(aa_m)){lines(x=x, y=aa_m[i,13:(ncol(aa_m)-2)], lwd=1, lty=2, col="#FF990050")} #"#FBC02D99"
  for(i in 1:nrow(aa_b)){lines(x=x, y=aa_b[i,13:(ncol(aa_b)-2)], lwd=1, lty=2, col="#3399FF50")} #"#0D47A199"

  # adding lines for individual spectra:
  for(i in 1:nrow(aaa_m)){lines(x=x, y=aaa_m[i,13:(ncol(aaa_m)-2)], lwd=1, lty=1,  col="#FF9900")}
  for(i in 1:nrow(aaa_b)){lines(x=x, y=aaa_b[i,13:(ncol(aaa_b)-2)], lwd=1, lty=1, col="#3399FF")}

  # blanking out noisy sections of the spectrum (straight lines connecting points)
  rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
  rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)

  box(bty="l")
  legend("topleft", legend=c("middle canopy - mean", "included spectrum", "excluded spectrum", "bottom canopy - mean", "included spectrum", "excluded spectrum"),
         col=c("#FBC02D","#FF9900", "#FF990050", "#0D47A1","#3399FF", "#3399FF50"), lwd=c(3,1,1, 3, 1, 1), lty=c(1,1,2,1,1,2), bty="n", cex=0.7)
}
# dev.off()
