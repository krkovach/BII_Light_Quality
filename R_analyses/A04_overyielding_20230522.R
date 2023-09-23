# Using mean transmittance calculated per plot (VIS-SWIR)
# Calculating overyielding in transmittance
# Overyielding in transmittance cf

#_______________________________________________________________________________
# Load packages
library(data.table)
library(stringr)
require(spectrolab)

#_______________________________________________________________________________
# Load data
file_dir <- "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance"
# all_spec <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230328.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed
all_spec <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230725.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed

all_spec[which(all_spec$Site=="Cloquet"),]$Plot <- substring(all_spec[which(all_spec$Site=="Cloquet"),]$ID,1,4)
all_spec[which(all_spec$Site=="FAB1"),]$ID <- with(all_spec[which(all_spec$Site=="FAB1"),], paste(Site,Plot,Mix,CanopyHt, sep="_"))
all_spec[which(all_spec$Site=="FAB2"),]$ID <- with(all_spec[which(all_spec$Site=="FAB2"),], paste(Site,Plot,Mix,CanopyHt, sep="_"))

#...community information: SR, diversity
#...would be nice to add PD, FD measures...
cl_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Cloquet_plots-composition-diversity.csv")
fab_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/FAB1&2_Composition of plots.csv")
fr_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Freiburg_plot-lookup.csv")

f1_comminfo <- fab_comminfo[which(fab_comminfo$Experiment=="FAB1"),]
f2_comminfo <- fab_comminfo[which(fab_comminfo$Experiment=="FAB2"),]

cl_spabund <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/Cloquet_plot species abundance matrix.csv") # excludes outermost row of trees
f1_spabund <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/FAB1_plot species abundance matrix2.csv") # excludes outermost row of trees
fr_spabund <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/Freiburg_plot species abundance matrix_inside.csv") # decide what makes more sense -- inside or whole of plot

#_______________________________________________________________________________
# Transmittance and ratios of transmittance at wavelengths of interest
# 360 (UV), 450 (blue), 660 (red), 730 (far red or red edge), 865 (NIR), 1610 (SWIR), 660/730 (R:FR), 660/865 (R:NIR), 660/1610 (R:SWIR)
# 290-310 for UVR8 (any better estimate of best wavelength to use?) -- UVR8 detects 280-315 nm -- average across this whole region? (BLK-C not very reliable below 300 nm)

all_spec$RtoFR <- all_spec$X660/all_spec$X730 # R:FR
all_spec$RtoNIR <- all_spec$X660/all_spec$X865 # R:NIR
all_spec$RtoSWIR <- all_spec$X660/all_spec$X1610 # R:SWIR
all_spec$FPAR <- rowSums(all_spec[,which(colnames(all_spec)=="X400"):which(colnames(all_spec)=="X700")])/300 # How to calculate FPAR -- does this make sense?

#_______________________________________________________________________________
# Adding community information and splitting by site
# add community info for whole site
comm_info <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/All-Sites_plot-composition-diversity-lookup.csv")
all_spec <- merge(all_spec, comm_info[,c(1,5:8)], by="ID", all.x=T, all.y=F)

cl_spec <- all_spec[which(all_spec$Site=="Cloquet"),]
fr_spec <- all_spec[which(all_spec$Site=="Freiburg"),]
f1_spec <- all_spec[which(all_spec$Site=="FAB1"),]
f2_spec <- all_spec[which(all_spec$Site=="FAB2"),]

# 
# #_______________________________________________________________________________
# # # Splitting by site and adding community information
# cl_spec_ <- all_spec[which(all_spec$Site=="Cloquet"),]
# fr_spec_ <- all_spec[which(all_spec$Site=="Freiburg"),]
# f1_spec_ <- all_spec[which(all_spec$Site=="FAB1"),]
# f2_spec_ <- all_spec[which(all_spec$Site=="FAB2"),]
# 
# cl_spec <- merge(cl_spec_,cl_comminfo, by.x="Plot", by.y="BlPl", all.x=T, all.y=F)
# 
# fr_spec_$Plot <- as.numeric(fr_spec_$Plot)
# fr_comminfo$Plot <- as.numeric(fr_comminfo$Plot)
# fr_spec <- merge(fr_spec_,fr_comminfo, by="Plot", all.x=T, all.y=F)
# 
# f1_spec_$Plot <- as.numeric(f1_spec_$Plot)
# f1_comminfo$Plot <- as.numeric(f1_comminfo$Plot)
# f1_spec <- merge(f1_spec_,f1_comminfo, by="Plot", all.x=T, all.y=F)
# 
# f2_spec_$Plot <- as.numeric(f2_spec_$Plot)
# f2_comminfo$Plot <- as.numeric(f2_comminfo$Plot)
# f2_spec <- merge(f2_spec_,f2_comminfo, by="Plot", all.x=T, all.y=F)

#_______________________________________________________________________________
# Calculate "overyielding" of (net biodiversity effect on) transmittance
# (observed transmittance minus expected transmittance [mean of transmittance in monocultures weighted by proportion planted])
# whether to calculate overyielding by block or averaging monocultures across blocks?
# whether to calculate in a neighborhood fashion based on location within plot?

#_______________________________________________________________________________
# Overyielding: Cloquet:
colnames(cl_spabund) <- c("plot_ID", "Block", "Plot", 
                          "PIGL", "PIST", "BEPE", "QURO", "PISY","BEPA", 
                          "LADE", "QURU", "ACSA", "ACPL", "LALA", "PIAB")
cl_monos <- cl_spec[cl_spec$SR==1,] # subsetting spectra for monocultures
which(colnames(cl_monos)=="X350"); which(colnames(cl_monos)=="X2200")
cl_monos_mean <- aggregate(cl_monos[,7:1593], by=list(Mix=cl_monos$Mix, CanopyHt=cl_monos$CanopyHt), FUN="mean")# 1587

cl_spec_expected_wl <- data.frame(matrix(data=NA, nrow=nrow(cl_spec), ncol=length(17:1589)))
for(j in 1:nrow(cl_spec)){
  a <- cl_spec[j,] # choosing plot
  a_comp <- data.frame(t(cl_spabund[which(cl_spabund$plot_ID==a$Plot),c(4:15)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  # monos_a <- monos[monos$BlMix%in%c(paste(substring(a$Plot,1,1),a_comp2$species, sep="_")) 
  #                & monos$CanopyHt==a$CanopyHt,] # subsetting monoculture spectra that match to plot, matching monocultures by BLOCK
  monos_a <- cl_monos_mean[cl_monos_mean$Mix%in%c(a_comp2$species) & cl_monos_mean$CanopyHt==a$CanopyHt,] # subsetting mono spectra (site level mean per species) that match to those spp on plot
  scaled_transmittance <- matrix(nrow=nrow(a_comp2), ncol=length(13:1585)) # rows are each species, columns are wavelengths
  for(i in 1:nrow(a_comp2)){
    scaled_transmittance[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),13:1585]) * a_comp2[i,1] # 1579
  }
  cl_spec_expected_wl[j,] <- colSums(scaled_transmittance)
}
colnames(cl_spec_expected_wl) <- colnames(cl_spec)[17:1589] # assigning column names as wavelengths (11:1583)

## calculating ratios -- calculate expected subtract from observed
cl_spec_expected_wl$RtoFR <- cl_spec_expected_wl$X660/cl_spec_expected_wl$X730 # R:FR
cl_spec_expected_wl$RtoNIR <- cl_spec_expected_wl$X660/cl_spec_expected_wl$X865 # R:NIR
cl_spec_expected_wl$RtoSWIR <- cl_spec_expected_wl$X660/cl_spec_expected_wl$X1610 # R:SWIR
cl_spec_expected_wl$FPAR <- rowSums(cl_spec_expected_wl[,which(colnames(cl_spec_expected_wl)=="X400"):which(colnames(cl_spec_expected_wl)=="X700")])/300

## Calculating overyielding
which(colnames(cl_spec)=="X350"); which(colnames(cl_spec)=="X2200") 
cl_spec_overyielding_wl <- data.frame(cl_spec[,17:1593]) - cl_spec_expected_wl # calculating overyielding as obs-exp (note ratios are in the same columns in each dataframe)
cl_spec_overyielding <- cbind(cl_spec[,1:16], cl_spec_overyielding_wl, cl_spec[,1594:1597]) # adding back plot information (with full plot information this could be to 1599)
cl_spec_overyielding <- cl_spec_overyielding[which(cl_spec_overyielding$SR!=1),]
# 
# ## Calculating vn overyielding
# # Vector normalise the simulated "expected" spectra
# cl_spec2_overyielding_wl2 <- cl_spec_overyielding_wl[,1:1573]
# colnames(cl_spec2_overyielding_wl2) <- as.numeric(substring(colnames(cl_spec_overyielding_wl)[1:1573],2,5))
# cl_specc_overyielding_wl <- as_spectra(cl_spec2_overyielding_wl2)
# cl_specc_overyielding_wl_vn <- normalize(cl_specc_overyielding_wl)
# 
# # Vector normalise the observed spectra
# cl_spec2 <- data.frame(cl_spec[,11:1583])
# colnames(cl_spec2) <- as.numeric(substring(colnames(cl_spec)[11:1583],2,5))
# cl_specc <- as_spectra(cl_spec2)
# cl_specc_vn <- normalize(cl_specc)
# cl_specc_vn2 <- data.frame(cl_specc_vn) # what does the normalisation magnitude mean?
# 
# # NBE in vn spectra
# cl_specc_vn_overyielding <- data.frame(cl_specc_vn)[,4:1575] - data.frame(cl_specc_overyielding_wl_vn)[,4:1575]
# cl_spec_vn_overyielding <- cbind(data.frame(cl_spec[,1:10]), cl_specc_vn_overyielding, data.frame(cl_spec[,1588:1591])) # adding back plot information
# cl_spec_vn_overyielding <- cl_spec_vn_overyielding[which(cl_spec_vn_overyielding$SR!=1),]

## Cloquet: overyielding of LAI
which(colnames(cl_spec)=="LAI_allloc"); which(colnames(cl_spec)=="LAI_ExclOuter_specloc_corrected") # now 10 variables

cl_lai_expected <- data.frame(matrix(data=NA, nrow=nrow(cl_spec), ncol=10))
for(j in 1:nrow(cl_spec)){
  a <- cl_spec[j,] # choosing plot
  a_comp <- data.frame(t(cl_spabund[which(cl_spabund$plot_ID==a$Plot),c(4:15)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  # monos_a <- monos[monos$BlMix%in%c(paste(substring(a$Plot,1,1),a_comp2$species, sep="_")) 
  #                  & monos$CanopyHt==a$CanopyHt,] # subsetting monoculture spectra that match to plot, matching monocultures by BLOCK
  monos_a <- cl_monos_mean[cl_monos_mean$Mix%in%c(a_comp2$species) & cl_monos_mean$CanopyHt==a$CanopyHt,] # subsetting mono spectra (site level mean per species) that match to those spp on plot
  scaled_transmittance <- matrix(nrow=nrow(a_comp2), ncol=10) # rows are each species, columns are wavelengths
  for(i in 1:nrow(a_comp2)){
    scaled_transmittance[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),3:12]) * a_comp2[i,1]
  }
  cl_lai_expected[j,] <- colSums(scaled_transmittance)
}

cl_lai_overyielding_only <- data.frame(cl_spec[,7:16]) - cl_lai_expected # calculating overyielding as obs-exp
colnames(cl_lai_overyielding_only) <- paste("NBE_", colnames(cl_lai_overyielding_only),  sep="")
cl_lai_overyielding <- cbind(cl_spec[,1:16], cl_lai_overyielding_only, cl_spec[,1594:1597]) # adding back plot information
cl_spec_overyielding2 <- merge(cl_spec_overyielding, cl_lai_overyielding[,c(1,17:26)], by="ID") # merging columns of OY in LAI
### CHECK TO SEE IF THIS IS ADDING MORE THAN WHAT WANT...?

#_______________________________________________________________________________
# Overyielding: FAB1:
# Note that we don't have mid canopy transmittance values for monocultures of all species in FAB1 (missing ACNE, ACRU, QUEL)
test <- f1_spabund[which(f1_spabund$Plot%in%f1_spec$Plot),]
test[which(test$ACNE==0&test$ACRU==0&test$QUEL==0),] # only two mixed species plots do not have ACNE, ACRU or QUEL

# Thus calculating overyielding solely for below canopy
f1_spec_b <- f1_spec[which(f1_spec$CanopyHt=="b"),]
colnames(f1_spabund)[6:17] <- c("ACNE", "ACRU", "BEPA", "JUVI", "PIBA", "PIRE", "PIST", "QUAL", "QUEL", "QUMA", "QURU", "TIAM")
f1_monos <- f1_spec_b[f1_spec_b$SR==1,] # subsetting spectra for monocultures
f1_monos_mean <- aggregate(f1_monos[,7:1593], by=list(Mix=f1_monos$Mix), FUN="mean") #7:1587 (this is )

f1_spec_expected_wl <- data.frame(matrix(data=NA, nrow=nrow(f1_spec_b), ncol=1573))
for(j in 1:nrow(f1_spec_b)){
  a <- f1_spec_b[j,] # choosing plot
  a_comp <- data.frame(t(f1_spabund[which(f1_spabund$Plot==a$Plot),c(6:17)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  # monos_a <- monos[monos$BlMix%in%c(paste(substring(a$Plot,1,1),a_comp2$species, sep="_"))
  #                & monos$CanopyHt==a$CanopyHt,] # subsetting monoculture spectra that match to plot, matching monocultures by BLOCK
  monos_a <- f1_monos_mean[f1_monos_mean$Mix%in%c(a_comp2$species),] # subsetting mono spectra (site level mean per species) that match to those spp on plot
  scaled_transmittance <- matrix(nrow=nrow(a_comp2), ncol=1573) # rows are each species, columns are wavelengths
  for(i in 1:nrow(a_comp2)){
    scaled_transmittance[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),12:1584]) * a_comp2[i,1] # 13:1585
  }
  f1_spec_expected_wl[j,] <- colSums(scaled_transmittance)
  print(j)
}
which(colnames(f1_spec)=="X350"); which(colnames(f1_spec)=="X2200")
colnames(f1_spec_expected_wl) <- colnames(f1_spec)[17:1589] # assigning column names as wavelengths

## calculating ratios -- calculate expected subtract from observed
f1_spec_expected_wl$RtoFR <- f1_spec_expected_wl$X660/f1_spec_expected_wl$X730 # R:FR
f1_spec_expected_wl$RtoNIR <- f1_spec_expected_wl$X660/f1_spec_expected_wl$X865 # R:NIR
f1_spec_expected_wl$RtoSWIR <- f1_spec_expected_wl$X660/f1_spec_expected_wl$X1610 # R:SWIR
f1_spec_expected_wl$FPAR <- rowSums(f1_spec_expected_wl[,which(colnames(f1_spec_expected_wl)=="X400"):which(colnames(f1_spec_expected_wl)=="X700")])/300

## Calculating spectral overyielding
which(colnames(f1_spec_b)=="X350"); which(colnames(f1_spec_b)=="FPAR")
f1_spec_overyielding_wl <- data.frame(f1_spec_b[,17:1593]) - f1_spec_expected_wl # calculating overyielding as obs-exp
f1_spec_overyielding <- cbind(f1_spec_b[,1:16], f1_spec_overyielding_wl, f1_spec_b[,1594:1597]) # adding back plot information
f1_spec_overyielding <- f1_spec_overyielding[which(f1_spec_overyielding$SR!=1),] # only 12 mixed species plots

## Overyielding of LAI
f1_lai_expected <- data.frame(matrix(data=NA, nrow=nrow(f1_spec_b), ncol=10))
for(j in 1:nrow(f1_spec_b)){
  a <- f1_spec_b[j,] # choosing plot
  a_comp <- data.frame(t(f1_spabund[which(f1_spabund$Plot==a$Plot),c(6:17)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  monos_a <- f1_monos_mean[f1_monos_mean$Mix%in%c(a_comp2$species),] # subsetting mono LAI (site level mean per species) that match to those spp on plot
  scaled_transmittance <- matrix(nrow=nrow(a_comp2), ncol=10) # rows are each species, columns are LAI metrics
  for(i in 1:nrow(a_comp2)){
    scaled_transmittance[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),2:11]) * a_comp2[i,1] # choosing columns for LAI metrics
  }
  f1_lai_expected[j,] <- colSums(scaled_transmittance)
}

f1_lai_overyielding_only <- data.frame(f1_spec_b[,7:16]) - f1_lai_expected # calculating overyielding as obs-exp
colnames(f1_lai_overyielding_only) <- paste("NBE_", colnames(f1_lai_overyielding_only),  sep="")
f1_lai_overyielding <- cbind(f1_spec_b[,1:6], f1_lai_overyielding_only, f1_spec_b[,1594:1597]) # adding back plot information
f1_spec_overyielding2 <- merge(f1_spec_overyielding, f1_lai_overyielding[,c(1,7:16)], by="ID") # merging columns of OY in LAI


#_______________________________________________________________________________
# Overyielding: Freiburg:
fr_monos <- fr_spec[fr_spec$SR==1,] # subsetting spectra for monocultures
fr_monos_mean <- aggregate(fr_monos[,7:1593], by=list(Mix=fr_monos$Mix, CanopyHt=fr_monos$CanopyHt), FUN="mean")

fr_spec_expected_wl <- data.frame(matrix(data=NA, nrow=nrow(fr_spec), ncol=1573))
for(j in 1:nrow(fr_spec)){
  a <- fr_spec[j,] # choosing plot
  a_comp <- data.frame(t(fr_spabund[which(fr_spabund$Plot==a$Plot),c(2:13)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  monos_a <- fr_monos_mean[fr_monos_mean$Mix%in%c(a_comp2$species) & fr_monos_mean$CanopyHt==a$CanopyHt,] # subsetting mono spectra (site level mean per species) that match to those spp on plot
  scaled_transmittance <- matrix(nrow=nrow(a_comp2), ncol=1573) # rows are each species, columns are wavelengths
  for(i in 1:nrow(a_comp2)){
    scaled_transmittance[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),13:1585]) * a_comp2[i,1]
  }
  fr_spec_expected_wl[j,] <- colSums(scaled_transmittance)
}
colnames(fr_spec_expected_wl) <- colnames(fr_spec)[17:1589] # assigning column names as wavelengths

## calculating ratios -- calculate expected subtract from observed
fr_spec_expected_wl$RtoFR <- fr_spec_expected_wl$X660/fr_spec_expected_wl$X730 # R:FR
fr_spec_expected_wl$RtoNIR <- fr_spec_expected_wl$X660/fr_spec_expected_wl$X865 # R:NIR
fr_spec_expected_wl$RtoSWIR <- fr_spec_expected_wl$X660/fr_spec_expected_wl$X1610 # R:SWIR
fr_spec_expected_wl$FPAR <- rowSums(fr_spec_expected_wl[,which(colnames(fr_spec_expected_wl)=="X400"):which(colnames(fr_spec_expected_wl)=="X700")])/300

## Calculating spectral overyielding
fr_spec_overyielding_wl <- data.frame(fr_spec[,17:1593]) - fr_spec_expected_wl # calculating overyielding as obs-exp (note ratios are in the same columns in each dataframe)
fr_spec_overyielding <- cbind(fr_spec[,1:10], fr_spec_overyielding_wl, fr_spec[,1594:1597]) # adding back plot information (with full plot information this could be to 1599)
fr_spec_overyielding <- fr_spec_overyielding[which(fr_spec_overyielding$SR!=1),]


## Overyielding of LAI
fr_lai_expected <- data.frame(matrix(data=NA, nrow=nrow(fr_spec), ncol=10))
for(j in 1:nrow(fr_spec)){
  a <- fr_spec[j,] # choosing plot
  a_comp <- data.frame(t(fr_spabund[which(fr_spabund$Plot==a$Plot),c(2:13)])) # subsetting composition matrix for plot
  a_comp$species <- rownames(a_comp)
  a_comp2 <- a_comp[which(a_comp[,1]>0),] # subsetting composition to only those species present in plot
  monos_a <- cl_monos_mean[fr_monos_mean$Mix%in%c(a_comp2$species) & fr_monos_mean$CanopyHt==a$CanopyHt,] # subsetting mono spectra (site level mean per species) that match to those spp on plot
  scaled_lai <- matrix(nrow=nrow(a_comp2), ncol=10) # rows are each species, columns are lai
    for(i in 1:nrow(a_comp2)){
      scaled_lai[i,] <- t(monos_a[which(monos_a$Mix==a_comp2$species[i]),3:12]) * a_comp2[i,1] # choosing columns for LAI metrics
  }
  fr_lai_expected[j,] <- colSums(scaled_lai)
}

fr_lai_overyielding_only <- data.frame(fr_spec[,7:16]) - fr_lai_expected # calculating overyielding as obs-exp
colnames(fr_lai_overyielding_only) <- paste("NBE_", colnames(fr_lai_overyielding_only),  sep="")
fr_lai_overyielding <- cbind(fr_spec[,1:6], fr_lai_overyielding_only, fr_spec[,1594:1597]) # adding back plot information
fr_spec_overyielding2 <- merge(fr_spec_overyielding, fr_lai_overyielding[,c(1,7:16)], by="ID") # merging columns of OY in LAI

###* CALCULATE STEM GROWTH AND NBE STEM GROWTH FOR AMAP (AT LEAST CLOQUET AND FAB1) -- does relationship to FPAR differ from whole spectrum or??
# 
# #_______________________________________________________________________________
# # Plotting overyielding of transmittance
# # (does this make sense -- if want to interpret this as the transmittance we would "expect" for a mixture, it assumes transmittance should be additive??)
# # (need to think about how to interpret this comparison)
# 
# #1. Look at overyielding in transmittance across spectrum
# # splitting dataframes:
# cl_spec_overyielding_m <- cl_spec_overyielding[which(cl_spec_overyielding$CanopyHt=="m"),]
# cl_spec_overyielding_b <- cl_spec_overyielding[which(cl_spec_overyielding$CanopyHt=="b"),]
# cl_spec_overyielding_m_angio <- cl_spec_overyielding_m[which(cl_spec_overyielding_m$AngioGymnoMix=="A"),]
# cl_spec_overyielding_m_gymno <- cl_spec_overyielding_m[which(cl_spec_overyielding_m$AngioGymnoMix=="G"),]
# cl_spec_overyielding_m_mix <- cl_spec_overyielding_m[which(cl_spec_overyielding_m$AngioGymnoMix=="M"),]
# cl_spec_overyielding_b_angio <- cl_spec_overyielding_b[which(cl_spec_overyielding_b$AngioGymnoMix=="A"),]
# cl_spec_overyielding_b_gymno <- cl_spec_overyielding_b[which(cl_spec_overyielding_b$AngioGymnoMix=="G"),]
# cl_spec_overyielding_b_mix <- cl_spec_overyielding_b[which(cl_spec_overyielding_b$AngioGymnoMix=="M"),] 
# 
# f1_spec_overyielding_b_angio <- f1_spec_overyielding[which(f1_spec_overyielding$AngioGymnoMix=="A"),] # note only below canopy values for FAB1
# f1_spec_overyielding_b_gymno <- f1_spec_overyielding[which(f1_spec_overyielding$AngioGymnoMix=="G"),]
# f1_spec_overyielding_b_mix <- f1_spec_overyielding[which(f1_spec_overyielding$AngioGymnoMix=="M"),] 
# 
# fr_spec_overyielding_m <- fr_spec_overyielding[which(fr_spec_overyielding$CanopyHt=="m"),]
# fr_spec_overyielding_b <- fr_spec_overyielding[which(fr_spec_overyielding$CanopyHt=="b"),]
# fr_spec_overyielding_m_angio <- fr_spec_overyielding_m[which(fr_spec_overyielding_m$AngioGymnoMix=="A"),]
# fr_spec_overyielding_m_gymno <- fr_spec_overyielding_m[which(fr_spec_overyielding_m$AngioGymnoMix=="G"),]
# fr_spec_overyielding_m_mix <- fr_spec_overyielding_m[which(fr_spec_overyielding_m$AngioGymnoMix=="M"),]
# fr_spec_overyielding_b_angio <- fr_spec_overyielding_b[which(fr_spec_overyielding_b$AngioGymnoMix=="A"),]
# fr_spec_overyielding_b_gymno <- fr_spec_overyielding_b[which(fr_spec_overyielding_b$AngioGymnoMix=="G"),]
# fr_spec_overyielding_b_mix <- fr_spec_overyielding_b[which(fr_spec_overyielding_b$AngioGymnoMix=="M"),] 
# 
# #### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-overyielding-full-range-angiogymno.pdf")
# color_palette1 <- colorRampPalette(c("#B71C1C", "#01579B", "#7E57C2"))(3)
# x=c(350:1350, 1450:1800, 1980:2200)
# par(mfrow=c(3,2))
# #...cloquet
# plot(x=x, cl_spec_overyielding_m[1,11:1583], type="n", lwd=2, bty="l", ylim=c(-0.5,1.2), 
#      xlab="Wavelength (nm)", ylab="NBE transmittance")
# for(i in 1:nrow(cl_spec_overyielding_m_angio)){lines(x=x, cl_spec_overyielding_m_angio[i,11:1583], col=color_palette1[1])} # Colorcode by angio/gymno
# for(i in 1:nrow(cl_spec_overyielding_m_gymno)){lines(x=x, cl_spec_overyielding_m_gymno[i,11:1583], col=color_palette1[2])}
# for(i in 1:nrow(cl_spec_overyielding_m_mix)){lines(x=x, cl_spec_overyielding_m_mix[i,11:1583], col=color_palette1[3])}
# rect(xleft = 1350, xright = 1450, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 1800, xright = 1980, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 2200, xright = 2505, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# box(bty="l")
# abline(h=0, lwd=2, lty=2, col="gray")
# legend("topleft", col=color_palette1, legend=c("Angio", "Gymno", "Mix"), bty="n", lwd=2)
# title(main="NBE transmittance - Cloquet - mid canopy")
# 
# plot(x=x, cl_spec_overyielding_b[1,11:1583], type="n", lwd=2, bty="l", ylim=c(-0.4,0.4), 
#      xlab="Wavelength (nm)", ylab="NBE transmittance")
# for(i in 1:nrow(cl_spec_overyielding_b_angio)){lines(x=x, cl_spec_overyielding_b_angio[i,11:1583], col=color_palette1[1])} # Colorcode by angio/gymno
# for(i in 1:nrow(cl_spec_overyielding_b_gymno)){lines(x=x, cl_spec_overyielding_b_gymno[i,11:1583], col=color_palette1[2])}
# for(i in 1:nrow(cl_spec_overyielding_b_mix)){lines(x=x, cl_spec_overyielding_b_mix[i,11:1583], col=color_palette1[3])}
# rect(xleft = 1350, xright = 1450, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 1800, xright = 1980, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 2200, xright = 2505, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# box(bty="l")
# abline(h=0, lwd=2, lty=2, col="gray")
# legend("topleft", col=color_palette1, legend=c("Angio", "Gymno", "Mix"), bty="n", lwd=2)
# title(main="NBE transmittance - Cloquet - below canopy")
# 
# #...freiburg
# plot(x=x, fr_spec_overyielding_m[1,11:1583], type="n", lwd=2, bty="l", ylim=c(-0.5,1.2), 
#      xlab="Wavelength (nm)", ylab="NBE transmittance")
# for(i in 1:nrow(fr_spec_overyielding_m_angio)){lines(x=x, fr_spec_overyielding_m_angio[i,11:1583], col=color_palette1[1])} # Colorcode by angio/gymno
# for(i in 1:nrow(fr_spec_overyielding_m_gymno)){lines(x=x, fr_spec_overyielding_m_gymno[i,11:1583], col=color_palette1[2])}
# for(i in 1:nrow(fr_spec_overyielding_m_mix)){lines(x=x, fr_spec_overyielding_m_mix[i,11:1583], col=color_palette1[3])}
# rect(xleft = 1350, xright = 1450, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 1800, xright = 1980, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 2200, xright = 2505, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# box(bty="l")
# abline(h=0, lwd=2, lty=2, col="gray")
# legend("topleft", col=color_palette1, legend=c("Angio", "Gymno", "Mix"), bty="n", lwd=2)
# title(main="NBE transmittance - Freiburg - mid canopy")
# 
# plot(x=x, fr_spec_overyielding_b[1,11:1583], type="n", lwd=2, bty="l", ylim=c(-0.4,0.4), 
#      xlab="Wavelength (nm)", ylab="NBE transmittance")
# for(i in 1:nrow(fr_spec_overyielding_b_angio)){lines(x=x, fr_spec_overyielding_b_angio[i,11:1583], col=color_palette1[1])} # Colorcode by angio/gymno
# for(i in 1:nrow(fr_spec_overyielding_b_gymno)){lines(x=x, fr_spec_overyielding_b_gymno[i,11:1583], col=color_palette1[2])}
# for(i in 1:nrow(fr_spec_overyielding_b_mix)){lines(x=x, fr_spec_overyielding_b_mix[i,11:1583], col=color_palette1[3])}
# rect(xleft = 1350, xright = 1450, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 1800, xright = 1980, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 2200, xright = 2505, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# box(bty="l")
# abline(h=0, lwd=2, lty=2, col="gray")
# legend("topleft", col=color_palette1, legend=c("Angio", "Gymno", "Mix"), bty="n", lwd=2)
# title(main="NBE transmittance - Freiburg - below canopy")
# 
# #...fab1
# # <-- FD would be more meaningful than A/G/M for OY (have only one con-group mixture)
# plot(x=0,y=1,type="n", axes=F, ylab="", xlab="") # note no mid canopy for fab1, so adding blank space to keep panels lined up
# 
# plot(x=x, f1_spec_overyielding[1,11:1583], type="n", lwd=2, bty="l", ylim=c(-0.5,0.4), 
#      xlab="Wavelength (nm)", ylab="NBE transmittance")
# for(i in 1:nrow(f1_spec_overyielding_b_angio)){lines(x=x, f1_spec_overyielding_b_angio[i,11:1583], col=color_palette1[1])} # Colorcode by angio/gymno
# for(i in 1:nrow(f1_spec_overyielding_b_gymno)){lines(x=x, f1_spec_overyielding_b_gymno[i,11:1583], col=color_palette1[2])}
# for(i in 1:nrow(f1_spec_overyielding_b_mix)){lines(x=x, f1_spec_overyielding_b_mix[i,11:1583], col=color_palette1[3])}
# rect(xleft = 1350, xright = 1450, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 1800, xright = 1980, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# rect(xleft = 2200, xright = 2505, ybottom = -0.6, ytop = 2000, col="white", border=NA)
# box(bty="l")
# abline(h=0, lwd=2, lty=2, col="gray")
# legend("topleft", col=color_palette1, legend=c("Angio", "Gymno", "Mix"), bty="n", lwd=2)
# title(main="NBE transmittance - FAB1 - below canopy")
# 
# 
# #### dev.off()
# 
# #2. OY transmittance vs LAI
# #looking at LAI vs transmittance at particular wavelengths color coded by composition
# 
# ### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-overyielding-select-wavelengths.pdf")
# 
# #...Cloquet
# par(mfrow=c(3,3))
# add_plot_info <- function(add_legend=F){
#   title("Cloquet")
#   if(add_legend==T){legend("topright", col=rep(color_palette1,2), pch=c(16,16,16,1,1,1), bty="n",
#          legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
#                   "Angio - below canopy", "Gymno - below canopy", "Mix - below canopy"))}
#   abline(h=0, lty=2, col="grey")
#   abline(v=0, lty=2, col="grey")}
# 
# cl_spec_overyielding$NBE_FPAR <- rowSums(cl_spec_overyielding[,61:361])/300 
# 
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$RtoFR, bty="l", 
#      ylab="NBE R:FR", xlab="LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$RtoNIR, bty="l",
#      ylab="NBE R:NIR", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$RtoSWIR, bty="l",
#      ylab="NBE R:SWIR", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$FPAR, bty="l",
#      ylab="NBE FPAR", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$X360, bty="l",
#      ylab="NBE UV (360 nm)", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$X660, bty="l",
#      ylab="NBE red (660 nm)", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$X730, bty="l",
#      ylab="NBE red edge (730 nm)", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$X865, bty="l",
#      ylab="NBE NIR (865 nm)", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding$LAI_allloc, cl_spec_overyielding$X1610, bty="l",
#      ylab="NBE SWIR (1610 nm)", xlab="LAI (average of all measurements)",
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info(add_legend=T)
# 
# 
# #2. OY transmittance vs OY LAI
# # cl_spec_overyielding2$NBE_FPAR <- rowMeans(cl_spec_overyielding2[,61:361]) 
# 
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$RtoFR, bty="l", 
#      ylab="NBE R:FR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$RtoNIR, bty="l", 
#      ylab="NBE R:NIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$RtoSWIR, bty="l", 
#      ylab="NBE R:SWIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(x=cl_spec_overyielding2$NBE_LAI_allloc, y=cl_spec_overyielding2$FPAR, bty="l", 
#      ylab="NBE FPAR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$X360, bty="l", 
#      ylab="NBE UV (360 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$X660, bty="l", 
#      ylab="NBE red (660 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$X730, bty="l", 
#      ylab="NBE red edge (730 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$X865, bty="l", 
#      ylab="NBE NIR (865 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(cl_spec_overyielding2$NBE_LAI_allloc, cl_spec_overyielding2$X1610, bty="l", 
#      ylab="NBE SWIR (1610 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(cl_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(cl_spec_overyielding2$CanopyHt)])
# add_plot_info(add_legend = T)
# 
# ## dev.off()
# 
# # Freiburg
# #2. OY transmittance vs OY LAI
# par(mfrow=c(3,3))
# add_plot_info <- function(add_legend=F){
#   title("Freiburg")
#   if(add_legend==T){legend("topright", col=rep(color_palette1,2), pch=c(16,16,16,1,1,1), bty="n",
#                            legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
#                                     "Angio - below canopy", "Gymno - below canopy", "Mix - below canopy"))}
#   abline(h=0, lty=2, col="grey")
#   abline(v=0, lty=2, col="grey")}
# 
# fr_spec_overyielding2$NBE_FPAR <- rowSums(fr_spec_overyielding2[,61:361])/300 
# 
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$RtoFR, bty="l", 
#      ylab="NBE R:FR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$RtoNIR, bty="l", 
#      ylab="NBE R:NIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$RtoSWIR, bty="l", 
#      ylab="NBE R:SWIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(x=fr_spec_overyielding2$NBE_LAI_allloc, y=fr_spec_overyielding2$FPAR, bty="l", 
#      ylab="NBE FPAR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$X360, bty="l", 
#      ylab="NBE UV (360 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$X660, bty="l", 
#      ylab="NBE red (660 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$X730, bty="l", 
#      ylab="NBE red edge (730 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$X865, bty="l", 
#      ylab="NBE NIR (865 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(fr_spec_overyielding2$NBE_LAI_allloc, fr_spec_overyielding2$X1610, bty="l", 
#      ylab="NBE SWIR (1610 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1)[as.factor(fr_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(fr_spec_overyielding2$CanopyHt)])
# add_plot_info(add_legend = T)
# 
# ## dev.off()
# 
# 
# 
# # FAB1
# #2. OY transmittance vs OY LAI
# par(mfrow=c(3,3))
# add_plot_info <- function(add_legend=F){
#   title("FAB1")
#   if(add_legend==T){legend("topright", col=rep(color_palette1,2), pch=c(1,1,1), bty="n",
#                            legend=c(
#                                     "Angio - below canopy", "Gymno - below canopy", "Mix - below canopy"))}
#   abline(h=0, lty=2, col="grey")
#   abline(v=0, lty=2, col="grey")}
# 
# f1_spec_overyielding2$NBE_FPAR <- rowSums(f1_spec_overyielding2[,61:361])/300 
# 
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$RtoFR, bty="l", 
#      ylab="NBE R:FR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$RtoNIR, bty="l", 
#      ylab="NBE R:NIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$RtoSWIR, bty="l", 
#      ylab="NBE R:SWIR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding$CanopyHt)])
# add_plot_info()
# plot(x=f1_spec_overyielding2$NBE_LAI_allloc, y=f1_spec_overyielding2$FPAR, bty="l", 
#      ylab="NBE FPAR", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$X360, bty="l", 
#      ylab="NBE UV (360 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$X660, bty="l", 
#      ylab="NBE red (660 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$X730, bty="l", 
#      ylab="NBE red edge (730 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$X865, bty="l", 
#      ylab="NBE NIR (865 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info()
# plot(f1_spec_overyielding2$NBE_LAI_allloc, f1_spec_overyielding2$X1610, bty="l", 
#      ylab="NBE SWIR (1610 nm)", xlab="NBE LAI (average of all measurements)", 
#      col=c(color_palette1[c(1,3)])[as.factor(f1_spec_overyielding2$AngioGymnoMix)], pch=c(1,16)[as.factor(f1_spec_overyielding2$CanopyHt)])
# add_plot_info(add_legend = T)

## dev.off()




