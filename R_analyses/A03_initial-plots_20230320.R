# Using mean transmittance calculated per plot (VIS-SWIR)
# Exploratory figures, analyses

#_______________________________________________________________________________
# Load packages
library(data.table)
library(stringr)

#_______________________________________________________________________________
# Load data
file_dir <- "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance"
# all_spec_ <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230328.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed
all_spec_ <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230616.txt")) # updated file with original (uncorrected) LAI and LAI corrected for within shoot-clumping

    #...community information: SR, diversity
    #...would be nice to add PD, FD measures...
cl_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Cloquet_plots-composition-diversity.csv")
fab_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/FAB1&2_Composition of plots.csv")
fr_comminfo <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/IDENT-Freiburg_plot-lookup.csv")

f1_comminfo <- fab_comminfo[which(fab_comminfo$Experiment=="FAB1"),]
f2_comminfo <- fab_comminfo[which(fab_comminfo$Experiment=="FAB2"),]

cl_spabund <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/Cloquet_plot species abundance matrix.csv")

# add community info for whole site
comm_info <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/All-Sites_plot-composition-diversity-lookup.csv")
all_spec <- merge(all_spec_, comm_info[,c(1,5:14)], by="ID", all.x=T, all.y=F)

#_______________________________________________________________________________
# Transmittance and ratios of transmittance at wavelengths of interest
# 360 (UV), 450 (blue), 660 (red), 730 (far red or red edge), 865 (NIR), 1610 (SWIR), 660/730 (R:FR), 660/865 (R:NIR), 660/1610 (R:SWIR)
# 290-310 for UVR8 (any better estimate of best wavelength to use?) -- UVR8 detects 280-315 nm -- average across this whole region? (BLK-C not very reliable below 300 nm)

all_spec$RtoFR <- all_spec$X660/all_spec$X730 # R:FR
all_spec$RtoNIR <- all_spec$X660/all_spec$X865 # R:NIR
all_spec$RtoSWIR <- all_spec$X660/all_spec$X1610 # R:SWIR
all_spec$FPAR <- rowSums(all_spec[,which(colnames(all_spec)=="X400"):which(colnames(all_spec)=="X700")])/300 # How to calculate FPAR -- does this make sense?

#_______________________________________________________________________________
# Splitting by site 
cl_spec <- all_spec[which(all_spec$Site=="Cloquet"),]
fr_spec <- all_spec[which(all_spec$Site=="Freiburg"),]
f1_spec <- all_spec[which(all_spec$Site=="FAB1"),]
f2_spec <- all_spec[which(all_spec$Site=="FAB2"),]

# # Adding community information (no longer needed as added above)
# cl_spec <- merge(cl_spec_,cl_comminfo[,7:17], by.x="Plot", by.y="BlPl", all.x=T, all.y=F)
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
color_palette1 <- colorRampPalette(c("#B71C1C", "#01579B", "#7E57C2"))(3)
color_palette2 <- colorRampPalette(c("#01579B", "#00BCD4", "#FBC02D"))(3) # "#0D47A1"

# this assigns colors according to rank not actual value
color_palette3 <- colorRampPalette(c("#01579B", "#00BCD4", "#FBC02D"))
all_spec$order = findInterval(all_spec$LAI_specloc, sort(all_spec$LAI_specloc))
all_spec$LAI_color <- color_palette3(nrow(all_spec))[all_spec$order]
color_palette4 <- colorRampPalette(c("#01579B", "#00BCD4", "#FBC02D","#B71C1C"))(4)
#_______________________________________________________________________________
### First: plotting full range transmittance

### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-full-range.pdf")

wl_cols <- c(which(colnames(all_spec)=="X350"):which(colnames(all_spec)=="X2200"))

#1. Transmittance across all plots
x=c(350:1350, 1450:1800, 1980:2200)
plot(x=x, all_spec[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
for(i in 1:nrow(all_spec)){
  lines(x=x, all_spec[i,..wl_cols])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
title(main="Full range transmittance across all plots")

#2. Transmittance, split by site and color code by canopy height
x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(2,2))
plot(x=x, all_spec[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet")
cl_spec_m <- cl_spec[cl_spec$CanopyHt=="m",]
cl_spec_b <- cl_spec[cl_spec$CanopyHt=="b",]
for(i in 1:nrow(cl_spec_m)){lines(x=x, cl_spec_m[i,..wl_cols], col=color_palette2[3])}
for(i in 1:nrow(cl_spec_b)){lines(x=x, cl_spec_b[i,..wl_cols], col=color_palette2[1])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", col=c(color_palette2[3],color_palette2[1]), lwd=2, legend=c("Mid canopy", "Below canopy"), bty="n", cex=0.9)

plot(x=x, all_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Freiburg")
fr_spec_m <- fr_spec[fr_spec$CanopyHt=="m",]
fr_spec_b <- fr_spec[fr_spec$CanopyHt=="b",]
for(i in 1:nrow(fr_spec_m)){lines(x=x, fr_spec_m[i,..wl_cols], col=color_palette2[3])}
for(i in 1:nrow(fr_spec_b)){lines(x=x, fr_spec_b[i,..wl_cols], col=color_palette2[1])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, all_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB1")
f1_spec_m <- f1_spec[f1_spec$CanopyHt=="m",]
f1_spec_b <- f1_spec[f1_spec$CanopyHt=="b",]
for(i in 1:nrow(f1_spec_m)){lines(x=x, f1_spec_m[i,..wl_cols], col=color_palette2[3])}
for(i in 1:nrow(f1_spec_b)){lines(x=x, f1_spec_b[i,..wl_cols], col=color_palette2[1])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, all_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB2")
f2_spec_m <- f2_spec[f2_spec$CanopyHt=="m",]
f2_spec_b <- f2_spec[f2_spec$CanopyHt=="b",]
for(i in 1:nrow(f2_spec_m)){lines(x=x, f2_spec_m[i,..wl_cols], col=color_palette2[3])}
for(i in 1:nrow(f2_spec_b)){lines(x=x, f2_spec_b[i,..wl_cols], col=color_palette2[1])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

#3. Transmittance, split by site and canopy height, color code by LAI
#  "#FBC02D", "#00579B" # yellows=high LAI, blues = low LAI
# Cloquet only:
cl_spec$order = findInterval(cl_spec$LAI_specloc, sort(cl_spec$LAI_specloc))
cl_spec$LAI_color <- color_palette3(nrow(cl_spec))[cl_spec$order]
cl_spec_m <- cl_spec[cl_spec$CanopyHt=="m",]
cl_spec_b <- cl_spec[cl_spec$CanopyHt=="b",]
i<-1
x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy (LAI colorcode)")
for(i in 1:nrow(cl_spec_m)){lines(x=x, cl_spec_m[i,..wl_cols], col=cl_spec_m[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - below canopy (LAI colorcode)")
for(i in 1:nrow(cl_spec_b)){lines(x=x, cl_spec_b[i,..wl_cols], col=cl_spec_b[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

# Freiburg only:
fr_spec$order = findInterval(fr_spec$LAI_specloc, sort(fr_spec$LAI_specloc))
fr_spec$LAI_color <- color_palette3(nrow(fr_spec))[fr_spec$order]
fr_spec_m <- fr_spec[fr_spec$CanopyHt=="m",]
fr_spec_b <- fr_spec[fr_spec$CanopyHt=="b",]

x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, fr_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Freiburg - mid canopy (LAI colorcode)")
for(i in 1:nrow(fr_spec_m)){lines(x=x, fr_spec_m[i,..wl_cols], col=fr_spec_m[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, fr_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Freiburg - bottom canopy (LAI colorcode)")
for(i in 1:nrow(fr_spec_b)){lines(x=x, fr_spec_b[i,..wl_cols], col=fr_spec_b[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

# FAB1 only:
f1_spec$order = findInterval(f1_spec$LAI_specloc, sort(f1_spec$LAI_specloc))
f1_spec$LAI_color <- color_palette3(nrow(f1_spec))[f1_spec$order]
f1_spec_m <- f1_spec[f1_spec$CanopyHt=="m",]
f1_spec_b <- f1_spec[f1_spec$CanopyHt=="b",]

x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, f1_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB1 - mid canopy")
for(i in 1:nrow(fr_spec_m)){lines(x=x, f1_spec_m[i,..wl_cols], col=f1_spec_m[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, f1_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB1 - bottom canopy")
for(i in 1:nrow(f1_spec_b)){lines(x=x, f1_spec_b[i,..wl_cols], col=f1_spec_b[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

# FAB2 only:
f2_spec$order = findInterval(f2_spec$LAI_specloc, sort(f2_spec$LAI_specloc))
f2_spec$LAI_color <- color_palette3(nrow(f2_spec))[f2_spec$order]
f2_spec_m <- f2_spec[f2_spec$CanopyHt=="m",]
f2_spec_b <- f2_spec[f2_spec$CanopyHt=="b",]

x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, f2_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB2 - mid canopy")
for(i in 1:nrow(fr_spec_m)){lines(x=x, f2_spec_m[i,..wl_cols], col=f2_spec_m[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, f2_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB2 - bottom canopy")
for(i in 1:nrow(f2_spec_b)){lines(x=x, f2_spec_b[i,..wl_cols], col=f2_spec_b[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

### dev.off()

#4. Transmittance, split by site and canopy height, color code by 
#               (a) angio/gymno
#               (b) species richness (?)
#               (c) functional diversity (?) 
#               (d) phylogenetic diversity (?) (at least, capture monoculture, congeners, beyond... especially where beyond is all angio or mixtures of angio/gymno)
#               (e) spectral diversity (???) 

### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-full-range-angiogymno.pdf")

# Cloquet only:
cl_spec_m <- cl_spec[cl_spec$CanopyHt=="m",]
cl_spec_b <- cl_spec[cl_spec$CanopyHt=="b",]

cl_spec_m_angio <- cl_spec_m[which(cl_spec_m$AngioGymnoMix=="A"),]
cl_spec_m_gymno <- cl_spec_m[which(cl_spec_m$AngioGymnoMix=="G"),]
cl_spec_m_mix <- cl_spec_m[which(cl_spec_m$AngioGymnoMix=="M"),]

cl_spec_b_angio <- cl_spec_b[which(cl_spec_b$AngioGymnoMix=="A"),]
cl_spec_b_gymno <- cl_spec_b[which(cl_spec_b$AngioGymnoMix=="G"),]
cl_spec_b_mix <- cl_spec_b[which(cl_spec_b$AngioGymnoMix=="M"),]

x=c(350:1350, 1450:1800, 1980:2200)
# par(mfrow=c(1,1))
plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy")
for(i in 1:nrow(cl_spec_m_angio)){lines(x=x, cl_spec_m_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(cl_spec_m_gymno)){lines(x=x, cl_spec_m_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(cl_spec_m_mix)){lines(x=x, cl_spec_m_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")

plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - below canopy")
for(i in 1:nrow(cl_spec_b_angio)){lines(x=x, cl_spec_b_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(cl_spec_b_gymno)){lines(x=x, cl_spec_b_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(cl_spec_b_mix)){lines(x=x, cl_spec_b_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")


# Freiburg only:
fr_spec_m <- fr_spec[fr_spec$CanopyHt=="m",]
fr_spec_b <- fr_spec[fr_spec$CanopyHt=="b",]

fr_spec_m_angio <- fr_spec_m[which(fr_spec_m$AngioGymnoMix=="A"),]
fr_spec_m_gymno <- fr_spec_m[which(fr_spec_m$AngioGymnoMix=="G"),]
fr_spec_m_mix <- fr_spec_m[which(fr_spec_m$AngioGymnoMix=="M"),]

fr_spec_b_angio <- fr_spec_b[which(fr_spec_b$AngioGymnoMix=="A"),]
fr_spec_b_gymno <- fr_spec_b[which(fr_spec_b$AngioGymnoMix=="G"),]
fr_spec_b_mix <- fr_spec_b[which(fr_spec_b$AngioGymnoMix=="M"),]

x=c(350:1350, 1450:1800, 1980:2200)
# par(mfrow=c(1,1))
plot(x=x, fr_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Freiburg - mid canopy")
for(i in 1:nrow(fr_spec_m_angio)){lines(x=x, fr_spec_m_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(fr_spec_m_gymno)){lines(x=x, fr_spec_m_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(fr_spec_m_mix)){lines(x=x, fr_spec_m_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")
plot(x=x, fr_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Freiburg - below canopy")
for(i in 1:nrow(fr_spec_b_angio)){lines(x=x, fr_spec_b_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(fr_spec_b_gymno)){lines(x=x, fr_spec_b_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(fr_spec_b_mix)){lines(x=x, fr_spec_b_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")


# FAB1 only:
f1_spec_m <- f1_spec[f1_spec$CanopyHt=="m",]
f1_spec_b <- f1_spec[f1_spec$CanopyHt=="b",]

f1_spec_m_angio <- f1_spec_m[which(f1_spec_m$AngioGymnoMix=="A"),]
f1_spec_m_gymno <- f1_spec_m[which(f1_spec_m$AngioGymnoMix=="G"),]
f1_spec_m_mix <- f1_spec_m[which(f1_spec_m$AngioGymnoMix=="M"),]

f1_spec_b_angio <- f1_spec_b[which(f1_spec_b$AngioGymnoMix=="A"),]
f1_spec_b_gymno <- f1_spec_b[which(f1_spec_b$AngioGymnoMix=="G"),]
f1_spec_b_mix <- f1_spec_b[which(f1_spec_b$AngioGymnoMix=="M"),]

x=c(350:1350, 1450:1800, 1980:2200)
# par(mf1ow=c(1,1))
plot(x=x, f1_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB1 - mid canopy")
for(i in 1:nrow(f1_spec_m_angio)){lines(x=x, f1_spec_m_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(f1_spec_m_gymno)){lines(x=x, f1_spec_m_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(f1_spec_m_mix)){lines(x=x, f1_spec_m_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")
plot(x=x, f1_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB1 - below canopy")
for(i in 1:nrow(f1_spec_b_angio)){lines(x=x, f1_spec_b_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(f1_spec_b_gymno)){lines(x=x, f1_spec_b_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(f1_spec_b_mix)){lines(x=x, f1_spec_b_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")


# FAB2 only:
f2_spec_m <- f2_spec[f2_spec$CanopyHt=="m",]
f2_spec_b <- f2_spec[f2_spec$CanopyHt=="b",]

f2_spec_m_angio <- f2_spec_m[which(f2_spec_m$AngioGymnoMix=="A"),]
f2_spec_m_gymno <- f2_spec_m[which(f2_spec_m$AngioGymnoMix=="G"),]
f2_spec_m_mix <- f2_spec_m[which(f2_spec_m$AngioGymnoMix=="M"),]

f2_spec_b_angio <- f2_spec_b[which(f2_spec_b$AngioGymnoMix=="A"),]
f2_spec_b_gymno <- f2_spec_b[which(f2_spec_b$AngioGymnoMix=="G"),]
f2_spec_b_mix <- f2_spec_b[which(f2_spec_b$AngioGymnoMix=="M"),]

x=c(350:1350, 1450:1800, 1980:2200)
# par(mfrow=c(1,1))
plot(x=x, f2_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB2 - mid canopy")
for(i in 1:nrow(f2_spec_m_angio)){lines(x=x, f2_spec_m_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(f2_spec_m_gymno)){lines(x=x, f2_spec_m_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(f2_spec_m_mix)){lines(x=x, f2_spec_m_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")
plot(x=x, f2_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("FAB2 - below canopy")
for(i in 1:nrow(f2_spec_b_angio)){lines(x=x, f2_spec_b_angio[i,..wl_cols], col=color_palette1[1])}
for(i in 1:nrow(f2_spec_b_gymno)){lines(x=x, f2_spec_b_gymno[i,..wl_cols], col=color_palette1[2])}
for(i in 1:nrow(f2_spec_b_mix)){lines(x=x, f2_spec_b_mix[i,..wl_cols], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")

### dev.off()

#5. Transmittance, split by site, canopy height, angio/gymno/mixture, and color code by LAI?
# Cloquet only:
cl_spec$order = findInterval(cl_spec$LAI_specloc, sort(cl_spec$LAI_specloc))
cl_spec$LAI_color <- color_palette3(nrow(cl_spec))[cl_spec$order]
cl_spec_m <- cl_spec[cl_spec$CanopyHt=="m",]
cl_spec_b <- cl_spec[cl_spec$CanopyHt=="b",]

x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy - angio")
for(i in 1:nrow(cl_spec_m_angio)){lines(x=x, cl_spec_m_angio[i,..wl_cols], col=cl_spec_m_angio[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy - gymno")
for(i in 1:nrow(cl_spec_m_gymno)){lines(x=x, cl_spec_m_gymno[i,..wl_cols], col=cl_spec_m_gymno[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

plot(x=x, cl_spec[i,..wl_cols], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy mix")
for(i in 1:nrow(cl_spec_m_mix)){lines(x=x, cl_spec_m_mix[i,..wl_cols], col=cl_spec_m_mix[i,LAI_color])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")

#_______________________________________________________________________________
# Second: vector normalised
require(spectrolab)

## Setting up vector normalized data:
cl_spec2 <- cl_spec
colnames(cl_spec2)[wl_cols] <- as.numeric(substring(colnames(cl_spec2),2,5)[wl_cols])
cl_specc <- as_spectra(cl_spec2, meta_idxs =c(1:12, 1586:1607))
cl_specc_vn <- normalize(cl_specc)

cl_spec_vn <- data.frame(cl_specc_vn) # putting metadata back in df

cl_spec_b_vn <- cl_spec_vn[which(cl_spec_vn$CanopyHt=="b"),]
cl_spec_m_vn <- cl_spec_vn[which(cl_spec_vn$CanopyHt=="m"),]

cl_spec_b_vn_angio <- cl_spec_b_vn[which(cl_spec_b_vn$AngioGymnoMix=="A"),]
cl_spec_b_vn_gymno <- cl_spec_b_vn[which(cl_spec_b_vn$AngioGymnoMix=="G"),]
cl_spec_b_vn_mix <- cl_spec_b_vn[which(cl_spec_b_vn$AngioGymnoMix=="M"),]

cl_spec_m_vn_angio <- cl_spec_m_vn[which(cl_spec_m_vn$AngioGymnoMix=="A"),]
cl_spec_m_vn_gymno <- cl_spec_m_vn[which(cl_spec_m_vn$AngioGymnoMix=="G"),]
cl_spec_m_vn_mix <- cl_spec_m_vn[which(cl_spec_m_vn$AngioGymnoMix=="M"),]

### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-full-range-angiogymno-vn.pdf")
x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, cl_spec_b_vn[i,37:1609], type="n", lwd=2, bty="l", ylim=c(0,0.06), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - below canopy - vn")
for(i in 1:nrow(cl_spec_b_vn_angio)){lines(x=x, cl_spec_b_vn_angio[i,37:1609], col=color_palette1[1])}
for(i in 1:nrow(cl_spec_b_vn_gymno)){lines(x=x, cl_spec_b_vn_gymno[i,37:1609], col=color_palette1[2])}
for(i in 1:nrow(cl_spec_b_vn_mix)){lines(x=x, cl_spec_b_vn_mix[i,37:1609], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")

x=c(350:1350, 1450:1800, 1980:2200)
par(mfrow=c(1,1))
plot(x=x, cl_spec_m_vn[i,37:1609], type="n", lwd=2, bty="l", ylim=c(0,0.06), 
     xlab="Wavelength (nm)", ylab="Transmittance")
title("Cloquet - mid canopy - vn")
for(i in 1:nrow(cl_spec_m_vn_angio)){lines(x=x, cl_spec_m_vn_angio[i,37:1609], col=color_palette1[1])}
for(i in 1:nrow(cl_spec_m_vn_gymno)){lines(x=x, cl_spec_m_vn_gymno[i,37:1609], col=color_palette1[2])}
for(i in 1:nrow(cl_spec_m_vn_mix)){lines(x=x, cl_spec_m_vn_mix[i,37:1609], col=color_palette1[3])}
rect(xleft = 1350, xright = 1450, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1980, ybottom = -0.2, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -0.2, ytop = 2000, col="white", border=NA)
box(bty="l")
legend("topleft", lwd=2, col=color_palette1, legend=c("Angiosperm", "Gymnosperm", "Mixture"), bty="n")

#### dev.off()




#_______________________________________________________________________________
### Third: transmittance at select wavelengths or ratios
#1. Transmittance vs LAI

### pdf("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Exploring data/Exploring transmittance/figures/Figures_Canopy-transmittance-select-wavelengths.pdf")

par(mfrow=c(1,1))
plot(all_spec$LAI_allloc, all_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
legend("bottomleft", pch=c(1,2,3,4), legend=c("Cloquet", "Freiburg", "FAB1", "FAB2"), bty="n")
legend("topright", pch=16, col=c(color_palette2[c(1,3)]), legend=c("Bottom canopy", "Top canopy"), bty="n")
plot(all_spec$LAI_allloc, all_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X450, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc, all_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])

## Corrected LAI:

par(mfrow=c(1,1))
plot(all_spec$LAI_allloc_corrected, all_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
legend("bottomleft", pch=c(1,2,3,4), legend=c("Cloquet", "Freiburg", "FAB1", "FAB2"), bty="n")
legend("topright", pch=16, col=c(color_palette2[c(1,3)]), legend=c("Bottom canopy", "Top canopy"), bty="n")
plot(all_spec$LAI_allloc_corrected, all_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements - corrected)",
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements - corrected)",
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X450, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])
plot(all_spec$LAI_allloc_corrected, all_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette2[c(1,3)])[as.factor(all_spec$CanopyHt)], pch=c(1,2,3,4)[as.factor(all_spec$Site)])


# For all sites colorcoded by angio/gymno
# b only, excluding FAB2
# all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"&all_spec$Site!="FAB2"),]

plot(all_spec_b$LAI_allloc, all_spec_b$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)], xlim=c(0,16))
points(all_spec_m$LAI_allloc, all_spec_m$FPAR, bty="l", 
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_b$Site)])
legend("top", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements) - corrected", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)], xlim=c(0,16))
points(all_spec_m$LAI_allloc_corrected, all_spec_m$FPAR, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_b$Site)])
legend("top", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")



par(mfrow=c(3,2))
plot(all_spec_b$LAI_allloc, all_spec_b$X360, bty="l", 
     ylab="Transmittance - UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X360, bty="l", 
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("topright", pch=c(1,2,5, 0, 16, 1), legend=c("Cloquet","FAB1","FAB2","Freiburg", "Mid canopy", "Below canopy"), bty="n")

plot(all_spec_b$LAI_allloc, all_spec_b$X450, bty="l", 
     ylab="Transmittance - blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X450, 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,0)[as.factor(all_spec_m$Site)])
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")

plot(all_spec_b$LAI_allloc, all_spec_b$X660, bty="l", 
     ylab="Transmittance - red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X660, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc, all_spec_b$X730, bty="l", 
     ylab="Transmittance - red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X730, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc, all_spec_b$X865, bty="l", 
     ylab="Transmittance - NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X865, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc, all_spec_b$X1610, bty="l", 
     ylab="Transmittance - SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc, all_spec_m$X1610, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])


par(mfrow=c(3,2))
plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X360, bty="l", 
     ylab="Transmittance - UV (360 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X360, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("topright", pch=c(1,2,5, 0, 16, 1), legend=c("Cloquet","FAB1","FAB2","Freiburg", "Mid canopy", "Below canopy"), bty="n")

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X450, bty="l", 
     ylab="Transmittance - blue (450 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X450, 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,0)[as.factor(all_spec_m$Site)])
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X660, bty="l", 
     ylab="Transmittance - red (660 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X660, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X730, bty="l", 
     ylab="Transmittance - red edge (730 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X730, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X865, bty="l", 
     ylab="Transmittance - NIR (865 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X865, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])

plot(all_spec_b$LAI_allloc_corrected, all_spec_b$X1610, bty="l", 
     ylab="Transmittance - SWIR (1610 nm)", xlab="LAI (average of all measurements - corrected)", 
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
points(all_spec_m$LAI_allloc_corrected, all_spec_m$X1610, bty="l", 
       col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])






par(mfrow=c(3,2))
all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc, all_spec_m$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc, all_spec_b$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")

all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc, all_spec_m$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc, all_spec_b$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")


all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc, all_spec_m$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc, all_spec_b$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")








par(mfrow=c(3,2))
all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc_corrected, all_spec_m$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc_corrected, all_spec_b$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")

all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc_corrected, all_spec_m$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc_corrected, all_spec_b$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")


all_spec_m <- all_spec[which(all_spec$CanopyHt=="m"),]
plot(all_spec_m$LAI_allloc_corrected, all_spec_m$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_m$AngioGymnoMix)], pch=c(16,17,18,15)[as.factor(all_spec_m$Site)])
legend("bottomleft", pch=c(16,17,18,15), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Mid canopy")

all_spec_b <- all_spec[which(all_spec$CanopyHt=="b"),]
plot(all_spec_b$LAI_allloc_corrected, all_spec_b$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements) - corrected", xlim=c(0,10), ylim=c(0,1),
     col=c(color_palette1)[as.factor(all_spec_b$AngioGymnoMix)], pch=c(1,2,5,0)[as.factor(all_spec_b$Site)])
legend("bottomleft", pch=c(1,2,5,0), legend=c("Cloquet","FAB1","FAB2","Freiburg"), bty="n")
legend("topright", lwd=2, col=c(color_palette1), legend=c("Angio", "Gymno","Mix"), bty="n")
title("Below canopy")





#.... just for cloquet, looking at LAI vs transmittance at particular wavelengths color coded by composition
add_plot_info <- function(){
  title("Cloquet")
  legend("topright", col=rep(color_palette1,2), pch=c(1,1,1,16,16,16), bty="n",
         legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
                  "Angio - bottom canopy", "Gymno - bottom canopy", "Mix - bottom canopy"))}
plot(cl_spec$LAI_allloc, cl_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$FPAR, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()

#.... just for Freiburg, looking at LAI vs transmittance at particular wavelengths color coded by composition
add_plot_info <- function(){
  title("Freiburg")
  legend("topright", col=rep(color_palette1,2), pch=c(1,1,1,16,16,16), bty="n",
         legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
                  "Angio - bottom canopy", "Gymno - bottom canopy", "Mix - bottom canopy"))}
plot(fr_spec$LAI_allloc, fr_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$FPAR, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()
plot(fr_spec$LAI_allloc, fr_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(fr_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(fr_spec$CanopyHt)])
add_plot_info()

#.... just for FAB1, looking at LAI vs transmittance at particular wavelengths color coded by composition
add_plot_info <- function(){
  title("FAB1")
  legend("topright", col=rep(color_palette1,2), pch=c(1,1,1,16,16,16), bty="n",
         legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
                  "Angio - bottom canopy", "Gymno - bottom canopy", "Mix - bottom canopy"))}
plot(f1_spec$LAI_allloc, f1_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$FPAR, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()
plot(f1_spec$LAI_allloc, f1_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f1_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f1_spec$CanopyHt)])
add_plot_info()

#.... just for FAB2, looking at LAI vs transmittance at particular wavelengths color coded by composition
add_plot_info <- function(){
  title("FAB2")
  legend("topright", col=rep(color_palette1,2), pch=c(1,1,1,16,16,16), bty="n",
         legend=c("Angio - mid canopy", "Gymno - mid canopy", "Mix - mid canopy", 
                  "Angio - bottom canopy", "Gymno - bottom canopy", "Mix - bottom canopy"))}
plot(f2_spec$LAI_allloc, f2_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$FPAR, bty="l", 
     ylab="FPAR", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$FPAR, bty="l", 
     ylab="blue (450 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()
plot(f2_spec$LAI_allloc, f2_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI (average of all measurements)", 
     col=c(color_palette1)[as.factor(f2_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(f2_spec$CanopyHt)])
add_plot_info()

### dev.off()


#2. Transmittance vs FPAR (??)
plot(cl_spec$FPAR, cl_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="FPAR", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$FPAR, cl_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="FPAR", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$FPAR, cl_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="FPAR", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$FPAR, cl_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="FPAR", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$FPAR, cl_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="FPAR", 
     col=c(color_palette1)[as.factor(cl_spec$AngioGymnoMix)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()


#3. Transmittance vs MonoCongenConspDiv (??)
add_plot_info <- function(){
  title("Cloquet")
  legend("topright", col=rep(color_palette4,2), pch=c(1,1,1,1,16,16,16,16), bty="n",
         legend=c("Mono - mid canopy", "Congen - mid canopy", "Con-group - mid canopy","Diff-group - mid canopy", 
                  "Mono - bottom canopy", "Congen - bottom canopy", "Con-group - bottom canopy","Diff-group - bottom canopy"), cex=0.8)}

plot(cl_spec$LAI_allloc, cl_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv.x)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$FPAR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()

#3b. Transmittance vs MonoCongenConspDiv -- ANGIOS ONLY
add_plot_info <- function(){
  title("Cloquet")
  legend("topright", col=rep(color_palette4,2), pch=c(1,1,1,1,16,16,16,16), bty="n",
         legend=c("Mono - mid canopy", "Congen - mid canopy", "Con-group - mid canopy","Diff-group - mid canopy", 
                  "Mono - bottom canopy", "Congen - bottom canopy", "Con-group - bottom canopy","Diff-group - bottom canopy"), cex=0.8)}

cl_spec_angio <- cl_spec[which(cl_spec$AngioGymnoMix=="A"),]
cl_spec_gymno <- cl_spec[which(cl_spec$AngioGymnoMix=="G"),]
cl_spec_mix <- cl_spec[which(cl_spec$AngioGymnoMix=="M"),]

plot(cl_spec_angio$LAI_allloc, cl_spec_angio$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$FPAR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()
plot(cl_spec_angio$LAI_allloc, cl_spec_angio$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_angio$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_angio$CanopyHt)])
add_plot_info()


#3c. Transmittance vs MonoCongenConspDiv -- Gymno ONLY
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$FPAR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()
plot(cl_spec_gymno$LAI_allloc, cl_spec_gymno$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_gymno$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_gymno$CanopyHt)])
add_plot_info()


#3d. Transmittance vs MonoCongenConspDiv -- Mixtures ONLY
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$FPAR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec_mix$MonoCongenCongpDiv)], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()
plot(cl_spec_mix$LAI_allloc, cl_spec_mix$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI", 
     col=c(color_palette4)[cl_spec_mix$MonoCongenCongpDiv], pch=c(16,1)[as.factor(cl_spec_mix$CanopyHt)])
add_plot_info()



#4. Transmittance vs mpd
cl_spec[is.na(cl_spec$mpd),]$mpd <- 0 # assigning monocultures 0
add_plot_info <- function(){
  title("Cloquet")
  legend("topright", col=rep(color_palette4,2), pch=c(1,1,1,1,16,16,16,16), bty="n",
         legend=c("Mono - mid canopy", "Congen - mid canopy", "Con-group - mid canopy","Diff-group - mid canopy", 
                  "Mono - bottom canopy", "Congen - bottom canopy", "Con-group - bottom canopy","Diff-group - bottom canopy"), cex=0.8)}

plot(cl_spec$LAI_allloc, cl_spec$RtoFR, bty="l", 
     ylab="R:FR", xlab="LAI (average of all measurements)", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoNIR, bty="l", 
     ylab="R:NIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$RtoSWIR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$FPAR, bty="l", 
     ylab="R:SWIR", xlab="LAI (average of all measurements)",
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X360, bty="l", 
     ylab="UV (360 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X660, bty="l", 
     ylab="red (660 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X730, bty="l", 
     ylab="red edge (730 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X865, bty="l", 
     ylab="NIR (865 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()
plot(cl_spec$LAI_allloc, cl_spec$X1610, bty="l", 
     ylab="SWIR (1610 nm)", xlab="LAI", 
     col=c(color_palette4)[(cl_spec$mpd)], pch=c(16,1)[as.factor(cl_spec$CanopyHt)])
add_plot_info()




# Predictor variables (x-axes) of diversity
# How to interpret overyielding of transmittance? 
# Expect mixtures would transmit less light than mean of monocultures (assuming mixtures tend to overyield in biomass and leaf area)

# Could look at relationships between transmittance and aboveground biomass (or proxy)? -- what would this tell us?


#_______________________________________________________________________________
# Analyses of relationships to formalize trends observed in figures
## Statistical model to formalize these figures? 
## Can relationship be linearized, or should we fit non-linear function?





