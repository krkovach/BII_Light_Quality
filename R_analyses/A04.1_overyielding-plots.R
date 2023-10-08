# Plots of overyielding spectra

#_______________________________________________________________________________
# Load data
file_dir <- "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance"
all_spec_overyielding <- fread(file.path(file_dir, "All-Sites_ASD-SVC_NetBiodiversityEffect_20230919.txt"))


#_______________________________________________________________________________
# color_palette2 <- colorRampPalette(c("#01579B", "#00BCD4", "#FBC02D"))(3) # "#0D47A1"
color_palette2 <- c("#164673","#2A9D8F","#E9C46A") # "#264653"

#_______________________________________________________________________________
# Splitting by canopy height
all_spec_overyielding_b <- all_spec_overyielding[which(all_spec_overyielding$CanopyHt=="b"),]
all_spec_overyielding_m <- all_spec_overyielding[which(all_spec_overyielding$CanopyHt=="m"),]

# Plotting each canopy height separately colorcoded by site
cl_spec_overyielding_b <- all_spec_overyielding_b[which(all_spec_overyielding_b$Site=="Cloquet"),]
f1_spec_overyielding_b <- all_spec_overyielding_b[which(all_spec_overyielding_b$Site=="FAB1"),]
fr_spec_overyielding_b <- all_spec_overyielding_b[which(all_spec_overyielding_b$Site=="Freiburg"),]

cl_spec_overyielding_m <- all_spec_overyielding_m[which(all_spec_overyielding_m$Site=="Cloquet"),]
f1_spec_overyielding_m <- all_spec_overyielding_m[which(all_spec_overyielding_m$Site=="FAB1"),]
fr_spec_overyielding_m <- all_spec_overyielding_m[which(all_spec_overyielding_m$Site=="Freiburg"),]

all_spec_overyielding_m[which(all_spec_overyielding_m$X350< -0.5),1:10] # what's going on here? For now just omitting 350 nm but check this
all_spec_overyielding_m[which(all_spec_overyielding_m$X2200< -0.5),]$X2200

x=c(351:1350, 1450:1800, 1980:2200)
wl_cols <- c(which(colnames(all_spec_overyielding)=="X351"):which(colnames(all_spec_overyielding)=="X2200"))

# to save plots:
# pdf(file="/Users/30062322/Documents/Research/UMN/Manuscripts/LightQuality/Figures/NBE_Transmittance_allsites.pdf", width=7.5, height=5)
# below canopy:
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,0.6), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(cl_spec_overyielding_b)){lines(x=x, cl_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[1],90))}
for(i in 1:nrow(f1_spec_overyielding_b)){lines(x=x, f1_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[2],90))}
for(i in 1:nrow(fr_spec_overyielding_b)){lines(x=x, fr_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[3],90))}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots all sites")

# mid canopy:
plot(x=x, all_spec_overyielding_m[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-1,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
# axis(1, at=seq(350,2200, by=50), labels=F, tck=-0.01)
for(i in 1:nrow(cl_spec_overyielding_m)){lines(x=x, cl_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[1],90))}
for(i in 1:nrow(f1_spec_overyielding_m)){lines(x=x, f1_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[2],90))}
for(i in 1:nrow(fr_spec_overyielding_m)){lines(x=x, fr_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[3],90))}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance mid canopy across all plots all sites")

# plot(x=x2, all_spec_overyielding_m[1,..wl_cols2], type="n", lwd=2, bty="l", ylim=c(0,1.5), 
#      xlab="Wavelength (nm)", ylab="Transmittance")
# legend("topleft", lwd=2, legend=c("Cloquet", "FAB1", "Freiburg"), col=color_palette2, bty="n")

# Now each site separately with top and bottom on same figure
# Cloquet
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(cl_spec_overyielding_b)){lines(x=x, cl_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[1],90))}
for(i in 1:nrow(cl_spec_overyielding_m)){lines(x=x, cl_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[3],90))}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Cloquet")

# FAB1
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(f1_spec_overyielding_b)){lines(x=x, f1_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[1],90))}
for(i in 1:nrow(f1_spec_overyielding_m)){lines(x=x, f1_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[3],90))}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-FAB1")

# Freiburg
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(fr_spec_overyielding_b)){lines(x=x, fr_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2[1],90))}
for(i in 1:nrow(fr_spec_overyielding_m)){lines(x=x, fr_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2[3],90))}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Freiburg")



# Now each site separately with top and bottom on separate figures
# Cloquet
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(cl_spec_overyielding_b)){lines(x=x, cl_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(cl_spec_overyielding_b$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Cloquet \n b-colorAGM")



# Now each site separately with top and bottom on separate figures
# Cloquet
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(cl_spec_overyielding_m)){lines(x=x, cl_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(cl_spec_overyielding_m$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Cloquet \n m-colorAGM")


# FAB1
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(f1_spec_overyielding_b)){lines(x=x, f1_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2,90)[c(1,2)][as.factor(f1_spec_overyielding_b$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-FAB1 \n b-colorAGM")


# FAB1
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(f1_spec_overyielding_m)){lines(x=x, f1_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2,90)[c(1,2)][as.factor(f1_spec_overyielding_m$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-FAB1 \n m-colorAGM")



# Freiburg
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(fr_spec_overyielding_b)){lines(x=x, fr_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(fr_spec_overyielding_b$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Freiburg \n b-colorAGM")


# Freiburg
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(fr_spec_overyielding_m)){lines(x=x, fr_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(fr_spec_overyielding_m$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots-Freiburg \n m-colorAGM")



# All-b
plot(x=x, all_spec_overyielding_b[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(all_spec_overyielding_b)){lines(x=x, all_spec_overyielding_b[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(all_spec_overyielding_b$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots&sites \n b-colorAGM")

# All-m
plot(x=x, all_spec_overyielding_m[1,..wl_cols], type="n", lwd=2, bty="l", ylim=c(-0.6,1), 
     xlab="Wavelength (nm)", ylab="", las=1)
axis(1, at = seq(400, 2200, by = 100), labels=F, tck=-0.01)
title(ylab="NBE canopy transmittance", line=3)
for(i in 1:nrow(all_spec_overyielding_m)){lines(x=x, all_spec_overyielding_m[i,..wl_cols], col=paste0(color_palette2,90)[c(1,3,2)][as.factor(all_spec_overyielding_m$AngioGymnoMix)][i])}
rect(xleft = 1350, xright = 1452, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 1800, xright = 1983, ybottom = -1, ytop = 2000, col="white", border=NA)
rect(xleft = 2200, xright = 2505, ybottom = -1, ytop = 2000, col="white", border=NA)
abline(h=0,lty=2, col="gray45")
box(bty="l")
title(main="Full range NBE transmittance below canopy across all plots&sites \n m-colorAGM")

# dev.off()





