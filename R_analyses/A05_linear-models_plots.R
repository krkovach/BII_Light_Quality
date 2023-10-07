#_______________________________________________________________________________
# Load packages
library(data.table)
library(stringr)
require(nlme)
require(MuMIn)
require(glmm.hp) # hierarchical partitioning contributions to marginal R2
require(emmeans)

#_______________________________________________________________________________
# Load data:
file_dir <- "/Users/30062322/Library/CloudStorage/GoogleDrive-will3972@umn.edu/Shared drives/Integration Institute/Themes & Bridging Goals/5. Light quality/Data processing"

all_spec <- fread(file.path(file_dir, "Transmittance/All-Sites_ASD-SVC_mean-transmittance_20230919.txt")) 
all_spec_overyielding <- fread(file.path(file_dir, "Transmittance/All-Sites_ASD-SVC_NetBiodiversityEffect_20230919.txt"))
comm_info <- read.csv(file.path(file_dir, "Metadata/All-Sites_plot-composition-diversity-lookup.csv"))
leaf_transmittance_cwm <- read.csv(file.path(file_dir, "Integrating Sphere/Integrating sphere/All-Sites_CWM_LeafTransmittance_PerPlot-MonoTrans-AllCanopyHts-BAweighted.csv")) # these are CWMs using monoculture transmittance values (exc ACNE at FAB1 and BEPA at Freiburg, where infilled) and weighted by relative BA (which correlates with LAI)

#_______________________________________________________________________________
# Setting up data frames 

# adding plot codes 
all_spec[which(all_spec$Site=="Cloquet"),]$Plot <- with(all_spec[which(all_spec$Site=="Cloquet"),], paste(substring(BlMix,1,1),Plot, sep="_"))
all_spec$SitePlot <- with(all_spec, paste(Site, Plot, sep="_"))

leaf_transmittance_cwm$SitePlot <- with(leaf_transmittance_cwm, paste(Site_names, Plot, sep="_"))
leaf_transmittance_cwm2 <- leaf_transmittance_cwm
colnames(leaf_transmittance_cwm2)[4:2155] <- paste("Leaf",colnames(leaf_transmittance_cwm2)[4:2155],sep="_") # adding "Leaf_" in front of all wavelength names

comm_info$SitePlot <- with(comm_info, paste(Site, Plot, sep="_"))
comm_info2 <- aggregate(comm_info[,c(7:14)], by=with(comm_info,list(Site=Site, SitePlot=SitePlot, SR=SR, AngioGymnoMix=AngioGymnoMix)), FUN="mean")

# merging canopy transmittance with leaf transmittance
all_spec2_ <- merge(all_spec, leaf_transmittance_cwm2[,which(colnames(leaf_transmittance_cwm2)%in% 
                                                               c("SitePlot", "Leaf_X360", "Leaf_X390", "Leaf_X447", 
                                                                 "Leaf_X440", "Leaf_X660", "Leaf_X730", "Leaf_X865", 
                                                                 "Leaf_X1610", "Leaf_RtoFR"))], by="SitePlot", all.x=T, all.y=F)

# merging canopy transmittance with plot info
all_spec2_info <- merge(all_spec2_, comm_info2[,c(2:12)], by="SitePlot", all.x=T, all.y=F)

# removing FAB2
all_spec2 <- all_spec2_info[which(all_spec2_info$Site!="FAB2"),]


# Adding red to far-red ratio and red to SWIR
all_spec2$RtoFR <- all_spec2$X660/all_spec2$X730
all_spec2$RtoSWIR <- all_spec2$X660/all_spec2$X1610
all_spec2$Diff865.660 <- all_spec2$X865-all_spec2$X660 # "Indirect transmittance" from Hovi and Rautiainen 2020
all_spec2$CanopyInterception660 <- 1-all_spec2$X660 # "Canopy interception" from Hovi and Rautiainen 2020
all_spec2$CanopyInterception400 <- 1-all_spec2$X400 # this is the wavelength we used for minimum leaf transmittance

# Removing Freiburg plot 72 -- don't have LAI data for this plot
all_spec2 <- all_spec2[which(all_spec2$SitePlot!="Freiburg_72"),]
all_spec_overyielding <- all_spec_overyielding[which(all_spec_overyielding$ID!="72_BEPE-PISY_b"),]
all_spec_overyielding <- all_spec_overyielding[which(all_spec_overyielding$ID!="72_BEPE-PISY_m"),]

# splitting by site
cl_spec <- all_spec2[which(all_spec2$Site=="Cloquet"),]
f1_spec <- all_spec2[which(all_spec2$Site=="FAB1"),]
fr_spec <- all_spec2[which(all_spec2$Site=="Freiburg"),]

# ensuring factors set as factors (might not b necessary)
all_spec2$AngioGymnoMix <- as.factor(all_spec2$AngioGymnoMix)
all_spec_overyielding$AngioGymnoMix <- as.factor(all_spec_overyielding$AngioGymnoMix)


#_______________________________________________________________________________
# RESPONSE VARIABLES:
# Wavelengths: 360, 390, 447, 440*, 660*, 730*, 865, 1610 
# Ratio: RtoFR*
# Consider this list. Definitely include those with asterisks.
# Ordination axes don't immediately look very meaningful?

# PREDICTOR VARIABLES:
# LAI: LAI_alllloc_corrected (this is average of all locations, conifer species are multiplied by clumping correction factor)
# Composition: AngioGymnoMix (this is whether plots are composed of angiosperms, gymnosperms or mixtures of the two)
# Leaf transmittance: community weighted means of leaf transmittance at relevant wavelengths

#_______________________________________________________________________________
# PLOTTING MODEL OUTPUTS:
#_______________________________________________________________________________

color_palette2 <- c("#164673","#2A9D8F","#E9C46A") # "#264653"
color_palette2_trans <- c("#16467395","#2A9D8F95","#E9C46A95") # "#264653"
color_palette2_trans2 <- c("#16467330","#2A9D8F30","#E9C46A30") # "#264653"

plot_lai.trans.agm <- function(m1 = uv_m1, 
                               y_var = "X360", 
                               x_var = "LAI_allloc_corrected",
                               y_lab = "Canopy transmittance at 360 nm",
                               x_lab = "Leaf area index",
                               dat= data.frame(all_spec2),
                               to_log="y"){
  
  par(mar=c(5, 6, 4, 4) + 0.1, mgp=c(2,1,0))
  options(scipen=5)
  
  newdat <- expand.grid(AngioGymnoMix=c("A","G","M"), x=seq(0,20,0.1))
  colnames(newdat)[2] <- x_var
  newdat$pred <- predict(m1, newdat, level=0)
  Designmat <- model.matrix(formula(m1)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(m1) %*% t(Designmat))
  newdat$SE <- sqrt(predvar)
  newdat$SE2 <- sqrt(predvar + m1$sigma^2)
  newdat$CI_hi <- newdat$pred+1.96*newdat$SE
  newdat$CI_lo <- newdat$pred-1.96*newdat$SE
  
  y <- dat[,which(colnames(dat)==y_var)]
  x <- dat[,which(colnames(dat)==x_var)]
  plot(y~x, 
       xlab=x_lab, ylab="",
       axes=F, las=1, log=to_log, ylim=c(0.0002,1), type="n")
  
  newdat_angio <- with(newdat, newdat[which(newdat$AngioGymnoMix=="A"),])
  newdat_gymno <- with(newdat, newdat[which(newdat$AngioGymnoMix=="G"),])
  newdat_mixture <- with(newdat, newdat[which(newdat$AngioGymnoMix=="M"),])
  
  axis(side=1, at=seq(0,16,2))
  
  if(to_log=="y"){
    axis(side=2, at=c(0.001,0.01,0.1,0.5,1), las=1)
    box(bty="l")
    title(ylab=y_lab, line=4)
    
    polygon(x=c(newdat_angio[,2],rev(newdat_angio[,2])), 
            y=10^c(newdat_angio$CI_lo, rev(newdat_angio$CI_hi)), 
            border=NA, col=color_palette2_trans2[1])
    polygon(x=c(newdat_gymno[,2],rev(newdat_gymno[,2])), 
            y=10^c(newdat_gymno$CI_lo, rev(newdat_gymno$CI_hi)), 
            border=NA, col=color_palette2_trans2[3])
    polygon(x=c(newdat_mixture[,2],rev(newdat_mixture[,2])), 
            y=10^c(newdat_mixture$CI_lo, rev(newdat_mixture$CI_hi)), 
            border=NA, col=color_palette2_trans2[2])
    lines(x=newdat_angio[,2], y=10^newdat_angio$pred, col=color_palette2_trans[1], lwd=2)
    lines(x=newdat_gymno[,2], y=10^newdat_gymno$pred, col=color_palette2_trans[3], lwd=2)
    lines(x=newdat_mixture[,2], y=10^newdat_mixture$pred, col=color_palette2_trans[2], lwd=2)
  }else{
    axis(side=2, at=c(0,0.2,0.4,0.6,0.8,1), las=1)
    box(bty="l")
    title(ylab=y_lab, line=4)
    
    polygon(x=c(newdat_angio[,2],rev(newdat_angio[,2])), 
            y=c(newdat_angio$CI_lo, rev(newdat_angio$CI_hi)), 
            border=NA, col=color_palette2_trans2[1])
    polygon(x=c(newdat_gymno[,2],rev(newdat_gymno[,2])), 
            y=c(newdat_gymno$CI_lo, rev(newdat_gymno$CI_hi)), 
            border=NA, col=color_palette2_trans2[3])
    polygon(x=c(newdat_mixture[,2],rev(newdat_mixture[,2])), 
            y=c(newdat_mixture$CI_lo, rev(newdat_mixture$CI_hi)), 
            border=NA, col=color_palette2_trans2[2])
    lines(x=newdat_angio[,2], y=newdat_angio$pred, col=color_palette2_trans[1], lwd=2)
    lines(x=newdat_gymno[,2], y=newdat_gymno$pred, col=color_palette2_trans[3], lwd=2)
    lines(x=newdat_mixture[,2], y=newdat_mixture$pred, col=color_palette2_trans[2], lwd=2)
  }
  points(y ~ x, 
         col=color_palette2_trans[c(1,3,2)][as.factor(dat$AngioGymnoMix)], 
         pch=16)
  
  
  
  
  
}

# Fitting models to plot
uv_m1 <- lme(log10(X360) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
lov390_m1 <- lme(log10(X390) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
blue_m1 <- lme(log10(X440) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
lov447_m1 <- lme(log10(X447) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
red_m1 <- lme(log10(X660) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
re_m1 <- lme(log10(X730) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
nir_m1 <- lme(log10(X865) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
swir_m1 <- lme(log10(X1610) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
rtofr_m1 <- lme(RtoFR ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
rtoswir_m1 <- lme(RtoSWIR ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
indirecttrans_inter_m1 <- lme(Diff865.660 ~ CanopyInterception660*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML") # Indirect transmittance at 865 nm (865-660 nm) vs canopy interception (1-660nm)
indirecttrans_inter2_m1 <- lme(Diff865.660 ~ CanopyInterception400*AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML") # Indirect transmittance at 865 nm (865-660 nm) vs canopy interception (1-660nm)

# pdf("/Users/30062322/Documents/Research/UMN/Manuscripts/LightQuality/Figures/Scatterplots_LAI-CanopyTrans-AGM_randsite.pdf", width=6, height=6)
plot_lai.trans.agm(m1=uv_m1, y_var="X360", y_lab="Canopy transmittance at 360 nm")
plot_lai.trans.agm(m1=lov390_m1, y_var="X390", y_lab="Canopy transmittance at 390 nm")
plot_lai.trans.agm(m1=blue_m1, y_var="X440", y_lab="Canopy transmittance at 440 nm")
plot_lai.trans.agm(m1=lov447_m1, y_var="X447", y_lab="Canopy transmittance at 447 nm")
plot_lai.trans.agm(m1=red_m1, y_var="X660", y_lab="Canopy transmittance at 660 nm")
plot_lai.trans.agm(m1=re_m1, y_var="X730", y_lab="Canopy transmittance at 730 nm")
plot_lai.trans.agm(m1=nir_m1, y_var="X865", y_lab="Canopy transmittance at 865 nm")
plot_lai.trans.agm(m1=swir_m1, y_var="X1610", y_lab="Canopy transmittance at 1610 nm")
plot_lai.trans.agm(m1=rtofr_m1, y_var="RtoFR", y_lab="Ratio of canopy transmittance (660/730 nm)", to_log="")
plot_lai.trans.agm(m1=rtofr_m1, y_var="RtoSWIR", y_lab="Ratio of canopy transmittance (660/1610 nm)", to_log="")
plot_lai.trans.agm(m1=indirecttrans_inter_m1, y_var="Diff865.660", x_var = "CanopyInterception660",  x_lab = "Canopy interception (1 - transmittance at 600nm)", y_lab="Ratio of canopy transmittance (660/1610 nm)", to_log="")
axis(side=1, at=seq(0.2,1,0.2))
plot_lai.trans.agm(m1=indirecttrans_inter2_m1, y_var="Diff865.660", x_var = "CanopyInterception400",  x_lab = "Canopy interception (1 - transmittance at 400nm)", y_lab="Ratio of canopy transmittance (660/1610 nm)", to_log="")
axis(side=1, at=seq(0.2,1,0.2))

# dev.off()


# # Add R2:
# df_uv_m1 <- data.frame(plot=names(predict(uv_m1)), AngioGymnoMix=uv_m1$data$AngioGymnoMix, obs=uv_m1$data$X360, pred=predict(uv_m1, level=0))
# r.squaredGLMM(uv_m1)
# cor(log10(df_uv_m1$obs), df_uv_m1$pred, method="pearson")^2
# cor(log10(df_uv_m1[df_uv_m1$AngioGymnoMix=="A",]$obs), df_uv_m1[df_uv_m1$AngioGymnoMix=="A",]$pred, method="pearson")^2
# cor(log10(df_uv_m1[df_uv_m1$AngioGymnoMix=="G",]$obs), df_uv_m1[df_uv_m1$AngioGymnoMix=="G",]$pred, method="pearson")^2
# cor(log10(df_uv_m1[df_uv_m1$AngioGymnoMix=="M",]$obs), df_uv_m1[df_uv_m1$AngioGymnoMix=="M",]$pred, method="pearson")^2


##______________________________________________________________________________
# Plotting net biodiversity effects:

plot_nbe.lai.trans.agm <- function(m1 = nbe_uv_m1, 
                                   y_var = "X360", 
                                   x_var = "NBE_LAI_allloc_corrected",
                                   x_lab = "NBE on leaf area index",
                                   y_lab = "NBE on canopy transmittance at 360 nm",
                                   ymin = -0.6, ymax = 0.8, xmin = -4, xmax = 8,
                                   dat = data.frame(all_spec_overyielding)
                                   ){
  
  par(mar=c(5, 6, 4, 4) + 0.1, mgp=c(2,1,0))
  options(scipen=5)
  
  newdat <- expand.grid(AngioGymnoMix=c("A","G","M"), NBE_LAI_allloc_corrected=seq(-20,20,0.1))
  newdat$pred <- predict(m1, newdat, level=0)
  Designmat <- model.matrix(formula(m1)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(m1) %*% t(Designmat))
  newdat$SE <- sqrt(predvar)
  newdat$SE2 <- sqrt(predvar + m1$sigma^2)
  newdat$CI_hi <- newdat$pred+1.96*newdat$SE
  newdat$CI_lo <- newdat$pred-1.96*newdat$SE
  
  y <- dat[,which(colnames(dat)==y_var)]
  x <- dat[,which(colnames(dat)==x_var)]
  
  plot(y~x, 
       xlab=x_lab, ylab="",
       axes=F,
       ylim=c(ymin, ymax), xlim=c(xmin,xmax),
       las=1, type="n")
  
  newdat_angio <- with(newdat, newdat[which(newdat$AngioGymnoMix=="A"),])
  newdat_gymno <- with(newdat, newdat[which(newdat$AngioGymnoMix=="G"),])
  newdat_mixture <- with(newdat, newdat[which(newdat$AngioGymnoMix=="M"),])
  
  axis(side=1, at=seq(-4,8,2))
  axis(side=2, at=seq(-1,1,0.2), las=1)
  title(ylab=y_lab, line=4)
  
  polygon(x=c(newdat_angio$NBE_LAI_allloc_corrected,rev(newdat_angio$NBE_LAI_allloc_corrected)), 
          y=c(newdat_angio$CI_lo, rev(newdat_angio$CI_hi)), 
          border=NA, col=color_palette2_trans2[1])
  polygon(x=c(newdat_gymno$NBE_LAI_allloc_corrected,rev(newdat_gymno$NBE_LAI_allloc_corrected)), 
          y=c(newdat_gymno$CI_lo, rev(newdat_gymno$CI_hi)), 
          border=NA, col=color_palette2_trans2[3])
  polygon(x=c(newdat_mixture$NBE_LAI_allloc_corrected,rev(newdat_mixture$NBE_LAI_allloc_corrected)), 
          y=c(newdat_mixture$CI_lo, rev(newdat_mixture$CI_hi)), 
          border=NA, col=color_palette2_trans2[2])
  with(newdat_angio,lines(x=NBE_LAI_allloc_corrected, y=pred, col=color_palette2_trans[1], lwd=2))
  with(newdat_gymno,lines(x=NBE_LAI_allloc_corrected, y=pred, col=color_palette2_trans[3], lwd=2))
  with(newdat_mixture,lines(x=NBE_LAI_allloc_corrected, y=pred, col=color_palette2_trans[2], lwd=2))
  
  abline(h=0,lty=2,col="grey")
  abline(v=0,lty=2,col="grey")
  
  box(bty="l")
  
  points(y ~ x, 
         col=color_palette2_trans[c(1,3,2)][as.factor(dat$AngioGymnoMix)], 
         pch=16)
}

nbe_uv_m1 <- lme(X360 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_lov390_m1 <- lme(X390 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_blue_m1 <- lme(X440 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_lov447_m1 <- lme(X447 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_red_m1 <- lme(X660 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_re_m1 <- lme(X730 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_nir_m1 <- lme(X865 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_swir_m1 <- lme(X1610 ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_rtofr_m1 <- lme(RtoFR ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")
nbe_rtoswir_m1 <- lme(RtoSWIR ~ NBE_LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML")


# pdf("/Users/30062322/Documents/Research/UMN/Manuscripts/LightQuality/Figures/Scatterplots_NBE_LAI-CanopyTrans-AGM_randsite.pdf", width=6, height=6)
plot_nbe.lai.trans.agm(m1=nbe_uv_m1, y_var="X360", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       x_lab = "NBE on leaf area index",
                       y_lab="NBE on canopy transmittance at 360 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_lov390_m1, y_var="X390", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 390 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_blue_m1, y_var="X440", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 440 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_lov447_m1, y_var="X447", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 447 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_red_m1, y_var="X660", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 660 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_re_m1, y_var="X730", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 730 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_nir_m1, y_var="X865", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 865 nm",
                       dat=data.frame(all_spec_overyielding))
plot_nbe.lai.trans.agm(m1=nbe_swir_m1, y_var="X1610", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE on canopy transmittance at 1610 nm",
                       dat=data.frame(all_spec_overyielding), ymin=-1, ymax=1)
plot_nbe.lai.trans.agm(m1=nbe_rtofr_m1, y_var="RtoFR", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE of ratio of canopy transmittance (660/730 nm)",
                       dat=data.frame(all_spec_overyielding), ymin=-1, ymax=1)
plot_nbe.lai.trans.agm(m1=nbe_rtoswir_m1, y_var="RtoFR", 
                       x_var = "NBE_LAI_allloc_corrected", 
                       y_lab="NBE of ratio of canopy transmittance (660/1610 nm)",
                       dat=data.frame(all_spec_overyielding), ymin=-1, ymax=1)
# dev.off()

# pdf("/Users/30062322/Documents/Research/UMN/Manuscripts/LightQuality/Figures/Scatterplot_legend.pdf", width=6, height=6)
plot(1,1, axes=F, type="n", xlab="", ylab="")
legend("topright", pch=16, col=color_palette2[c(1,3,2)], legend=c("Angiosperms", "Gymnosperms", "Angio-gymno mixtures"), bty="n")
legend("bottomright", pch=16, col=color_palette2_trans[c(1,3,2)], legend=c("Angiosperms", "Gymnosperms", "Angio-gymno mixtures"), bty="n")
# dev.off()

# A_E8_PIGL-PIAB_m --- this has very high NBE transmittance value... error somewhere here?

wls <- as.numeric(substring(colnames(all_spec2[,18:1590]),2,20))
plot(x=wls, y=all_spec2[which(all_spec2$ID=="A_E8_PIGL-PIAB_m"),18:1590], type="l", ylab="transmittance", ylim=c(0,1.2))
plot(x=wls, y=all_spec_overyielding[which(all_spec_overyielding$ID=="A_E8_PIGL-PIAB_m"),18:1590], type="l", ylab="transmittance", ylim=c(-1,1.2))







