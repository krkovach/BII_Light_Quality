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
# FITTING MODELS:
#_______________________________________________________________________________
# UV: 360
boxplot(all_spec2$X360) # log
plot(X360~LAI_allloc_corrected, data=all_spec2)
plot(X360~LAI_allloc_corrected, log="y", data=all_spec2, col=as.factor(AngioGymnoMix), pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X360~LAI_allloc_corrected, log="xy", data=all_spec2, col=as.factor(AngioGymnoMix), pch=16)

# Fit and compare models with ML
uv_mfull_ml <- lme(log10(X360) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="ML")
uv_m1_ml <- lme(log10(X360) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
uv_m2_ml <- lme(log10(X360) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="ML") # singularity warning
uv_m3_ml <- lme(log10(X360) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML") # singularity warning
uv_m4_ml <- lme(log10(X360) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
uv_m5_ml <- lme(log10(X360) ~ LAI_allloc_corrected + Leaf_X360, random=~1|Site, data=all_spec2, method="ML")
uv_null_ml <- lme(log10(X360) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(uv_mfull_ml, uv_m1_ml) # dropping LOP has no effect
anova(uv_mfull_ml, uv_m2_ml) # dropping AGM-LAI interaction is very significant
anova(uv_m2_ml, uv_m3_ml)
anova(uv_m3_ml, uv_m4_ml)
anova(uv_m4_ml, uv_m5_ml) # Is AGM soaking up variation that would otherwise be explained by LOP: without having AGM in model, does adding Leaf_X360 do anything? -- very little
anova(uv_m4_ml, uv_null_ml)

r.squaredGLMM(uv_m4_ml)
r.squaredGLMM(uv_m5_ml)

# Models with REML: 
# Full model
uv_mfull <- lme(log10(X360) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
uv_m1 <- lme(log10(X360) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
uv_m2 <- lme(log10(X360) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
uv_lmfull <- lm(log10(X360) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X360, data=all_spec2)
anova(uv_mfull, uv_lmfull) # doesn't make a difference here (could use ranova() from lmerTest package)
anova(uv_mfull_ml, uv_lmfull) # technically should use ML (though automatically refitted when don't)

summary(uv_mfull)
summary(uv_m1)
r.squaredGLMM(uv_mfull) # with all terms
r.squaredGLMM(uv_m1) # with LOP dropped -- R2m drop is extremely small 
r.squaredGLMM(uv_mfull)-r.squaredGLMM(uv_m1)

predict(uv_mfull, level=0) # population model
fitted(uv_mfull, level=0) # population model (another way of getting this) -- matches manual estimates of fitted values (from just fixed effects)
fitted(uv_mfull, level=1)

# Calculating R2 values for terms from the full model:
df_uv_mfull <- data.frame(plot=names(predict(uv_mfull)), AngioGymnoMix=uv_mfull$data$AngioGymnoMix, obs=uv_mfull$data$X360, pred=predict(uv_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
# We are modelling log transmittance values, so I presume we should report model with logged transmittance values:
cor(log10(df_uv_mfull$obs), df_uv_mfull$pred, method="pearson")^2
cor(log10(df_uv_mfull[df_uv_mfull$AngioGymnoMix=="A",]$obs), df_uv_mfull[df_uv_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_uv_mfull[df_uv_mfull$AngioGymnoMix=="G",]$obs), df_uv_mfull[df_uv_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_uv_mfull[df_uv_mfull$AngioGymnoMix=="M",]$obs), df_uv_mfull[df_uv_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2


# Hierarchical partitioning: to assess independent and joint explanatory power of each variable
# Not sure if/how to add interaction terms? (i.e., LAI:AngioGymnoMix)
# Assessing without interactions -- I am not sure if this makes sense:
uv_m2 <- lme(log10(X360) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="REML")
glmm.hp(uv_m2, commonality=T)
glmm.hp(uv_m2)
glmm.hp(uv_m5_ml)

#_______________________________________________________________________________
# blue: 440
boxplot(all_spec2$X440) # log
plot(X440~LAI_allloc_corrected, data=all_spec2)
plot(X440~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X440~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
blue_mfull_ml <- lme(log10(X440) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X440, random=~1|Site, data=all_spec2, method="ML")
blue_m1_ml <- lme(log10(X440) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
blue_m2_ml <- lme(log10(X440) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X440, random=~1|Site, data=all_spec2, method="ML")
blue_m3_ml <- lme(log10(X440) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
blue_m4_ml <- lme(log10(X440) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
blue_null_ml <- lme(log10(X440) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(blue_mfull_ml, blue_m1_ml) # dropping LOP has no effect
anova(blue_mfull_ml, blue_m2_ml) # dropping AGM-LAI interaction is very significant
anova(blue_m2_ml, blue_m3_ml)
anova(blue_m3_ml, blue_m4_ml)
anova(blue_m4_ml, blue_null_ml)

# Models with REML: 
# Full model
blue_mfull <- lme(log10(X440) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X440, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
blue_m1 <- lme(log10(X440) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
blue_m2 <- lme(log10(X440) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X440, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
blue_lmfull <- lm(log10(X440) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X440, data=all_spec2)
anova(blue_mfull, blue_lmfull)

summary(blue_mfull)
summary(blue_m1)
r.squaredGLMM(blue_mfull) # with all terms
r.squaredGLMM(blue_m1) 
r.squaredGLMM(blue_mfull)[1]-r.squaredGLMM(blue_m1)[1]# minimal r2 drop by dropping LOP

glmm.hp(blue_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

# Calculating R2 values for terms from the full model:
r.squaredGLMM(blue_mfull)
df_blue_mfull <- data.frame(plot=names(predict(blue_mfull)), AngioGymnoMix=blue_mfull$data$AngioGymnoMix, obs=blue_mfull$data$X440, pred=predict(blue_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
cor(log10(df_blue_mfull$obs), df_blue_mfull$pred, method="pearson")^2
cor(log10(df_blue_mfull[df_blue_mfull$AngioGymnoMix=="A",]$obs), df_blue_mfull[df_blue_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_blue_mfull[df_blue_mfull$AngioGymnoMix=="G",]$obs), df_blue_mfull[df_blue_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_blue_mfull[df_blue_mfull$AngioGymnoMix=="M",]$obs), df_blue_mfull[df_blue_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2

#_______________________________________________________________________________
# red: 660
boxplot(all_spec2$X660) # log
plot(X660~LAI_allloc_corrected, data=all_spec2)
plot(X660~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X660~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
red_mfull_ml <- lme(log10(X660) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X660, random=~1|Site, data=all_spec2, method="ML")
red_m1_ml <- lme(log10(X660) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
red_m2_ml <- lme(log10(X660) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X660, random=~1|Site, data=all_spec2, method="ML")
red_m3_ml <- lme(log10(X660) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
red_m4_ml <- lme(log10(X660) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
red_null_ml <- lme(log10(X660) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(red_mfull_ml, red_m1_ml) # dropping LOP has no effect
anova(red_mfull_ml, red_m2_ml) # dropping AGM-LAI interaction is very significant
anova(red_m2_ml, red_m3_ml)
anova(red_m3_ml, red_m4_ml)
anova(red_m4_ml, red_null_ml)

# Models with REML: 
# Full model
red_mfull <- lme(log10(X660) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X660, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
red_m1 <- lme(log10(X660) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
red_m2 <- lme(log10(X660) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X660, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
red_lmfull <- lm(log10(X660) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X660, data=all_spec2)
anova(red_mfull, red_lmfull)

summary(red_mfull)
summary(red_m1)
r.squaredGLMM(red_mfull) # with all terms
r.squaredGLMM(red_m1) 
r.squaredGLMM(red_mfull)[1]-r.squaredGLMM(red_m1)[1]# minimal r2 drop by dropping LOP

glmm.hp(red_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

# Calculating R2 values for terms from the full model:
r.squaredGLMM(red_mfull)
df_red_mfull <- data.frame(plot=names(predict(red_mfull)), AngioGymnoMix=red_mfull$data$AngioGymnoMix, obs=red_mfull$data$X660, pred=predict(red_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
cor(log10(df_red_mfull$obs), df_red_mfull$pred, method="pearson")^2
cor(log10(df_red_mfull[df_red_mfull$AngioGymnoMix=="A",]$obs), df_red_mfull[df_red_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_red_mfull[df_red_mfull$AngioGymnoMix=="G",]$obs), df_red_mfull[df_red_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_red_mfull[df_red_mfull$AngioGymnoMix=="M",]$obs), df_red_mfull[df_red_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2


#_______________________________________________________________________________
# fr: 730
boxplot(all_spec2$X730) # log
plot(X730~LAI_allloc_corrected, data=all_spec2)
plot(X730~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X730~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
fr_mfull_ml <- lme(log10(X730) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X730, random=~1|Site, data=all_spec2, method="ML")
fr_m1_ml <- lme(log10(X730) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
fr_m2_ml <- lme(log10(X730) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X730, random=~1|Site, data=all_spec2, method="ML")
fr_m3_ml <- lme(log10(X730) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
fr_m4_ml <- lme(log10(X730) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
fr_null_ml <- lme(log10(X730) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(fr_mfull_ml, fr_m1_ml) # dropping LOP has no effect
anova(fr_mfull_ml, fr_m2_ml) # dropping AGM-LAI interaction is very significant
anova(fr_m2_ml, fr_m3_ml)
anova(fr_m3_ml, fr_m4_ml)
anova(fr_m4_ml, fr_null_ml)

# Models with REML: 
# Full model
fr_mfull <- lme(log10(X730) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X730, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
fr_m1 <- lme(log10(X730) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
fr_m2 <- lme(log10(X730) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X730, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
fr_lmfull <- lm(log10(X730) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X730, data=all_spec2)
anova(fr_mfull, fr_lmfull)

summary(fr_mfull)
summary(fr_m1)
r.squaredGLMM(fr_mfull) # with all terms
r.squaredGLMM(fr_m1) 
r.squaredGLMM(fr_mfull)[1]-r.squaredGLMM(fr_m1)[1]# minimal r2 drop by dropping LOP

glmm.hp(fr_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

# Calculating R2 values for terms from the full model:
r.squaredGLMM(fr_mfull)
df_fr_mfull <- data.frame(plot=names(predict(fr_mfull)), AngioGymnoMix=fr_mfull$data$AngioGymnoMix, obs=fr_mfull$data$X730, pred=predict(fr_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
cor(log10(df_fr_mfull$obs), df_fr_mfull$pred, method="pearson")^2
cor(log10(df_fr_mfull[df_fr_mfull$AngioGymnoMix=="A",]$obs), df_fr_mfull[df_fr_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_fr_mfull[df_fr_mfull$AngioGymnoMix=="G",]$obs), df_fr_mfull[df_fr_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_fr_mfull[df_fr_mfull$AngioGymnoMix=="M",]$obs), df_fr_mfull[df_fr_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2



#_______________________________________________________________________________
# nir: 865
boxplot(all_spec2$X865) # log
plot(X865~LAI_allloc_corrected, data=all_spec2)
plot(X865~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X865~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
nir_mfull_ml <- lme(log10(X865) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X865, random=~1|Site, data=all_spec2, method="ML")
nir_m1_ml <- lme(log10(X865) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
nir_m2_ml <- lme(log10(X865) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X865, random=~1|Site, data=all_spec2, method="ML")
nir_m3_ml <- lme(log10(X865) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
nir_m4_ml <- lme(log10(X865) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
nir_null_ml <- lme(log10(X865) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(nir_mfull_ml, nir_m1_ml) # dropping LOP has no effect
anova(nir_mfull_ml, nir_m2_ml) # dropping AGM-LAI interaction is very significant
anova(nir_m2_ml, nir_m3_ml)
anova(nir_m3_ml, nir_m4_ml)
anova(nir_m4_ml, nir_null_ml)

# Models with REML: 
# Full model
nir_mfull <- lme(log10(X865) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X865, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
nir_m1 <- lme(log10(X865) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
nir_m2 <- lme(log10(X865) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X865, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
nir_lmfull <- lm(log10(X865) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X865, data=all_spec2)
anova(nir_mfull, nir_lmfull)

summary(nir_mfull)
summary(nir_m1)
r.squaredGLMM(nir_mfull) # with all terms
r.squaredGLMM(nir_m1) 
r.squaredGLMM(nir_mfull)[1]-r.squaredGLMM(nir_m1)[1]# minimal r2 drop by dropping LOP

glmm.hp(nir_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

# Calculating R2 values for terms nirom the full model:
r.squaredGLMM(nir_mfull)
df_nir_mfull <- data.frame(plot=names(predict(nir_mfull)), AngioGymnoMix=nir_mfull$data$AngioGymnoMix, obs=nir_mfull$data$X865, pred=predict(nir_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
cor(log10(df_nir_mfull$obs), df_nir_mfull$pred, method="pearson")^2
cor(log10(df_nir_mfull[df_nir_mfull$AngioGymnoMix=="A",]$obs), df_nir_mfull[df_nir_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_nir_mfull[df_nir_mfull$AngioGymnoMix=="G",]$obs), df_nir_mfull[df_nir_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_nir_mfull[df_nir_mfull$AngioGymnoMix=="M",]$obs), df_nir_mfull[df_nir_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2


#_______________________________________________________________________________
# swir: 1610
boxplot(all_spec2$X1610) # log
plot(X1610~LAI_allloc_corrected, data=all_spec2)
plot(X1610~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X1610~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
swir_mfull_ml <- lme(log10(X1610) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X1610, random=~1|Site, data=all_spec2, method="ML")
swir_m1_ml <- lme(log10(X1610) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
swir_m2_ml <- lme(log10(X1610) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X1610, random=~1|Site, data=all_spec2, method="ML")
swir_m3_ml <- lme(log10(X1610) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
swir_m4_ml <- lme(log10(X1610) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
swir_null_ml <- lme(log10(X1610) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(swir_mfull_ml, swir_m1_ml) # dropping LOP has no effect
anova(swir_mfull_ml, swir_m2_ml) # dropping AGM-LAI interaction is very significant
anova(swir_m2_ml, swir_m3_ml)
anova(swir_m3_ml, swir_m4_ml)
anova(swir_m4_ml, swir_null_ml)

# Models with REML: 
# Full model
swir_mfull <- lme(log10(X1610) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X1610, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
swir_m1 <- lme(log10(X1610) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
swir_m2 <- lme(log10(X1610) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X1610, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
swir_lmfull <- lm(log10(X1610) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X1610, data=all_spec2)
anova(swir_mfull, swir_lmfull)

summary(swir_mfull)
summary(swir_m1)
r.squaredGLMM(swir_mfull) # with all terms
r.squaredGLMM(swir_m1) 
r.squaredGLMM(swir_mfull)[1]-r.squaredGLMM(swir_m1)[1]# minimal r2 drop by dropping LOP

glmm.hp(swir_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

# Calculating R2 values for terms from the full model:
r.squaredGLMM(swir_mfull)
df_swir_mfull <- data.frame(plot=names(predict(swir_mfull)), AngioGymnoMix=swir_mfull$data$AngioGymnoMix, obs=swir_mfull$data$X1610, pred=predict(swir_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
cor(log10(df_swir_mfull$obs), df_swir_mfull$pred, method="pearson")^2
cor(log10(df_swir_mfull[df_swir_mfull$AngioGymnoMix=="A",]$obs), df_swir_mfull[df_swir_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_swir_mfull[df_swir_mfull$AngioGymnoMix=="G",]$obs), df_swir_mfull[df_swir_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_swir_mfull[df_swir_mfull$AngioGymnoMix=="M",]$obs), df_swir_mfull[df_swir_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2


#_______________________________________________________________________________
# LOV: 390
boxplot(all_spec2$X390) # log
plot(X390~LAI_allloc_corrected, data=all_spec2)
plot(X390~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X390~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
lov390_mfull_ml <- lme(log10(X390) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X390, random=~1|Site, data=all_spec2, method="ML")
lov390_m1_ml <- lme(log10(X390) ~ LAI_allloc_corrected*AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
lov390_m2_ml <- lme(log10(X390) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X390, random=~1|Site, data=all_spec2, method="ML")
lov390_m3_ml <- lme(log10(X390) ~ LAI_allloc_corrected + AngioGymnoMix, random=~1|Site, data=all_spec2, method="ML")
lov390_m4_ml <- lme(log10(X390) ~ LAI_allloc_corrected, random=~1|Site, data=all_spec2, method="ML")
lov390_null_ml <- lme(log10(X390) ~ 1, random=~1|Site, data=all_spec2, method="ML")
anova(lov390_mfull_ml, lov390_m1_ml) # dropping LOP has no effect
anova(lov390_mfull_ml, lov390_m2_ml) # dropping AGM-LAI interaction is very significant
anova(lov390_m2_ml, lov390_m3_ml)
anova(lov390_m3_ml, lov390_m4_ml)
anova(lov390_m4_ml, lov390_null_ml)

# Models with REML: 
# Full model
lov390_mfull <- lme(log10(X390) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X390, random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
lov390_m1 <- lme(log10(X390) ~ LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
lov390_m2 <- lme(log10(X390) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X390, random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
lov390_lmfull <- lm(log10(X390) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X390, data=all_spec2)
anova(lov390_mfull, lov390_lmfull) # doesn't make a difference here (could use ranova() from lmerTest package)
anova(lov390_mfull_ml, lov390_lmfull) # technically should use ML (though automatically refitted when don't)

summary(lov390_mfull)
summary(lov390_m1)
r.squaredGLMM(lov390_mfull) # with all terms
r.squaredGLMM(lov390_m1) # with LOP dropped -- R2m drop is extremely small 
glmm.hp(lov390_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

library(lattice)
randoms <- ranef(lov390_mfull)
dotplot(randoms)
fixef(lov390_mfull)

# Calculating R2 values for terms from the full model:
r.squaredGLMM(lov390_mfull)
df_lov390_mfull <- data.frame(plot=names(predict(lov390_mfull)), AngioGymnoMix=lov390_mfull$data$AngioGymnoMix, obs=lov390_mfull$data$X390, pred=predict(lov390_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
# We are modelling log transmittance values, so I presume we should report model with logged transmittance values:
cor(log10(df_lov390_mfull$obs), df_lov390_mfull$pred, method="pearson")^2
cor(log10(df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="A",]$obs), df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="G",]$obs), df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="M",]$obs), df_lov390_mfull[df_lov390_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2

#_______________________________________________________________________________
# LOV: 447
boxplot(all_spec2$X447) # log
plot(X447~LAI_allloc_corrected, data=all_spec2)
plot(X447~LAI_allloc_corrected, log="y", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16) # check residuals to see whether logging both x and y helps or just y
plot(X447~LAI_allloc_corrected, log="xy", data=all_spec2, col=color_palette2_trans[as.factor(AngioGymnoMix)], pch=16)

# Fit and compare models with ML
lov447_mfull_ml <- lme(log10(X447) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X447 , random=~1|Site, data=all_spec2, method="ML")
lov447_m1_ml <- lme(log10(X447) ~ LAI_allloc_corrected*AngioGymnoMix , random=~1|Site, data=all_spec2, method="ML") # singularity warning
lov447_m2_ml <- lme(log10(X447) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X447 , random=~1|Site, data=all_spec2, method="ML") # singularity warning
lov447_m3_ml <- lme(log10(X447) ~ LAI_allloc_corrected + AngioGymnoMix , random=~1|Site, data=all_spec2, method="ML") # singularity warning
lov447_m4_ml <- lme(log10(X447) ~ LAI_allloc_corrected , random=~1|Site, data=all_spec2, method="ML")
lov447_m5_ml <- lme(log10(X447) ~ 1 , random=~1|Site, data=all_spec2, method="ML")
anova(lov447_mfull_ml, lov447_m1_ml) # dropping LOP has no effect
anova(lov447_mfull_ml, lov447_m2_ml) # dropping AGM-LAI interaction is very significant
anova(lov447_m2_ml, lov447_m3_ml)
anova(lov447_m3_ml, lov447_m4_ml)
anova(lov447_m4_ml, lov447_m5_ml)

# Models with REML: 
# Full model
lov447_mfull <- lme(log10(X447) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X447 , random=~1|Site, data=all_spec2, method="REML")
# Model without LOP
lov447_m1 <- lme(log10(X447) ~ LAI_allloc_corrected * AngioGymnoMix , random=~1|Site, data=all_spec2, method="REML")
# Model without interactions (needed for hierarchical partitioning -- not sure how to code interaction as a separate variable?)
lov447_m2 <- lme(log10(X447) ~ LAI_allloc_corrected + AngioGymnoMix + Leaf_X447 , random=~1|Site, data=all_spec2, method="REML")
# Checking effect of including random effect for site 
lov447_lmfull <- lm(log10(X447) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X447, data=all_spec2)
anova(lov447_mfull, lov447_lmfull) # doesn't make a difference here (could use ranova() from lmeTest package)
anova(lov447_mfull_ml, lov447_lmfull) # technically should use ML (though automatically refitted when don't)

summary(lov447_mfull)
summary(lov447_m1)
r.squaredGLMM(lov447_mfull) # with all terms
r.squaredGLMM(lov447_m1) # with LOP dropped -- R2m drop is extremely small 
glmm.hp(lov447_m2) # not sure if this is meaningful without having LAI:AngioGymnoMix accounted for?

library(lattice)
randoms <- ranef(lov447_mfull)
dotplot(randoms)
fixef(lov447_mfull)
intervals(lov447_mfull)

coef(summary(lov447_mfull))

# Calculating R2 values for terms from the full model:
r.squaredGLMM(lov447_mfull)
df_lov447_mfull <- data.frame(plot=names(predict(lov447_mfull)), AngioGymnoMix=lov447_mfull$data$AngioGymnoMix, obs=lov447_mfull$data$X447, pred=predict(uv_mfull, level=0)) # AngioGymnoMix=rep(c("A","G","M"),200) # not quite sure how to represent this
# We are modelling log transmittance values, so I presume we should report model with logged transmittance values:
cor(log10(df_lov447_mfull$obs), df_lov447_mfull$pred, method="pearson")^2
cor(log10(df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="A",]$obs), df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="A",]$pred, method="pearson")^2
cor(log10(df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="G",]$obs), df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="G",]$pred, method="pearson")^2
cor(log10(df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="M",]$obs), df_lov447_mfull[df_lov447_mfull$AngioGymnoMix=="M",]$pred, method="pearson")^2


#_______________________________________________________________________________
# NET BIODIVERSITY EFFECTS:
summary(nbe_uv_mfull <- lme(X360 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_blue_mfull <- lme(X440 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_red_mfull <- lme(X660 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_re_mfull <- lme(X730 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_nir_mfull <- lme(X865 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_swir_mfull <- lme(X1610 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_lov390_mfull <- lme(X390 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_lov447_mfull <- lme(X447 ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))
summary(nbe_rtofr_mfull <- lme(RtoFR ~ NBE_LAI_allloc_corrected * AngioGymnoMix, random=~1|Site, data=all_spec_overyielding, method="REML"))

r.squaredGLMM(nbe_uv_mfull)
r.squaredGLMM(nbe_blue_mfull)
r.squaredGLMM(nbe_red_mfull)
r.squaredGLMM(nbe_re_mfull)
r.squaredGLMM(nbe_nir_mfull)
r.squaredGLMM(nbe_swir_mfull)
r.squaredGLMM(nbe_lov390_mfull)
r.squaredGLMM(nbe_lov447_mfull)
r.squaredGLMM(nbe_rtofr_mfull)





#_______________________________________________________________________________
# Pulling out key models and stats:
#Wavelengths: 360, 390, 447, 440*, 660*, 730*, 865, 1610 
summary(uv_mfull <- lme(log10(X360) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X360, random=~1|Site, data=all_spec2, method="REML"))
summary(blue_mfull <- lme(log10(X440) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X440, random=~1|Site, data=all_spec2, method="REML"))
summary(red_mfull <- lme(log10(X660) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X660, random=~1|Site, data=all_spec2, method="REML"))
summary(re_mfull <- lme(log10(X730) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X730, random=~1|Site, data=all_spec2, method="REML"))
summary(nir_mfull <- lme(log10(X865) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X865, random=~1|Site, data=all_spec2, method="REML"))
summary(swir_mfull <- lme(log10(X1610) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X1610, random=~1|Site, data=all_spec2, method="REML"))
summary(lov390_mfull <- lme(log10(X390) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X390, random=~1|Site, data=all_spec2, method="REML"))
summary(lov447_mfull <- lme(log10(X447) ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_X447, random=~1|Site, data=all_spec2, method="REML"))
summary(rtofr_mfull <- lme(RtoFR ~ LAI_allloc_corrected * AngioGymnoMix + Leaf_RtoFR, random=~1|Site, data=all_spec2, method="REML"))

r.squaredGLMM(uv_mfull)
r.squaredGLMM(blue_mfull)
r.squaredGLMM(red_mfull)
r.squaredGLMM(re_mfull)
r.squaredGLMM(nir_mfull)
r.squaredGLMM(swir_mfull)
r.squaredGLMM(lov390_mfull)
r.squaredGLMM(lov447_mfull)
r.squaredGLMM(rtofr_mfull)



###_____________________________________________________________________________
### TO DO:
### FIGURE OUT HOW TO ASSESS CONTRIBUTION OF ADDING LOP (comparing model with and without this term? and looking at increase in marginal R2?; or if wanted to run as a second model looking at residuals, how to treat AGM [should AGM and its interactions be included in the first model, second model, none, both?]?)
## Look at LOP differences between t, m, b of canopy














