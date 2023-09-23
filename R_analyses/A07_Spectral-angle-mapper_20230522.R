# Applying spectral angle mapper (SAM) to canopy transmittance spectra
# Performing ordinations (NMDS, PCoA) on the resulting distance matrix

#_______________________________________________________________________________
# Load packages
library(data.table)
library(stringr)
require(spectrolab)
require(hsdar)
require(vegan)

#_______________________________________________________________________________
# Load data

file_dir <- "/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Transmittance"
# all_spec_ <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230328.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed
all_spec_ <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230725.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed

# making some updates so that the community info easily merges:
all_spec_[which(all_spec_$Site=="Cloquet"),]$Plot <- with(all_spec_[which(all_spec_$Site=="Cloquet"),], paste(substring(BlMix,1,1), Plot, sep="_"))
all_spec_[which(all_spec_$Site=="FAB1"),]$ID <- with(all_spec_[which(all_spec_$Site=="FAB1"),], paste(Site, Plot, Mix ,CanopyHt,sep="_"))
all_spec_[which(all_spec_$Site=="FAB2"),]$ID <- with(all_spec_[which(all_spec_$Site=="FAB2"),], paste(Site, Plot, Mix ,CanopyHt,sep="_"))

comm_info <- read.csv("/Users/30062322/Documents/Research/UMN/IDENT/Light quality/Data processing/Metadata/All-Sites_plot-composition-diversity-lookup.csv")

all_spec <- merge(all_spec_, comm_info[,c(1,5:8)], by="ID", all.x=T, all.y=F)

#_______________________________________________________________________________
# Setting up spectral data as speclib
# (Spectral Angle Mapper is in hsdar package, which requires data in this format)

sub_all_spec <- all_spec[which(all_spec$Site!="FAB2"),] # excluding FAB2 (some plots have bad transmittance values)

wls <- as.numeric(c(350:1350, 1450:1800, 1980:2200)) # wavelengths
specc <- as.matrix(sub_all_spec[,17:1589]) # spectral data 11:1583
ids <- with(sub_all_spec, paste(Site, Plot, Mix, CanopyHt, sep="_")) # IDs for each spectrum
si <- data.frame(sub_all_spec[,c(1:16,1590:1593)]) # other supplementary info
  
spec_all_s <- speclib(specc, wls) # setting up speclib

idSpeclib(spec_all_s) <- as.character(ids) # adding IDs for each spectrum
SI(spec_all_s) <- si # adding supplementary information

spec_b_s <- subset(spec_all_s, CanopyHt=="b") # subsetting below canopy spectra
spec_m_s <- subset(spec_all_s, CanopyHt=="m") # subsetting mid canopy spectra

# # Now for different parts of the spectrum (note that masking gets rid of all values BETWEEN those specified, not inclusive)
# Note this approach doesn't seem to work when using objects in SAM -- tf set up manually
# spec_UV <- spec_all_s
# mask(spec_UV) <- c(399,2201) 
# spec_VIS <- spec_all_s
# mask(spec_VIS) <- c(349,400)
# mask(spec_VIS) <- c(700,2201)
# spec_RE <- spec_all_s
# mask(spec_RE) <- c(349,670) 
# mask(spec_RE) <- c(780,2201) 
# spec_NIR <- spec_all_s
# mask(spec_NIR) <- c(349,700) 
# mask(spec_NIR) <- c(1100,2201) 
# spec_SWIR <- spec_all_s
# mask(spec_SWIR) <- c(349,1100) 
#------
# manually setting:
wls_uv <- as.numeric(c(350:400)) # wavelengths
col_numbers <- which(colnames(sub_all_spec)%in%c("X350","X400"))
specc_uv <- as.matrix(sub_all_spec[,col_numbers[1]:col_numbers[2]]) # spectral data
spec_uv <- speclib(specc_uv, wls_uv) # setting up speclib
idSpeclib(spec_uv) <- as.character(ids) # adding IDs for each spectrum
SI(spec_uv) <- si # adding supplementary information

wls_vis <- as.numeric(c(400:700)) # wavelengths
col_numbers <- which(colnames(sub_all_spec)%in%c("X400","X700"))
specc_vis <- as.matrix(sub_all_spec[,col_numbers[1]:col_numbers[2]]) # spectral data
spec_vis <- speclib(specc_vis, wls_vis) # setting up speclib
idSpeclib(spec_vis) <- as.character(ids) # adding IDs for each spectrum
SI(spec_vis) <- si # adding supplementary information

wls_re <- as.numeric(c(670:780)) # wavelengths
col_numbers <- which(colnames(sub_all_spec)%in%c("X670","X780"))
specc_re <- as.matrix(sub_all_spec[,col_numbers[1]:col_numbers[2]]) # spectral data
spec_re <- speclib(specc_re, wls_re) # setting up speclib
idSpeclib(spec_re) <- as.character(ids) # adding IDs for each spectrum
SI(spec_re) <- si # adding supplementary information

wls_nir <- as.numeric(c(701:1100)) # wavelengths
col_numbers <- which(colnames(sub_all_spec)%in%c("X701","X1100"))
specc_nir <- as.matrix(sub_all_spec[,col_numbers[1]:col_numbers[2]]) # spectral data
spec_nir <- speclib(specc_nir, wls_nir) # setting up speclib
idSpeclib(spec_nir) <- as.character(ids) # adding IDs for each spectrum
SI(spec_nir) <- si # adding supplementary information

wls_swir <- as.numeric(c(1100:1350, 1450:1800, 1980:2200)) # wavelengths
col_numbers <- which(colnames(sub_all_spec)%in%c("X1100","X2200"))
specc_swir <- as.matrix(sub_all_spec[,col_numbers[1]:col_numbers[2]]) # spectral data
spec_swir <- speclib(specc_swir, wls_swir) # setting up speclib
idSpeclib(spec_swir) <- as.character(ids) # adding IDs for each spectrum
SI(spec_swir) <- si # adding supplementary information

wls_biol <- as.numeric(c(360,450,600,730)) # wavelengths of biological importance (TO BE DETERMINED)
col_numbers <- which(colnames(sub_all_spec)%in%c("X360","X450","X600","X730")) 
specc_biol <- as.matrix(sub_all_spec[,..col_numbers]) # spectral data
spec_biol <- speclib(specc_biol, wls_biol) # setting up speclib
idSpeclib(spec_biol) <- as.character(ids) # adding IDs for each spectrum
SI(spec_biol) <- si # adding supplementary information

spec_b_uv <- subset(spec_uv, CanopyHt=="b") # subsetting below canopy spectra
spec_m_uv <- subset(spec_uv, CanopyHt=="m") # subsetting mid canopy spectra
spec_b_vis <- subset(spec_vis, CanopyHt=="b") 
spec_m_vis <- subset(spec_vis, CanopyHt=="m") 
spec_b_re <- subset(spec_re, CanopyHt=="b") 
spec_m_re <- subset(spec_re, CanopyHt=="m") 
spec_b_nir <- subset(spec_nir, CanopyHt=="b") 
spec_m_nir <- subset(spec_nir, CanopyHt=="m") 
spec_b_swir <- subset(spec_swir, CanopyHt=="b") 
spec_m_swir <- subset(spec_swir, CanopyHt=="m") 
spec_b_biol <- subset(spec_biol, CanopyHt=="b") 
spec_m_biol <- subset(spec_biol, CanopyHt=="m") 

spec_b_s_cl <- subset(spec_all_s, CanopyHt=="b"&Site=="Cloquet") # subsetting below canopy spectra
spec_m_s_cl <- subset(spec_all_s, CanopyHt=="m"&Site=="Cloquet") # subsetting mid canopy spectra; presumably need to add site throughout
spec_b_uv_cl <- subset(spec_uv, CanopyHt=="b"&Site=="Cloquet") # subsetting below canopy spectra
spec_m_uv_cl <- subset(spec_uv, CanopyHt=="m"&Site=="Cloquet") # subsetting mid canopy spectra
spec_b_vis_cl <- subset(spec_vis, CanopyHt=="b"&Site=="Cloquet") 
spec_m_vis_cl <- subset(spec_vis, CanopyHt=="m"&Site=="Cloquet") 
spec_b_re_cl <- subset(spec_re, CanopyHt=="b"&Site=="Cloquet") 
spec_m_re_cl <- subset(spec_re, CanopyHt=="m"&Site=="Cloquet") 
spec_b_nir_cl <- subset(spec_nir, CanopyHt=="b"&Site=="Cloquet") 
spec_m_nir_cl <- subset(spec_nir, CanopyHt=="m"&Site=="Cloquet") 
spec_b_swir_cl <- subset(spec_swir, CanopyHt=="b"&Site=="Cloquet") 
spec_m_swir_cl <- subset(spec_swir, CanopyHt=="m"&Site=="Cloquet") 
spec_b_biol_cl <- subset(spec_biol, CanopyHt=="b"&Site=="Cloquet") 
spec_m_biol_cl <- subset(spec_biol, CanopyHt=="m"&Site=="Cloquet") 

#_______________________________________________________________________________
# Calculating spectral angle mapper distance

#...All transmittance spectra (excluding FAB2)
dist_sam_all <- dist.speclib(spec_all_s, method = "sam") # dist object
# dist_sam_all <- sam_distance(spec_all_s) # matrix


#...Transmittance of parts of the spectrum 
dist_sam_uv <- dist.speclib(spec_uv, method = "sam") # dist object
dist_sam_vis <- dist.speclib(spec_vis, method = "sam") # dist object
dist_sam_re <- dist.speclib(spec_re, method = "sam") # dist object
dist_sam_nir <- dist.speclib(spec_nir, method = "sam") # dist object
dist_sam_swir <- dist.speclib(spec_swir, method = "sam") # dist object

#... Transmittance of parts of spectrum split by canopy height
dist_sam_uv_b <- dist.speclib(spec_b_uv, method = "sam")
dist_sam_vis_b <- dist.speclib(spec_b_vis, method = "sam")
dist_sam_re_b <- dist.speclib(spec_b_re, method = "sam")
dist_sam_nir_b <- dist.speclib(spec_b_nir, method = "sam")
dist_sam_swir_b <- dist.speclib(spec_b_swir, method = "sam")

dist_sam_uv_m <- dist.speclib(spec_m_uv, method = "sam")
dist_sam_vis_m <- dist.speclib(spec_m_vis, method = "sam")
dist_sam_re_m <- dist.speclib(spec_m_re, method = "sam")
dist_sam_nir_m <- dist.speclib(spec_m_nir, method = "sam")
dist_sam_swir_m <- dist.speclib(spec_m_swir, method = "sam")


#_______________________________________________________________________________
# Transmittance across whole spectrum and all canopy heights:
# NMDS:
mds_all <- metaMDS(dist_sam_all, autotransform=T, k=2, trymax=50)  
mds_all2 <- metaMDS(dist_sam_all, previous.best=mds_all, autotransform=F, k=2, trymax=50)

stressplot(mds_all2)
mds_all2$stress
plot(mds_all2, type="n")
points(mds_all2, pch=16, col=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv)))

points(mds_all2, pch=16, col=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv)))
symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc, inches=0.1,
        ann=F, bg=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv)), fg=NULL) # size of points varying by LAI


# wascores -- weighted averages scores for ordination configuration
all_points <- postMDS(mds_all2$points, dist_sam_all)
# all_wa <- wascores(all_points, specc) # doesn't like negative spectral values?
scores(all_points)
# wascores

# env: 2, 5, 11, 13, 17,18,19,20
fit_all <- envfit(mds_all2, SI(spec_all_s)[,c(2,5,11,13,17:20)], na.rm=T, perm=999) # fitting all predictors
# fit_sub <- envfit(mds_all2, SI(spec_all_s)[,c(2,5,12)], na.rm=T, perm=999) # fitting subset of predictors
# fit_sub2 <- envfit(mds_all2, SI(spec_all_s)[,c(2,5,14)], na.rm=T, perm=999) # fitting subset of predictors
print(fit_sub_vectors <- envfit(mds_all2, SI(spec_all_s)[,c(7,11,14)], na.rm=T, perm=999)) # fitting subset of predictors
print(fit_sub_factors <- envfit(mds_all2, SI(spec_all_s)[,c(2,5,12)], na.rm=T, perm=999)) # fitting subset of predictors

# plot(mds_all2)
# plot(mds_all2, type="n") 
# points(mds_all2, pch=16, col=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv)))
# symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc, inches=0.1, ## plotting point size by LAI
#         ann=F, fg=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv))+1, asp=1) ## plotting color by diversity
# legend("topleft",col=c(2,3,4,5), legend=c("Mono", "Congeners", "All angio or all gymno", "Mixtures"), bty="n", pch=1)
# color_palette1 <- colorRampPalette(c("#B71C1C", "#01579B", "#7E57C2"))(3)
symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc, inches=0.1, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_all_s)$AngioGymnoMix))+1, asp=1) ## plotting color by diversity
plot(fit_all, col="blue", cex=0.8, asp=1)
plot(fit_sub_vectors, col="black", cex=0.8, asp=1)
legend("topleft",col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n", pch=1)
legend("bottomright", legend=paste("Stress:", round(mds_all2$stress,3)), bty="n")
title("All sites - full spectrum")
# plot(mds_all2, type="n") 
# # points(mds_all2, pch=16, col=c("gray", "skyblue", "","","blue","darkblue", "","","","","","red")[SI(spec_all_s)$SR])
# # symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc, inches=0.1,
#         # ann=F, bg=as.numeric(as.factor(SI(spec_all_s)$MonoCongenCongpDiv)), fg=NULL)
# symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc, inches=0.1,
#         ann=F, bg=c("gray", "skyblue", "","","blue","darkblue", "","","","","","red")[SI(spec_all_s)$SR], fg=NULL, asp=1)
# legend("topleft",col=c("gray", "skyblue", "blue", "darkblue", "red"), legend=c("Mono", "Two spp", "Five spp", "Six spp", "12 spp"), bty="n", pch=16)
# legend("bottomright", legend=paste("Stress:", round(mds_all2$stress,3)), bty="n")
# plot(fit_sub_vectors, col="gray30", cex=0.5, asp=1)

# top vs bottom of canopy (show that these differ -- but not of much ecological interest)
# LAI -- note not available for all communities 
# community composition 

#_______________________________________________________________________________
#.... Transmittance different parts of the spectrum
dist_sam_uv <- dist.speclib(spec_uv, method = "sam") # dist object
dist_sam_vis <- dist.speclib(spec_vis, method = "sam") # dist object
dist_sam_re <- dist.speclib(spec_re, method = "sam") # dist object
dist_sam_nir <- dist.speclib(spec_nir, method = "sam") # dist object
dist_sam_swir <- dist.speclib(spec_swir, method = "sam") # dist object
dist_sam_biol <- dist.speclib(spec_biol, method = "sam") # dist object

# NMDS:
mds_uv <- metaMDS(dist_sam_uv, autotransform=T, k=2, trymax=50)  
mds_vis <- metaMDS(dist_sam_vis, autotransform=T, k=2, trymax=50)  
mds_re <- metaMDS(dist_sam_re, autotransform=T, k=2, trymax=50)  
mds_nir <- metaMDS(dist_sam_nir, autotransform=T, k=2, trymax=50)  
mds_swir <- metaMDS(dist_sam_swir, autotransform=T, k=2, trymax=50)  
mds_biol <- metaMDS(dist_sam_biol, autotransform=T, k=2, trymax=50)  

## 2, 5, 11, 13, 17, 18, 20
print(fit_uv <- envfit(mds_uv, SI(spec_uv)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors (2,5,7,11,12,14)
print(fit_uv_vectors <- envfit(mds_uv, SI(spec_uv)[,c(11,13,17,20)], na.rm=T, perm=999)) # LAI, SR, MonoCongenCongpDiv
print(fit_uv_factors <- envfit(mds_uv, SI(spec_uv)[,c(2,5,18)], na.rm=T, perm=999)) # Site, CanopyHt, AngioGymnoMix

print(fit_vis <- envfit(mds_vis, SI(spec_vis)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors
print(fit_vis_vectors <- envfit(mds_uv, SI(spec_uv)[,c(11,13,17,20)], na.rm=T, perm=999))
print(fit_vis_factors <- envfit(mds_uv, SI(spec_uv)[,c(2,5,18)], na.rm=T, perm=999))

print(fit_re <- envfit(mds_re, SI(spec_re)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors
print(fit_re_vectors <- envfit(mds_re, SI(spec_re)[,c(11,13,17,20)], na.rm=T, perm=999))
print(fit_re_factors <- envfit(mds_re, SI(spec_re)[,c(2,5,18)], na.rm=T, perm=999))

print(fit_nir <- envfit(mds_nir, SI(spec_nir)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors
print(fit_nir_vectors <- envfit(mds_nir, SI(spec_nir)[,c(11,13,17,20)], na.rm=T, perm=999))
print(fit_nir_factors <- envfit(mds_nir, SI(spec_nir)[,c(2,5,18)], na.rm=T, perm=999))

print(fit_swir <- envfit(mds_swir, SI(spec_swir)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors
print(fit_swir_vectors <- envfit(mds_swir, SI(spec_swir)[,c(11,13,17,20)], na.rm=T, perm=999))
print(fit_swir_factors <- envfit(mds_swir, SI(spec_swir)[,c(2,5,18)], na.rm=T, perm=999))

print(fit_biol <- envfit(mds_biol, SI(spec_biol)[,c(2,5,11,13,17,18,20)], na.rm=T, perm=999)) # fitting all predictors
print(fit_biol_vectors <- envfit(mds_biol, SI(spec_biol)[,c(11,13,17,20)], na.rm=T, perm=999))
print(fit_biol_factors <- envfit(mds_biol, SI(spec_biol)[,c(2,5,18)], na.rm=T, perm=999))


## Plotting all:
par(mfrow=c(3,2))
symbols(scores(mds_all2), circles=SI(spec_all_s)$LAI_allloc_corrected, inches=0.1, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_all_s)$AngioGymnoMix))+1, asp=1) ## plotting color by diversity
title("All sites - full spectrum")
legend("topleft",col=c(2,3,4), pch=1, legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_all2$stress,3)), bty="n")
plot(fit_sub_vectors, col="black", cex=0.8, asp=1)

symbols(scores(mds_uv), circles=SI(spec_uv)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_uv)$AngioGymnoMix))+1,  asp=1) ## plotting color by site
title("All sites - UV (350-400 nm)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_uv$stress,3)), bty="n")
plot(fit_uv_vectors, col="black", cex=0.8)

symbols(scores(mds_vis), circles=SI(spec_vis)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_vis)$AngioGymnoMix))+1, asp=1) ## plotting color by site
title("All sites - VIS (400-700 nm)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_vis$stress,3)), bty="n")
plot(fit_vis_vectors, col="black", cex=0.8)

symbols(scores(mds_re), circles=SI(spec_re)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_re)$AngioGymnoMix))+1, asp=1) ## plotting color by site
title("All sites - red edge (670-780 nm)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_re$stress,3)), bty="n")
plot(fit_re_vectors, col="black", cex=0.8)

symbols(scores(mds_nir), circles=SI(spec_nir)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_nir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
title("All sites - NIR (700-1100 nm)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_nir$stress,3)), bty="n")
plot(fit_nir_vectors, col="black", cex=0.8)

symbols(scores(mds_swir), circles=SI(spec_swir)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_swir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
title("All sites - SWIR (1100-2200 nm)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_swir$stress,3)), bty="n")
plot(fit_swir_vectors, col="black", cex=0.8)

symbols(scores(mds_biol), circles=SI(spec_biol)$LAI_allloc_corrected, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_biol)$AngioGymnoMix))+1, asp=1) ## plotting color by site
title("All sites - biol (360, 450, 660, 730)")
legend("topright", pch=1, col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_biol$stress,3)), bty="n")
plot(fit_biol_vectors, col="black", cex=0.8)

## Look at ordinations split by site? Or for a given canopy height to try to more clearly see diversity effects?


#_______________________________________________________________________________
#...Transmittance spectra below canopy
dist_sam_b <- dist.speclib(spec_b_s, method="sam") # dist object

# NMDS:
mds_b <- metaMDS(dist_sam_b, autotransform=T, k=2, trymax=50) 
mds_b2 <- metaMDS(dist_sam_b, previous.best=mds_b, autotransform=F, k=2, trymax=50)
mds_uv_b <- metaMDS(dist_sam_uv_b, autotransform=T, k=2, trymax=50) 
mds_vis_b <- metaMDS(dist_sam_vis_b, autotransform=T, k=2, trymax=50) 
mds_re_b <- metaMDS(dist_sam_re_b, autotransform=T, k=2, trymax=50) 
mds_nir_b <- metaMDS(dist_sam_nir_b, autotransform=T, k=2, trymax=50) 
mds_swir_b <- metaMDS(dist_sam_swir_b, autotransform=T, k=2, trymax=50) 

# Fitting vectors/factors
print(fit_b <- envfit(mds_b, SI(spec_b_s)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) # fitting all predictors
print(fit_uv_b <- envfit(mds_uv_b, SI(spec_b_uv)[,c(2,5,7,11,12,14)], na.rm=T, perm=999))
print(fit_vis_b <- envfit(mds_vis_b, SI(spec_b_vis)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 
print(fit_re_b <- envfit(mds_re_b, SI(spec_b_re)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 
print(fit_nir_b <- envfit(mds_nir_b, SI(spec_b_nir)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 
print(fit_swir_b <- envfit(mds_swir_b, SI(spec_b_swir)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 

print(fit_b_vectors <- envfit(mds_b, SI(spec_b_s)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_b_factors <- envfit(mds_b, SI(spec_b_s)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_uv_b_vectors <- envfit(mds_uv_b, SI(spec_b_uv)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_uv_b_factors <- envfit(mds_uv_b, SI(spec_b_uv)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_vis_b_vectors <- envfit(mds_vis_b, SI(spec_b_vis)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_vis_b_factors <- envfit(mds_vis_b, SI(spec_b_vis)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_re_b_vectors <- envfit(mds_re_b, SI(spec_b_re)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_re_b_factors <- envfit(mds_re_b, SI(spec_b_re)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_nir_b_vectors <- envfit(mds_nir_b, SI(spec_b_nir)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_nir_b_factors <- envfit(mds_nir_b, SI(spec_b_nir)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_swir_b_vectors <- envfit(mds_swir_b, SI(spec_b_swir)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_swir_b_factors <- envfit(mds_swir_b, SI(spec_b_swir)[,c(2,5,12)], na.rm=T, perm=999))


# # PCoA
# pcoa_b <- capscale(dist_sam_b ~ 1, data=dist_sam_b)
# screeplot(pcoa_b, type="line")
# summary(eigenvals(pcoa_b))
# plot(pcoa_b)

# Plotting
par(mfrow=c(3,2))
symbols(scores(mds_b2), circles=SI(spec_b_s)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_s)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n", pch=1)
legend("bottomright", legend=paste("Stress:", round(mds_b2$stress,3)), bty="n")
title("All sites - below canopy - full spectrum")
plot(fit_b_vectors, col="black", cex=0.5)

symbols(scores(mds_uv_b), circles=SI(spec_b_uv)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_uv)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_uv_b$stress,3)), bty="n")
title("All sites - below canopy - UV")
plot(fit_uv_b_vectors, col="black", cex=0.5)

symbols(scores(mds_vis_b), circles=SI(spec_b_vis)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_vis)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_vis_b$stress,3)), bty="n")
title("All sites - below canopy - VIS")
plot(fit_vis_b_vectors, col="black", cex=0.5)

symbols(scores(mds_re_b), circles=SI(spec_b_re)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_re)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_re_b$stress,3)), bty="n")
title("All sites - below canopy - RE")
plot(fit_re_b_vectors, col="black", cex=0.5)

symbols(scores(mds_nir_b), circles=SI(spec_b_nir)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_nir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_nir_b$stress,3)), bty="n")
title("All sites - below canopy - NIR")
plot(fit_nir_b_vectors, col="black", cex=0.5)

symbols(scores(mds_swir_b), circles=SI(spec_b_swir)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_b_swir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_swir_b$stress,3)), bty="n")
title("All sites - below canopy - SWIR")
plot(fit_swir_b_vectors, col="black", cex=0.5)


#_______________________________________________________________________________
#...Transmittance spectra mid canopy
dist_sam_m <- dist.speclib(spec_m_s, method="sam") # dist object

# NMDS:
mds_m <- metaMDS(dist_sam_m, autotransform=T, k=2, trymax=50) 
mds_m2 <- metaMDS(dist_sam_m, previous.best=mds_m, autotransform=F, k=2, trymax=50)
mds_uv_m <- metaMDS(dist_sam_uv_m, autotransform=T, k=2, trymax=50) 
mds_vis_m <- metaMDS(dist_sam_vis_m, autotransform=T, k=2, trymax=50) 
mds_re_m <- metaMDS(dist_sam_re_m, autotransform=T, k=2, trymax=50) 
mds_nir_m <- metaMDS(dist_sam_nir_m, autotransform=T, k=2, trymax=50) 
mds_swir_m <- metaMDS(dist_sam_swir_m, autotransform=T, k=2, trymax=50) 

# Fitting vectors/factors
print(fit_m <- envfit(mds_m, SI(spec_m_s)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) # fitting all predictors
print(fit_uv_m <- envfit(mds_uv_m, SI(spec_m_uv)[,c(2,5,7,11,12,14)], na.rm=T, perm=999))
print(fit_vis_m <- envfit(mds_vis_m, SI(spec_m_vis)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 
print(fit_re_m <- envfit(mds_re_m, SI(spec_m_re)[,c(2,5,7,11,12,14)], na.rm=T, perm=999))
print(fit_nir_m <- envfit(mds_nir_m, SI(spec_m_nir)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 
print(fit_swir_m <- envfit(mds_swir_m, SI(spec_m_swir)[,c(2,5,7,11,12,14)], na.rm=T, perm=999)) 

print(fit_m_vectors <- envfit(mds_m, SI(spec_m_s)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_m_factors <- envfit(mds_m, SI(spec_m_s)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_uv_m_vectors <- envfit(mds_uv_m, SI(spec_m_uv)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_uv_m_factors <- envfit(mds_uv_m, SI(spec_m_uv)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_vis_m_vectors <- envfit(mds_vis_m, SI(spec_m_vis)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_vis_m_factors <- envfit(mds_vis_m, SI(spec_m_vis)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_re_m_vectors <- envfit(mds_re_m, SI(spec_m_re)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_re_m_factors <- envfit(mds_re_m, SI(spec_m_re)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_nir_m_vectors <- envfit(mds_nir_m, SI(spec_m_nir)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_nir_m_factors <- envfit(mds_nir_m, SI(spec_m_nir)[,c(2,5,12)], na.rm=T, perm=999))

print(fit_swir_m_vectors <- envfit(mds_swir_m, SI(spec_m_swir)[,c(7,11,14)], na.rm=T, perm=999))
print(fit_swir_m_factors <- envfit(mds_swir_m, SI(spec_m_swir)[,c(2,5,12)], na.rm=T, perm=999))


# # PCoA
# pcoa_m <- capscale(dist_sam_m ~ 1, data=dist_sam_m)
# screeplot(pcoa_m, type="line")
# summary(eigenvals(pcoa_m))
# plot(pcoa_m)

# Plotting
par(mfrow=c(3,2))
symbols(scores(mds_m2), circles=SI(spec_m_s)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_s)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), bty="n", pch=1)
legend("bottomright", legend=paste("Stress:", round(mds_m2$stress,3)), bty="n")
title("All sites - mid canopy - full spectrum")
plot(fit_m_vectors, col="black", cex=0.5)

symbols(scores(mds_uv_m), circles=SI(spec_m_uv)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_uv)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_uv_m$stress,3)), bty="n")
title("All sites - mid canopy - UV")
plot(fit_uv_m_vectors, col="black", cex=0.5)

symbols(scores(mds_vis_m), circles=SI(spec_m_vis)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_vis)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_vis_m$stress,3)), bty="n")
title("All sites - mid canopy - VIS")
plot(fit_vis_m_vectors, col="black", cex=0.5)

symbols(scores(mds_re_m), circles=SI(spec_m_re)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_re)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_re_m$stress,3)), bty="n")
title("All sites - mid canopy - RE")
plot(fit_re_m_vectors, col="black", cex=0.5)

symbols(scores(mds_nir_m), circles=SI(spec_m_nir)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_nir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_nir_m$stress,3)), bty="n")
title("All sites - mid canopy - NIR")
plot(fit_nir_m_vectors, col="black", cex=0.5)

symbols(scores(mds_swir_m), circles=SI(spec_m_swir)$LAI_allloc, inches=0.08, ## plotting point size by LAI
        ann=F, fg=as.numeric(as.factor(SI(spec_m_swir)$AngioGymnoMix))+1, asp=1) ## plotting color by site
legend("topright", col=c(2,3,4), legend=c("Angio", "Gymno", "Mix"), pch=1, bty="n")
legend("bottomright", legend=paste("Stress:", round(mds_swir_m$stress,3)), bty="n")
title("All sites - mid canopy - SWIR")
plot(fit_swir_m_vectors, col="black", cex=0.5)




