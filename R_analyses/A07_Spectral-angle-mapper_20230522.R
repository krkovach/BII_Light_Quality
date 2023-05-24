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
all_spec <- fread(file.path(file_dir, "All-Sites_ASD-SVC_mean-transmittance_20230328.txt"))  # note this is mean transmittance, outlying individual transmittance spectra have been removed

#_______________________________________________________________________________
# Setting up spectral data as speclib
# (Spectral Angle Mapper is in hsdar package, which requires data in this format)

sub_all_spec <- all_spec[which(all_spec$Site!="FAB2"),] # excluding FAB2 (some plots have bad transmittance values)

wls <- as.numeric(c(350:1350,1450:1800, 1980:2200)) # wavelengths
specc <- as.matrix(sub_all_spec[,11:1583]) # spectral data
ids <- with(sub_all_spec, paste(Site, Plot, Mix, CanopyHt, sep="_")) # IDs for each spectrum
si <- data.frame(sub_all_spec[,1:10]) # other supplementary info
  
spec_all_s <- speclib(specc, wls) # setting up speclib

idSpeclib(spec_all_s) <- as.character(ids) # adding IDs for each spectrum
SI(spec_all_s) <- si # adding supplementary information

spec_b_s <- subset(spec_all_s, CanopyHt=="b") # subsetting below canopy spectra
spec_m_s <- subset(spec_all_s, CanopyHt=="m") # subsetting mid canopy spectra



#_______________________________________________________________________________
# Calculating spectral angle mapper distance

#...All transmittance spectra (excluding FAB2)
dist_sam_all <- dist.speclib(spec_all_s, method = "sam") # dist object
# dist_sam_all <- sam_distance(spec_all_s) # matrix
# NMDS:
mds_all <- metaMDS(dist_sam_all, autotransform=F, k=2, trymax=50) # no repeated best solution; perhaps try with three dimensions?
mds_all2 <- metaMDS(dist_sam_all, previous.best=mds_all, autotransform=F, k=2, trymax=50)

stressplot(mds_all2)
mds_all2$stress
plot(mds_all2, type="n")
points(mds_all2, pch=16, col=as.numeric(as.factor(SI(spec_all_s)$CanopyHt)))


#...Transmittance spectra below canopy
dist_sam_b <- dist.speclib(spec_b_s, method = "sam") # dist object

# NMDS:
mds_b <- metaMDS(dist_sam_b, autotransform=F, k=3, trymax=50) # no repeated best solution; perhaps try with three dimensions?
mds_b2 <- metaMDS(dist_sam_b, previous.best=mds_b, autotransform=F, k=3, trymax=50)
plot(mds_b2, type="n")
points(mds_b2, pch=16, col=as.numeric(as.factor(SI(spec_b_s)$Site)))
legend("bottomright", col=c(1,3,2), legend=c("Cloquet", "Freiburg", "FAB1"), pch=16, bty="n")

# PCoA
pcoa_b <- capscale(dist_sam_b ~ 1, data=dist_sam_b)
screeplot(pcoa_b, type="line")
summary(eigenvals(pcoa_b))
plot(pcoa_b)


#...Transmittance spectra mid canopy
dist_sam_m <- dist.speclib(spec_m_s, method = "sam") # dist object

# NMDS:
mds_m <- metaMDS(dist_sam_m, autotransform=F, k=2, trymax=50) # no repeated best solution; perhaps try with three dimensions?
mds_m2 <- metaMDS(dist_sam_m, previous.best=mds_m, autotransform=F, k=2, trymax=50)
plot(mds_m2, type="n")
points(mds_m2, pch=16, col=as.numeric(as.factor(SI(spec_m_s)$Site)))
legend("bottomright", col=c(1,3,2), legend=c("Cloquet", "Freiburg", "FAB1"), pch=16, bty="n")

# PCoA
pcoa_m <- capscale(dist_sam_m ~ 1, data=dist_sam_m)
screeplot(pcoa_m, type="line")
summary(eigenvals(pcoa_m))
plot(pcoa_m)



