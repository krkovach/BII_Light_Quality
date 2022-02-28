################################################################################
##### Read and Process Integrating Sphere files
################################################################################

# Read and process .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(spectrolab)

#-------------------------------------------------------------------------------
# Load source

source("R/sig_files.R")

#-------------------------------------------------------------------------------
# Step 1 - Organize .sig files based on root path

# Different in your machine
root_path <- "/media/antonio/antonio_ssd/ligth_quality/Integrating Sphere"

# Load files
frame_organized <- sig_files(root_path)

#-------------------------------------------------------------------------------
# Step 2 - Read files based on path

# Create paths
files_path <- paste0(root_path, "/", frame_organized$folder, "/", frame_organized$file)

# Read spectra
spectra_library <- read_spectra(path = files_path,
                                format = "sig",
                                type = "target_reflectance",
                                extract_metadata = TRUE)

metadata_spectra <- meta(spectra_library)

# Integrate meta to files
frame_organized$time <- metadata_spectra$time

#-------------------------------------------------------------------------------
# Step 3 - Match spectra at sensor transitions

spectra_match <- match_sensors(x = spectra_library, 
                               splice_at =  c(990, 1900), 
                               fixed_sensor = 2, 
                               interpolate_wvl = c(10, 2))

#-------------------------------------------------------------------------------
# Step 4 - Spectra library as frame for manual modifications

frame_spectra <- as.data.table(as.data.frame(spectra_match)[, 27:1019])

frame_organized <- cbind(frame_organized, frame_spectra)

#-------------------------------------------------------------------------------
# Step 5 - Export organized spectra

# Select your own path and name
export_pathname <- paste0(root_path, "/",
                          "spectra_organized.txt")

fwrite(frame_organized, export_pathname, sep = "\t")
# Please remember check for suspicions values
