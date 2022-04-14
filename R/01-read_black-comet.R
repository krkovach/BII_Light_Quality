################################################################################
##### Read and stack black comet files
################################################################################

# Read and process .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)

#-------------------------------------------------------------------------------
# Load source

source("R/ssm_files.R")
source("R/read_ssm.R")

#-------------------------------------------------------------------------------
# Step 1 - Organize .SSM files based on root path

# Different in your machine
ssm_folder_path <- "E:/ligth_quality/IDENT-Cloquet/BLK-C"

# Load files
frame_ssm <- ssm_files(ssm_folder_path)

#-------------------------------------------------------------------------------
# Step 2 - Read files based on path and organize them

# Create paths
files_ssm <- paste0(ssm_folder_path, "/", frame_ssm$folder, "/", frame_ssm$file)

# Read spectra
ssm_library <- read_ssm(path = files_ssm)

# Integrate meta and spectra to frame_ssm
frame_ssm <- cbind(frame_ssm, ssm_library)

#-------------------------------------------------------------------------------
# Step 3 - Export spectra

fwrite(frame_ssm, 
       paste0(ssm_folder_path, "/", "FAB_BLK-C.txt"),
       sep = "\t")
