################################################################################
##### Read files from spectrometers: ASD, PSM, and Spectra Vista
################################################################################

# Read and process .sig files based on a root path

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(asdreader)
library(lubridate)

#-------------------------------------------------------------------------------
# Load source

source("R/asd_files.R")
source("R/psm_files.R")
source("R/svc_files.R")
source("R/read_psm.R")
source("R/read_svc.R")

#-------------------------------------------------------------------------------
# Root file directory

path <- "E:/ligth_quality/IDENT-Cloquet"

################################################################################
#### ASD folders ####
################################################################################

#-------------------------------------------------------------------------------
# Step 1 - Organize .asd files based on root path

# Different in your machine
asd_folder_path <- paste0(path, "/ASD")

# Load files
frame_asd <- asd_files(asd_folder_path)

#-------------------------------------------------------------------------------
# Step 2 - Read files based on path and organize them

# Create paths
files_asd <- paste0(asd_folder_path, "/", frame_asd$folder, "/", frame_asd$file)

# Read spectra
asd_library <- as.data.table(get_spectra(f = files_asd,
                                         type = "radiance"))

asd_metadata <- do.call(rbind, lapply(files_asd, get_metadata))
asd_time <- as.vector(t(setDT(asd_metadata[,3])))

#Get time and date
datetime <- as_datetime(asd_time, tz = NULL)
date <- as.Date(datetime)
time <- format(datetime, format = "%H:%M:%S")

# Integrate meta and spectra to frame_asd
frame_asd$date <- date
frame_asd$time <- time
frame_asd <- cbind(frame_asd, asd_library)

#-------------------------------------------------------------------------------
# Step 3 - Export spectra
# NOTE: this files need to be converted to irradiance

fwrite(frame_asd, 
       paste0(asd_folder_path, "/", "asd_files.txt"),
       sep = "\t")



################################################################################
#### PSM folders ####
################################################################################

#-------------------------------------------------------------------------------
# Step 1 - Organize .psm files based on root path

# Different in your machine
psm_folder_path <- paste0(path, "/PSM")

# Load files
frame_psm <- psm_files(psm_folder_path)

#-------------------------------------------------------------------------------
# Step 2 - Read files based on path

# Create paths
files_psm <- paste0(psm_folder_path, "/", frame_psm$folder, "/", frame_psm$file)

# Read spectra
psm_library <- read_psm(path = files_psm,
                        column = 2) #1 Irrad. (Ref.), 2, Irrad. (Target), 3 Tgt./Ref.

# Integrate meta and spectra to frame_psm
frame_psm <- cbind(frame_psm, psm_library)

#-------------------------------------------------------------------------------
# Step 3 - Export spectra
# NOTE: some of these rows need to be filtered

fwrite(frame_psm, 
       paste0(psm_folder_path, "/", "psm_files.txt"),
       sep = "\t")


################################################################################
#### SVC folders ####
################################################################################

#-------------------------------------------------------------------------------
# Step 1 - Organize .sig files based on root path

# Different in your machine
svc_folder_path <- paste0(path, "/SVC")

# Load files
frame_svc <- svc_files(svc_folder_path)

#-------------------------------------------------------------------------------
# Step 2 - Read files based on path

# Create paths
files_svc <- paste0(svc_folder_path, "/", frame_svc$folder, "/", frame_svc$file)

# Read spectra
svc_library <- read_svc(paths = files_svc,
                        column = 2, #1 Irrad. (Ref.), 2, Irrad. (Target), 3 Tgt./Ref.
                        new_bands = 283:2500)

#Note on read spectra: .sig files have different bands, these are re sampled to a
# standard band span for storage.

# Integrate meta and spectra to frame_svc
frame_svc <- cbind(frame_svc, svc_library)

#-------------------------------------------------------------------------------
# Step 3 - Export spectra
# NOTE: some of these rows need to be filtered

fwrite(frame_svc, 
       paste0(svc_folder_path, "/", "svc_files.txt"),
       sep = "\t")
