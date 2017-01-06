# Prep data
# Load SRDB; filter for 'good' data (unmanipulated ecosystems, IRGA/GC only, etc);
# match with CRU climate data and Max Planck GPP dataset
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "1-prepdata.R"
PROBLEM       <- FALSE


# -----------------------------------------------------------------------------
# Extract CRU and Max Planck data given vectors of lon/lat/time info
extract_ncdf_data <- function(filename, lon, lat, midyear, nyears, file_startyear,
                              baseline = c(1961, 1990), print_every = 100) {
  
  assert_that(length(lon) == length(lat))
  assert_that(length(lon) == length(midyear))
  assert_that(length(midyear) == length(nyears))
  
  # Load file, decompressing first if necessary
  compressed <- grepl("gz$", filename)
  if(compressed) {
    printlog("Decompressing", filename)
    ncfile <- R.utils::gunzip(filename, remove = FALSE, overwrite = TRUE)
  } else {
    ncfile <- filename
  }
  library(raster) # 2.5.8
  nc <- brick(ncfile)
  
  # Results vectors: variable and variable normal (1961-1990 by default)  
  x <- rep(NA_real_, length(lon))
  normx <- rep(NA_real_, length(lon))
  varname <- attributes(nc@data)$zvar
  
  # Find nearest neighbors for all lon/lat pairs
  for(i in seq_along(lon)) {
    sp <- SpatialPoints(cbind(lon[i], lat[i]))
    midyear_layer <- (midyear[i] - file_startyear) * 12
    start_layer <- midyear_layer - nyears[i] / 2 * 12
    nlayers <- nyears[i] * 12
    
    printit <- print_every & i %% print_every == 0
    if(printit) {
      printlog("Extracting", i, lon[i], lat[i], midyear[i], nyears[i], "- layers", start_layer, 
               "to", start_layer + nlayers - 1, "...")
    }
    
    # Weirdly, raster::extract does not throw an error if we pass it a negative start layer
    # It just rolls merrily along, returning from the beginning of the file
    if(start_layer < 0 | start_layer + nlayers > nc@data@nlayers) {
      printlog(varname, "layer out of bounds #", i)
      next
    }
    
    # Extract the information for this point, over as many years as needed, then average
    nc %>%
      extract(sp, layer = start_layer, nl = nlayers) %>%
      mean(na.rm = TRUE) ->
      x[i]
    
    if(printit) {
      printlog(varname, "value:", x[i])
    }
    
    # Calculate baseline data
    start_layer <- max(1, (baseline[1] - file_startyear) * 12 + 1)
    nlayers <- (baseline[2] - baseline[1] + 1) * 12
    nc %>%
      extract(sp, layer = start_layer, nl = nlayers) %>%
      mean(na.rm = TRUE) ->
      normx[i]
    if(printit) {
      printlog(varname, "normal:", normx[i])
    }
  }
  
  # Clean up
  if(compressed) {
    printlog("Removing", ncfile)
    file.remove(ncfile)
  }
  
  out <- tibble(x = x, normx = normx)
  names(out) <- c(varname, paste0(varname, "_norm"))
  out
}


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# 1. Get SRDB data and filter

fn <- "/Users/d3x290/Documents/Work/Data-ongoing/Soil respiration database/srdb-data.xlsx"
printlog("Reading", fn)
srdb <- readxl::read_excel(fn)
print_dims(srdb)

printlog("Filtering...")
srdb %>%
  filter(!is.na(Longitude), !is.na(Latitude), 
         !is.na(Study_midyear), !is.na(YearsOfData),
         Ecosystem_state != "Managed", 
         Manipulation == "None",
         Meas_method %in% c("IRGA", "Gas chromatography")) %>%
  dplyr::select(Quality_flag, Study_midyear, YearsOfData, Longitude, Latitude, 
                Rs_annual, Rh_annual,
                GPP) ->
  d
print_dims(d)

#d <- data.frame(Longitude=-97, Latitude=55.5, Study_midyear=2001, YearsOfData=2)

# 2. Match with CRU climate data

fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.tmp.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/tmp/cru_ts3.24.1901.2015.tmp.dat.nc.gz
tmp <- extract_ncdf_data(fn, d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = 1901)
fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.pre.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/pre/cru_ts3.24.1901.2015.pre.dat.nc.gz
pre <- extract_ncdf_data(fn, d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = 1901)
fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.pet.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/pet/cru_ts3.24.1901.2015.pet.dat.nc.gz
pet <- extract_ncdf_data(fn, d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = 1901)

# 3. Match with GPP data
fn <- "/Users/d3x290/Data/MaxPlanck/201715151429EnsembleGPP_GL.nc.gz"
# Downloaded 5 Jan 2017 from https://www.bgc-jena.mpg.de/geodb/tmpdnld/201715151429EnsembleGPP_GL.nc
# See https://www.bgc-jena.mpg.de/bgi/index.php/Services/Overview
gpp <- extract_ncdf_data(fn, d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = 1982)
gpp <- gpp * 1000 * 60 * 60 * 24 * 365  # Convert from kgC/m2/s to gC/m2/yr

srdb_filtered <- bind_cols(d, tmp, pre, pet, gpp)

save_data(srdb_filtered)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
