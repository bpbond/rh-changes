# Prep data for heterotrophic respiration analysis
#
# Load SRDB; filter for 'good' data (unmanipulated ecosystems, IRGA/GC only, etc);
# spatially match with CRU climate, Max Planck GPP, and MODIS GPP datasets
#
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "1-prepdata.R"
PROBLEM       <- FALSE

library(raster) # 2.5.8



# -----------------------------------------------------------------------------
# For every SRDB record, find distance and ID of nearest FLUXNET tower 
# Slow but steady...
match_fluxnet <- function(d, fluxnet) {
  library(fossil)  # 0.3.7
  
  printlog("Starting FLUXNET nearest-neigbor matching...")
  y <- fluxnet[c("LOCATION_LONG", "LOCATION_LAT")]
  
  for(i in seq_len(nrow(d))) {
    if(i %% 10 == 0) {
      printlog(i, "of", nrow(d))
    }
    x <- d[i, c("Longitude", "Latitude")]
    names(y) <- names(x)
    z <- rbind(x, y)
    coordinates(z) <- ~ Longitude + Latitude
    dists <- fossil::earth.dist(z, dist = FALSE)[1,][-1]
    d$FLUXNET_DIST[i] <- dists[which.min(dists)]
    d$FLUXNET_SITE_ID[i] <- fluxnet$SITE_ID[which.min(dists)]
  }
  d  
}


# -----------------------------------------------------------------------------
# Extract data from a raster brick or raster stack, given vectors of lon/lat/time info
# This is general-purpose and called by both extract_ncdf_data and extract_geotiff_data below
extract_data <- function(rasterstack, varname, lon, lat, midyear, nyears, 
                         file_startyear, file_layers,
                         baseline = c(1961, 1990), print_every = 100) {
  
  printlog(SEPARATOR)
  printlog("Starting extraction for varname:", varname)
  
  # Results vectors: variable and variable normal (1961-1990 by default)  
  x <- rep(NA_real_, length(lon))
  normx <- rep(NA_real_, length(lon))
  
  # Find nearest neighbors for all lon/lat pairs
  for(i in seq_along(lon)) {
    sp <- SpatialPoints(cbind(lon[i], lat[i]))
    
    if(is.null(file_startyear)) {
      start_layer <- nlayers <- 1  # no time dimension
    } else {
      midyear_layer <- (midyear[i] - file_startyear) * 12
      start_layer <- midyear_layer - nyears[i] / 2 * 12
      nlayers <- nyears[i] * 12
    }
    
    printit <- print_every & i %% print_every == 0
    if(printit) {
      printlog("Extracting", i, lon[i], lat[i], midyear[i], nyears[i], "- layers", start_layer, 
               "to", start_layer + nlayers - 1, "...")
    }
    
    # Weirdly, raster::extract does not throw an error if we pass it a negative start layer
    # It just rolls merrily along, returning from the beginning of the file
    if(start_layer < 0 | start_layer + nlayers - 1 > file_layers) {
      printlog(varname, "layer out of bounds #", i)
      next
    }
    
    # Extract the information for this point, over as many years as needed, then average
    rasterstack %>%
      raster::extract(sp, layer = start_layer, nl = nlayers) %>%
      mean(na.rm = TRUE) ->
      x[i]
    
    if(printit) {
      printlog(varname, "value:", x[i])
    }
    
    # Calculate baseline (normal) data
    if(!is.null(baseline)) {
      start_layer <- max(1, (baseline[1] - file_startyear) * 12 + 1)
      nlayers <- (baseline[2] - baseline[1] + 1) * 12
      rasterstack %>%
        raster::extract(sp, layer = start_layer, nl = nlayers) %>%
        mean(na.rm = TRUE) ->
        normx[i]
      if(printit) {
        printlog(varname, "normal:", normx[i])
      }
    }
  }
  
  if(is.null(baseline)) {
    out <- tibble(x = x)
    names(out) <- c(varname)
  } else {
    out <- tibble(x = x, normx = normx)
    names(out) <- c(varname, paste0(varname, "_norm"))
  }
  out
}

# -----------------------------------------------------------------------------
# Extract MODIS NPP data given vectors of lon/lat/time info
extract_geotiff_data <- function(directory, varname, lon, lat, midyear, nyears, file_startyear,
                                 print_every = 100) {
  
  # Decompress if necessary
  zipfiles <- list.files(directory, pattern = "*.tif.gz$", full.names = TRUE)
  for(f in zipfiles) {
    printlog("Decompressing", basename(f))
    R.utils::gunzip(f, remove = FALSE, overwrite = TRUE)
  }
  
  files <- list.files(directory, pattern = "*.tif$", full.names = TRUE)
  printlog("Creating raster stack from", length(files), "files in", directory)
  
  # NB: extracting data from a a stack is **much** slower than from a brick(), as
  # used below. I don't anticipate running this program often enough to care, but FYI.
  nc <- stack(as.list(files))
  
  out <- extract_data(nc, varname, lon, lat, midyear, nyears, 
                      file_startyear = file_startyear, file_layers = length(files), 
                      baseline = NULL, print_every)
  
  # Clean up if we decompressed anything
  for(f in zipfiles) {
    printlog("Removing", f)
    file.remove(gsub(".gz$", "", f))  # remove the unzipped file
  }
  
  out
}

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
  nc <- brick(ncfile)
  varname <- attributes(nc@data)$zvar
  
  out <- extract_data(nc, varname, lon, lat, midyear, nyears, 
                      file_startyear = file_startyear, file_layers = nc@data@nlayers, 
                      baseline, print_every)
  
  # Clean up
  if(compressed) {
    printlog("Removing", ncfile)
    file.remove(ncfile)
  }
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
         is.na(Duplicate_record),
         Ecosystem_state != "Managed", 
         Manipulation == "None",
         Meas_method %in% c("IRGA", "Gas chromatography")) %>%
  dplyr::select(Quality_flag, Study_midyear, YearsOfData, Longitude, Latitude, 
                Rs_annual, Rh_annual, Stage,
                GPP) ->
  d
print_dims(d)

#d <- d[1,]


# 1.5. FLUXNET
fluxnet <- read_csv("outputs/fluxnet.csv")
d <- match_fluxnet(d, fluxnet)

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

# 3. Match with Max Planck GPP data
fn <- "/Users/d3x290/Data/MaxPlanck/201715151429EnsembleGPP_GL.nc.gz"
# Downloaded 5 Jan 2017 from https://www.bgc-jena.mpg.de/geodb/tmpdnld/201715151429EnsembleGPP_GL.nc
# See https://www.bgc-jena.mpg.de/bgi/index.php/Services/Overview
gpp <- extract_ncdf_data(fn, d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, baseline = NULL, file_startyear = 1982)
gpp <- gpp * 1000 * 60 * 60 * 24 * 365  # Convert from kgC/m2/s to gC/m2/yr

# 4. Match with MODIS GPP data
dir <- "/Users/d3x290/Data/MODIS_GPP/"
# Downloaded 6 Jan 2017 from http://www.ntsg.umt.edu/project/mod17
modisgpp <- extract_geotiff_data(dir, "modisgpp", d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = 2000)
modisgpp <- modisgpp * 0.1 # scale factor, per README file; results in gC/m2
modisgpp <- modisgpp * 12 # from mean monthly value to annual sum

# There are some crazy (>30,000 gC/m2) values in MODIS GPP. Remove those
modisgpp$modisgpp[modisgpp$modisgpp > 10000] <- NA

# 4. Match with SoilGrids1km data
# Downloaded 9 Jan 2017 from ftp://ftp.soilgrids.org/data/archive/12.Apr.2014/
dir <- "/Users/d3x290/Data/soilgrids1km/BLD/"
bd <- extract_geotiff_data(dir, "BD", d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = NULL)
dir <- "/Users/d3x290/Data/soilgrids1km/ORCDRC/"
orc <- extract_geotiff_data(dir, "ORC", d$Longitude, d$Latitude, d$Study_midyear, d$YearsOfData, file_startyear = NULL)

soc <- tibble(SOC = bd$BD * orc$ORC / 1000)  # kg C in top 1 m

# Done! Combine the various spatial data with the SRDB data and save
srdb_filtered <- bind_cols(d, tmp, pre, pet, gpp, modisgpp, soc)
save_data(srdb_filtered, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
